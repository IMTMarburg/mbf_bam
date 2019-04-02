/// Read counting functions
/// 
/// These are smart, fast read counters
/// that avoid the more common pitfalls of htseq data
/// in terms of handling reads matching multiple genes
/// or being multi mapped
///

use bio::data_structures::interval_tree::IntervalTree;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::prelude::*;
use pyo3::prelude::*;

use pyo3::types::{PyList, PyObjectRef, PyTuple};
use crate::bam_ext::{BamRecordExtensions, open_bam};
use crate::BamError;

use std::collections::{HashMap, HashSet};
fn add_hashmaps(mut a: HashMap<String, u32>, b: HashMap<String, u32>) -> HashMap<String, u32> {
    for (k, v) in b.iter() {
        let x = a.entry(k.to_string()).or_insert(0);
        *x += v;
    }
    a
}

/// OurTree stores an interval tree 
/// a gene_no (ie. an index into a vector of gene_ids)
/// and the strand (+1/ -1)
pub type OurTree = IntervalTree<u32, (u32, i8)>;

/// build_tree converts an list of python intervals
/// into an OurTree and a Vec of gene_ids
///
/// python intervals are a list of tuples
/// (gene_stable_id, strand, start, stop)
/// - each reference sequence has it's own list.
pub fn build_tree(iv_obj: &PyObjectRef) -> Result<(OurTree, Vec<String>), PyErr> {
    let iv_list: &PyList = iv_obj.extract()?;
    let mut tree = IntervalTree::new();
    let mut gene_ids = Vec::new();
    for (gene_no, iv_entry_obj) in iv_list.iter().enumerate() {
        let iv_tuple: &PyTuple = iv_entry_obj.extract()?;
        let lgene_id: String = iv_tuple.get_item(0).extract()?;
        gene_ids.push(lgene_id);
        let lstrand: i8 = iv_tuple.get_item(1).extract()?;
        let lstart: &PyList = iv_tuple.get_item(2).extract()?;
        let lend: &PyList = iv_tuple.get_item(3).extract()?;
        for (ls, le) in lstart.iter().zip(lend.iter()) {
            let ls: u32 = ls.extract()?;
            let le: u32 = le.extract()?;
            tree.insert(ls..le, (gene_no as u32, lstrand))
        }
    }
    Ok((tree, gene_ids))
}

///count_reads_in_region_unstranded
///
///counts the unstranded reads in a region, 
///matching the to the tree entries.
fn count_reads_in_region_unstranded(
    mut bam: bam::IndexedReader,
    tree: &OurTree,
    tid: u32,
    start: u32,
    stop: u32,
    gene_count: u32,
) -> Result<Vec<u32>, BamError> {
    let mut result = vec![0; gene_count as usize];
    let mut multimapper_dedup: HashMap<u32, HashSet<Vec<u8>>> = HashMap::new();
    let mut gene_nos_seen = HashSet::<u32>::new();
    let mut read: bam::Record = bam::Record::new();
    bam.fetch(tid, start, stop)?;
    while let Ok(_) = bam.read(&mut read) {
        let blocks = read.blocks();
        // do not count multiple blocks matching in one gene multiple times
        gene_nos_seen.clear();
        for iv in blocks.iter() {
            if (iv.1 < start) || iv.0 >= stop || ((iv.0 < start) && (iv.1 >= start)) {
                // if this block is outside of the region
                // don't count it at all.
                // if it is on a block boundary
                // only count it for the left side.
                // which is ok, since we place the blocks to the right
                // of our intervals.
                continue;
            }
            for r in tree.find(iv.0..iv.1) {
                let entry = r.data();
                let gene_no = (*entry).0;
                let nh = read.aux(b"NH");
                let nh = nh.map_or(1, |aux| aux.integer());
                if nh == 1 {
                    gene_nos_seen.insert(gene_no);
                } else {
                    let hs = multimapper_dedup
                        .entry(gene_no)
                        .or_insert_with(HashSet::new);
                    hs.insert(read.qname().to_vec());
                }
                /*if gene_ids[gene_no as usize] == "FBgn0037275" {
                println!(
                    "{}, {}, {}",
                    start,
                    stop,
                    std::str::from_utf8(read.qname()).unwrap()
                );
                }*/
            }
        }
        for gene_no in gene_nos_seen.iter() {
            result[*gene_no as usize] += 1;
        }
    }
    for (gene_no, hs) in multimapper_dedup.iter() {
        result[*gene_no as usize] += hs.len() as u32;
    }
    Ok(result)
}

struct ChunkedGenome {
    trees: HashMap<String, (OurTree, Vec<String>)>,
    bam: bam::IndexedReader,
}

impl ChunkedGenome {
    fn new(
        trees: HashMap<String, (OurTree, Vec<String>)>,
        bam: bam::IndexedReader,
    ) -> ChunkedGenome {
        ChunkedGenome { trees, bam }
    }

    fn iter(&self) -> ChunkedGenomeIterator {
        ChunkedGenomeIterator {
            cg: &self,
            it: self.trees.keys(),
            last_start: 0,
            last_tid: 0,
            last_chr_length: 0,
            last_chr: "".to_string(),
        }
    }
}

struct ChunkedGenomeIterator<'a> {
    cg: &'a ChunkedGenome,
    it: std::collections::hash_map::Keys<'a, String, (OurTree, Vec<String>)>,
    last_start: u32,
    last_chr: String,
    last_tid: u32,
    last_chr_length: u32,
}
struct Chunk<'a> {
    chr: String,
    tid: u32,
    tree: &'a OurTree,
    gene_ids: &'a Vec<String>,
    start: u32,
    stop: u32,
}

impl<'a> Iterator for ChunkedGenomeIterator<'a> {
    type Item = Chunk<'a>;
    fn next(&mut self) -> Option<Chunk<'a>> {
        let chunk_size = 1_000_000;
        if self.last_start >= self.last_chr_length {
            let next_chr = match self.it.next() {
                Some(x) => x,
                None => return None,
            };
            let tid = self.cg.bam.header().tid(next_chr.as_bytes()).unwrap();
            let chr_length = self.cg.bam.header().target_len(tid).unwrap();
            self.last_tid = tid;
            self.last_chr_length = chr_length;
            self.last_chr = next_chr.to_string();
            self.last_start = 0;
        }

        let (next_tree, next_gene_ids) = self.cg.trees.get(&self.last_chr).unwrap();
        let mut stop = self.last_start + chunk_size;
        loop {
            //TODO: migh have to extend this for exon based counting to not
            //cut gene in half?
            //option 0 for that is to pass in the gene intervals as well
            //just for constructing the chunks
            //option 1 is to get the immediate left/right entrys (from the tree?)
            //and if they have the same gene_no -> advance...
            //right is easy, just find (stop..length) and take only the next()
            //left is more difficult.
            let overlapping = next_tree.find(stop..stop+1).next();
            match overlapping {
                None => break,
                Some(entry) => {
                    let iv = entry.interval();
                    stop = iv.end + 1;
                }
            }
        }
        let c = Chunk {
            chr: self.last_chr.clone(),
            tid: self.last_tid,
            tree: next_tree,
            gene_ids: next_gene_ids,
            start: self.last_start,
            stop,
        };
        self.last_start = stop;
        Some(c)
    }
}

/// python wrapper for py_count_reads_unstranded
pub fn py_count_reads_unstranded(
    filename: &str,
    index_filename: Option<&str>,
    trees: HashMap<String, (OurTree, Vec<String>)>,
) -> Result<HashMap<String, u32>, BamError> {
    //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;
    
    //perform the counting
    let cg = ChunkedGenome::new(trees, bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let result = it
        .into_par_iter()
        .map(|chunk| {
          let bam = open_bam(filename, index_filename).unwrap();
            
            let counts = count_reads_in_region_unstranded(
                bam,
                &chunk.tree,
                chunk.tid,
                chunk.start,
                chunk.stop,
                chunk.gene_ids.len() as u32,
            );
            let mut total = 0;
            let mut result: HashMap<String, u32> = match counts {
                Ok(counts) => {
                    let mut res = HashMap::new();
                    for (gene_no, cnt) in counts.iter().enumerate() {
                        let gene_id = &chunk.gene_ids[gene_no];
                        res.insert(gene_id.to_string(), *cnt);
                        total += cnt;
                    }
                    res
                }
                _ => HashMap::new(),
            };
            result.insert("_total".to_string(), total);
            result.insert(format!("_{}", chunk.chr), total);
            result
        })
        .reduce(HashMap::<String, u32>::new, add_hashmaps);
    //.fold(HashMap::<String, u32>::new(), add_hashmaps);
    Ok(result)
}

