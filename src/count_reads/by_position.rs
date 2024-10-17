
use super::chunked_genome::{Chunk, ChunkedGenome};
use crate::bam_ext::open_bam;
use crate::rust_htslib::bam::Read;
use crate::BamError;
use rayon::prelude::*;
use rust_htslib::bam;
use std::collections::HashMap;


pub type PositionCountResult = HashMap<String, u32>;

pub fn py_count_positions(
    filename: &str,
    index_filename: Option<&str>,
) -> Result<PositionCountResult, BamError> {
  //check whether the bam file can be openend
    //and we need it for the chunking
    let bam = open_bam(filename, index_filename)?;

    //perform the counting
    let cg = ChunkedGenome::new_without_tree(bam); // can't get the ParallelBridge to work with our lifetimes.
    let it: Vec<Chunk> = cg.iter().collect();
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    pool.install(|| {
        let result = it
            .into_par_iter()
            .map(|chunk| {
                let bam = open_bam(filename, index_filename).unwrap();

                let pos_count = count_positions(bam, chunk.tid, chunk.start, chunk.stop);
                let mut result: PositionCountResult = HashMap::new();
                match pos_count {
                    Ok(pos_count) => {
                        result.insert(chunk.chr.to_string(), pos_count);
                    }
                    Err(_) => (),
                };
                result
            })
            .reduce(|| PositionCountResult::new(), combine_position_results);
        //.fold(HashMap::<String, u32>::new(), add_hashmaps);
        Ok(result)
    })
}

fn count_positions(
    mut bam: bam::IndexedReader,
    tid: u32,
    start: u32,
    stop: u32,
) -> Result<u32, BamError> {
    bam.fetch((tid, start as u64, stop as u64))?;
    let mut read: bam::Record = bam::Record::new();
    let mut last_pos: i64 = -1;
    let mut pos_count = 0;
    while let Some(bam_result) = bam.read(&mut read) {
        bam_result?;
        // do not count multiple blocks matching in one gene multiple times
        if ((read.pos() as u32) < start) || ((read.pos() as u32) >= stop) {
            continue;
        }
        if read.pos() as i64 != last_pos {
            pos_count += 1;
            last_pos = read.pos();
        }
    }
    Ok(pos_count)
}



fn combine_position_results(mut a: PositionCountResult, b: PositionCountResult) -> PositionCountResult {
    for (key, count) in b.iter() {
        let target = &mut a
            .entry(key.to_string())
            .or_insert(0);
        **target += count;
    }
    a
}


