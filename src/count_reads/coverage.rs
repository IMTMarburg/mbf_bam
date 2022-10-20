/*
 * Calculate coverage vectors over intervals
 *
 */
use crate::bam_ext::{open_bam};
use crate::rust_htslib::bam::Read;
use crate::BamError;
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions as htslib_record_extensions;
use std::path::{Path, PathBuf};

pub struct Interval<'a> {
    chr: &'a str,
    start: i64,
    stop: i64,
    flip: bool,
}

impl<'a> Interval<'a> {
    pub fn new(chr: &'a str, start: i64, stop: i64, flip: bool) -> Self {
        Interval {
            chr,
            start,
            stop,
            flip,
        }
    }
}

type Counts = Vec<u32>;

/// Calculate a coverage vector over a list of intervals
/// Intervals may be of unequal size.
/// one point per abse actually in a read
/// no extension etc.
pub fn calculate_coverage(
    bam_filename: impl AsRef<Path>,
    index_filename: Option<impl AsRef<Path>>,
    intervals: &[Interval],
    extend_reads: u32,
) -> Result<Vec<Counts>, BamError> {
    let pool = rayon::ThreadPoolBuilder::new().build().unwrap();
    let bam_filename = bam_filename.as_ref().to_owned();
    let index_filename = match index_filename {
        Some(x) => Some(x.as_ref().to_owned()),
        None => None,
    };
    let _bam = open_bam(&bam_filename, (&index_filename).as_ref())?; // just for error checking
    let mut chunked: Vec<_> = pool.install(|| {
        intervals
            .par_chunks(1000)
            .map(|chunk| {
                coverage_in_intervals(
                    chunk,
                    &bam_filename,
                    (&index_filename).as_ref(),
                    extend_reads,
                )
            })
            .collect()
    });
    let mut res = Vec::new();
    for el in chunked.iter_mut() {
        match el {
            Ok(el) => res.append(el),
            Err(e) => return Err(e.clone()),
        }
    }
    Ok(res)
}

fn coverage_in_intervals(
    chunk: &[Interval],
    bam_filename: &PathBuf,
    index_filename: Option<&PathBuf>,
    extend_reads: u32,
) -> Result<Vec<Counts>, BamError> {
    let mut bam = open_bam(bam_filename, index_filename).unwrap();
    let mut res: Vec<Vec<u32>> = Vec::new();
    for iv in chunk {
        let mut cov: Vec<u32> =
            coverage_in_interval(&mut bam, iv.chr, iv.start, iv.stop, extend_reads)?;
        if iv.flip {
            cov.reverse();
        }
        res.push(cov)
    }
    Ok(res)
}

fn coverage_in_interval(
    bam: &mut bam::IndexedReader,
    chr: &str,
    start: i64,
    stop: i64,
    extend_reads: u32,
) -> Result<Vec<u32>, BamError> {
    let mut res = vec![0; (stop - start) as usize];
    let mut read: bam::Record = bam::Record::new();
    let tid = bam
        .header()
        .tid(chr.as_bytes())
        .ok_or_else(|| BamError::UnknownError {
            msg: format!("Chromosome {} not found", &chr),
        })?;

    bam.fetch((tid, start, stop))?;

    while let Some(result) = bam.read(&mut read) {
        result?;

        let mut first: i64 = i64::MAX;
        let mut last: i64 = i64::MIN;
        for pos in read.reference_positions() {
            if (start <= pos) && (pos < stop) {
                first = first.min(pos);
                last = last.max(pos);
                res[((pos as i64) - start) as usize] += 1;
            }
        }
        if !read.is_reverse() {
            for pos in 0.max(first - start - (extend_reads as i64))..(first - start) {
                res[pos as usize] += 1;
            }
        } else {
            for pos in (stop - start)..(stop - start + (extend_reads as i64)).min(stop - start -1) {
                res[pos as usize] += 1;
            }
        }
    }
    Ok(res)
}
