use crate::BamError;
use bio::io::fastq;
use flate2::read::GzDecoder;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{index, record::Aux, Format, Header, Read, Reader, Record, Writer};
use std::collections::{HashMap, HashSet};
use ex::fs::File;
use std::io::BufReader;
use std::str;

/// substract all reads from the subtrahend bam file
/// from the minuend bam file,
/// writing to output bam file
/// and indexing it
pub fn py_substract_bam(
    output_filename: &str,
    minuend_filename: &str,
    subtrahend_filename: &str,
) -> Result<(), BamError> {
    let mut minuend = Reader::from_path(minuend_filename)?;
    let mut subtrahend = Reader::from_path(subtrahend_filename)?;
    let header = Header::from_template(minuend.header());

    let mut seen = HashSet::new();

    let mut read: Record = Record::new();
    while let Some(bam_result) = subtrahend.read(&mut read) {
        bam_result?;
        if !read.is_unmapped() {
            let q = read.qname().to_owned();
            seen.insert(q);
        }
    }

    {
        let mut output = Writer::from_path(output_filename, &header, Format::Bam)?;
        while let Some(bam_result) = minuend.read(&mut read) {
            bam_result?;
            if !seen.contains(read.qname()) {
                output.write(&read)?;
            }
        }
    } // output.drop will be called
    index::build(output_filename, None, index::Type::Bai, 4)?; //I see four threads
    Ok(())
}

fn read_gz_or_not(input_filename: &str) -> Result<Box<dyn std::io::Read>, BamError> {
    let file = File::open(input_filename)?;
    if input_filename.ends_with(".gz") {
        return Ok(Box::new(GzDecoder::new(file)));
    } else {
        return Ok(Box::new(file));
    }
}

pub fn py_annotate_barcodes_from_fastq(
    output_filename: &str,
    input_filename: &str,
    fastq2_filenames: Vec<&str>,
    barcodes: Vec<(String, usize, usize)>,
) -> Result<(), BamError> {
    let mut input = Reader::from_path(input_filename)?;
    let header = Header::from_template(input.header());

    let mut qname_to_tags: HashMap<Vec<u8>, Vec<(Vec<u8>, Vec<u8>)>> = HashMap::new();
    for filename in fastq2_filenames {
        let buf = BufReader::new(read_gz_or_not(filename)?);
        let reader = fastq::Reader::new(buf);
        for record in reader.records() {
            let record = record?;
            let mut tags: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
            for (tag, start, end) in barcodes.iter() {
                let barcode = (record.seq()[*start..*end]).to_vec();
                tags.push((tag.as_bytes().to_vec(), barcode));
            }
            qname_to_tags.insert(record.id().as_bytes().to_vec(), tags);
        }
    }
    {
        let mut output = Writer::from_path(output_filename, &header, Format::Bam)?;
        let mut read: Record = Record::new();
        while let Some(bam_result) = input.read(&mut read) {
            bam_result?;
            let tags = qname_to_tags
                .get(read.qname())
                .ok_or_else(|| BamError::UnknownError {
                    msg: format!(
                        "Read not found in fastq: {}",
                        std::str::from_utf8(read.qname()).unwrap()
                    )
                    .to_string(),
                })?;
            for (tag, value) in tags {
                read.push_aux(tag, Aux::String(std::str::from_utf8(value).unwrap())).unwrap();
            }
            output.write(&read)?;
        }
    } // output.drop will be called
    index::build(output_filename, None, index::Type::Bai, 4)?; //I see four threads
    Ok(())
}

pub fn bam_to_fastq(output_filename: &str, input_filename: &str) -> Result<(), BamError> {
    let mut input = Reader::from_path(input_filename)?;
    let mut output = fastq::Writer::to_file(output_filename)?;
    let mut read: Record = Record::new();
    while let Some(bam_result) = input.read(&mut read) {
        bam_result?;
        let q: Vec<u8> = match read.is_reverse() {
            true => read.qual().iter().map(|x| x + 33).rev().collect(),
            false => read.qual().iter().map(|x| x + 33).collect(),
        };
        let seq: Vec<u8> = match read.is_reverse() {
            true => bio::alphabets::dna::revcomp(read.seq().as_bytes()),
            false => read.seq().as_bytes(),
        };
        output.write(std::str::from_utf8(read.qname()).unwrap(), None, &seq, &q)?;
    }
    Ok(())
}


pub fn filter_and_rename_references(output_filename: &str, input_filename: &str, 
                                    reference_lookup: HashMap<String, Option<String>>) -> Result<(), BamError> {
    let mut input = Reader::from_path(input_filename)?;
    let old_header = Header::from_template(input.header());
    let mut new_header = Header::new();

    let mut tid_lookup = HashMap::new();
    //remove all SQ from header, keep the rest.
    let mut new_tid = 0;
    for (key, records) in old_header.to_hashmap() {
        dbg!(&key);
        if key != "SQ" {
            for record in records {
                let mut rec = HeaderRecord::new(key.as_bytes());
                for (k,v) in record.iter() {
                    rec.push_tag(k.as_bytes(),v);
                }
                new_header.push_record(&rec);
            }
        } else {
            let mut tid = 0;
            for record in records.iter() {
                let name = record.get("SN").ok_or_else(||BamError::UnknownError{msg:"SN in @SQ missing".to_string()})?;
                let length = record.get("LN")
                    .ok_or_else(||BamError::UnknownError{msg:"LN in @SQ missing".to_string()})?
                    .parse::<u32>().or_else(|_|Err(BamError::UnknownError{msg:"LN not an integere".to_string()}))?;
                let str_len = length.to_string();
                let new_name = reference_lookup.get(name).and_then(Option::as_ref);
                if let Some(new_name) = new_name {
                    let mut rec = HeaderRecord::new(b"SQ");
                    rec.push_tag(b"SN", new_name.to_string());
                    rec.push_tag(b"LN", (&str_len).to_string());
                    new_header.push_record(&rec);
                    tid_lookup.insert(tid, new_tid);
                    new_tid += 1;
                }
                tid += 1;
            }

        }
    }
    for comment in old_header.comments() {
        new_header.push_comment(comment.as_bytes());
    }

    let mut rec = HeaderRecord::new(b"PG");
    rec.push_tag(b"ID", "mbf_bam.filter_and_rename_references");
    new_header.push_record(&rec);

    if new_tid == 0 {
        return Err(BamError::UnknownError{msg:"No references left after filtering".to_string()});
    }
    //
    let mut output = Writer::from_path(output_filename, &new_header, Format::Bam)?;
    let mut read = Record::new();
    while let Some(bam_result) = input.read(&mut read) {
        bam_result?;
        let tid = read.tid();
        match tid_lookup.get(&tid) {
            Some(new_tid) => {
                read.set_tid(*new_tid);
                output.write(&read)?;
            }
            None => {}
        }
    } 
    Ok(())

}


#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use rust_htslib::bam::{Read, Reader, Record};

    fn count_reads(bam_filename: &str) -> Result<u32, crate::BamError> {
        let mut input = Reader::from_path(bam_filename)?;
        let mut count = 0;
        let mut read = Record::new();
        while let Some(bam_result) = input.read(&mut read) {
            bam_result?;
            count += 1;
        }
        Ok(count)
    }


    #[test]
    fn test_filter_and_rename(){
         let input = "sample_data/ex2.bam";
         let lookup: HashMap<_,_> = vec![("chr1".to_string(), Some("1".to_string())),
         ].into_iter().collect();

        let td = tempfile::Builder::new()
            .prefix("mbf_bam_test")
            .tempdir()
            .expect("could not create tempdir");
        let output = td.path().join("test.bam");
        let output_str = output.as_os_str().to_str().unwrap();
         crate::bam_manipulation::filter_and_rename_references(&output_str, input, lookup).unwrap();
         assert!(output.exists());
         assert!(count_reads(output_str).unwrap() == 1);
    }
    #[test]
    fn test_filter_and_rename2(){
         let input = "sample_data/ex2.bam";
         let lookup: HashMap<_,_> = vec![("chr2".to_string(), Some("2".to_string())),
            ("chr1".to_string(), None),
            ("chrX".to_string(), Some("X".to_string())), //not in ther.
         ].into_iter().collect();

        let td = tempfile::Builder::new()
            .prefix("mbf_bam_test")
            .tempdir()
            .expect("could not create tempdir");
        let output = td.path().join("test.bam");
        let output_str = output.as_os_str().to_str().unwrap();
         crate::bam_manipulation::filter_and_rename_references(&output_str, input, lookup).unwrap();
         assert!(output.exists());
         assert!(count_reads(output_str).unwrap() == 7);
    }

}
