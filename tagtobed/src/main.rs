use std::env;
extern crate bam;
use crate::bam::RecordReader;
use crate::bam::record::tags::TagViewer;
use crate::bam::record::tags::TagValue;

fn main() {
    let mut record = bam::Record::new();
    let mut reader = bam::IndexedReader::from_path("test/test.bam").unwrap();
    let mut r  = reader.full();

    loop {
        match r.read_into(&mut record) {
            Ok(true) => {
                // retrieve tag and sequence
                match record.tags().get(b"MM") {
                    Some(bam::record::tags::TagValue::String(mm_tag_u8,  bam::record::tags::StringType::String))  => {
                        let mut mm_tag = std::str::from_utf8(mm_tag_u8).unwrap();
                        println!("{:?}", mm_tag);
                        // nothing to do if there are no base modifications
                        if mm_tag.len() == 0 {
                            break;
                        }
                        mm_tag = &mm_tag[0..mm_tag.len() - 1];

                        let seq = String::from_utf8(record.sequence().to_vec()).unwrap();
                        // parse mm tag and Cs in sequence
                        // assumes only one type of modification
                        let c_vec : Vec<usize> = seq.match_indices("C").map(|c| c.0).collect();
                        let mm_indices : Vec<i32> = mm_tag.split(',').skip(1).
                            map(|i| { i.parse::<i32>().unwrap() }).collect();

                        println!("{}", seq);
                        println!("{:?}", c_vec);
                        println!("{:?}", mm_indices);
                    }
                    _ => { panic!("No MM:Z tag."); }
                }

            },
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
    }
}
