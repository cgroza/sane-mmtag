
extern crate bam;
use crate::bam::RecordReader;

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
                        // println!("{:?}", mm_tag);
                        // nothing to do if there are no base modifications
                        if mm_tag.len() == 0 {
                            break;
                        }
                        mm_tag = &mm_tag[0..mm_tag.len() - 1];

                        let seq : String;
                        // base modifications coordinates are on the original strand
                        if record.flag().is_reverse_strand() {
                            seq = String::from_utf8(record.sequence().rev_compl(std::ops::RangeFull).collect()).unwrap();
                        } else {
                            seq = String::from_utf8(record.sequence().to_vec()).unwrap();
                        }

                        // parse mm tag and Cs in sequence
                        // assumes only one type of modification
                        let c_vec : Vec<usize> = seq.match_indices("C").map(|c| c.0).collect();
                        let mm_indices : Vec<usize> = mm_tag.split(',').skip(1).
                            map(|i| { i.parse::<usize>().unwrap() }).collect();
                        let max_c_vec_index : usize = mm_indices.iter().map(|a| { if *a == 0 { 1 } else { *a } } ).sum();

                        // make sure the base modification coordinates do not go past the last potentially modified bases
                        assert!(max_c_vec_index < c_vec.len());

                        // calculate position of modified bases
                        let mut mod_bases : Vec<usize> = Vec::new();
                        let mut i = 0;
                        for mc in &mm_indices {
                            // the current base is modified
                            if *mc == 0 {
                                mod_bases.push(c_vec[i]);
                            // skip mc bases, then take the modified base
                            } else {
                                i = i + mc;
                                mod_bases.push(c_vec[i]);
                            }
                            // we consumed the modified base
                            i = i + 1;
                        }

                        println!("{:?}", c_vec);
                        println!("{:?}", mm_indices);
                        println!("{:?}", mod_bases);
                        println!("----------------");
                    }
                    _ => { panic!("No MM:Z tag!"); }
                }

            },
            Ok(false) => break,
            Err(e) => panic!("{}", e),
        }
    }
}
