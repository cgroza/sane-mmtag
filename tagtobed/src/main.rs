#![feature(split_array)]

extern crate bam;
use getopts::Options;
use std::io::Write;
use std::env;
use std::process;
use std::iter::zip;

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    opts.optopt("o", "", "set output file name", "NAME");
    opts.optopt("b", "", "set input bam", "NAME");
    opts.optopt("p", "", "probability threshold", "0-255");
    // opts.optopt("T", "", "tag name", "NAME");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
    };

    // write to stdout if no output provided
    let mut out_file : Box<dyn Write> = if !matches.opt_present("o") {
        Box::new(std::io::stdout())
    } else {
         Box::new(std::fs::File::create(matches.opt_str("o").unwrap()).unwrap())
    };

    // let tag_name = if matches.opt_present("T") {
    //     let tag = matches.opt_str("T").unwrap();
    //     *tag.as_bytes().split_array_ref::<2>().0
    // } else {
    //     *"MM".as_bytes().split_array_ref::<2>().0
    // };

    if !matches.opt_present("b") {
        eprintln!("No input provided.");
        process::exit(1);
    }

    if !matches.opt_present("p") {
        eprintln!("No probability threshold (ML) provided.");
        process::exit(1);
    }

    let p_threshold = matches.opt_str("p").unwrap().parse::<i64>().unwrap();

    // let mut record = bam::Record::new();
    let mut reader = bam::IndexedReader::from_path(matches.opt_str("b").unwrap()).unwrap();
    let r  = reader.full();

    for nr  in r {
        match nr {
            Ok(record) => {
                // retrieve tag and sequence
                match record.tags().get(b"MM").or_else(| | record.tags().get(b"Mm") ) {
                    Some(bam::record::tags::TagValue::String(mm_tag_u8,  bam::record::tags::StringType::String))  => {
                        match record.tags().get(b"ML").or_else(| | record.tags().get(b"Ml") ) {
                            Some(bam::record::tags::TagValue::IntArray(ml_tag))  => {
                        // skip secondary and supplementary alignments
                        if record.flag().is_supplementary() || record.flag().is_secondary() {
                            continue
                        }

                        let mut mm_tag = std::str::from_utf8(mm_tag_u8).unwrap();
                        // let ml_tag = std::str::from_utf8(ml_tag_array).unwrap();

                        // nothing to do if there are no base modifications
                        if mm_tag.len() == 0 {
                            continue;
                        }
                        mm_tag = &mm_tag[0..mm_tag.len() - 1];


                        // base modifications coordinates are on the original strand
                        let seq : String = if record.flag().is_reverse_strand() {
                            String::from_utf8(record.sequence().rev_compl(std::ops::RangeFull).collect()).unwrap()
                        } else {
                            String::from_utf8(record.sequence().to_vec()).unwrap()
                        };

                        // parse mm tag and Cs in sequence
                        // assumes only one type of modification
                        let c_vec : Vec<usize> = seq.match_indices("C").map(|c| c.0).collect();
                        let mm_indices : Vec<usize> = mm_tag.split(',').skip(1).
                                    map(|i| { i.parse::<usize>().unwrap() }).collect();


                        // let mm_probs : Vec<usize> = ml_tag.split(',').skip(1).
                        //             map(|i| { i.parse::<usize>().unwrap() }).collect();

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

                        let _ = out_file.write_fmt(format_args!("{name}\t{bases}\n",
                                 name = std::str::from_utf8(record.name()).unwrap(),
                                 bases = zip(mod_bases, ml_tag.iter()).filter(|pb| pb.1 >= p_threshold).map(|b| b.0.to_string()).collect::<Vec<String>>().join(",")));
                            }
                            _ => {} // for ml
                        }
                    }
                    _ => {}     // for mm
                }
            },
            Err(e) => {
                eprintln!("{}", e);
                process::exit(1);
            }
        }
    }
}
