use clap::{App, Arg};
use rayon;
use rayon::prelude::*;
use std::fs::File;
use std::io::{self, BufRead};
// use std::path::Path;
use std::sync::{Arc, Mutex};

fn main() {
    let (fastq, _format) = arguments().unwrap_or_else(|err| panic!("Argument error: {}", err));
    rayon::scope(|s| {
        let seq = Arc::new(Mutex::new(None));
        let seq_clone = Arc::clone(&seq);

        let finished = Arc::new(Mutex::new(false));
        let finished_clone = Arc::clone(&finished);
        s.spawn(move |_| {
            let fastq_file =
                File::open(fastq).unwrap_or_else(|err| panic!("Error reading fastq file: {}", err));
            let mut line_num = 1;
            for line in io::BufReader::new(fastq_file).lines() {
                if line_num == 2 {
                    if let Ok(sequence_data) = line {
                        while seq_clone.lock().unwrap().is_some() {}
                        *seq_clone.lock().unwrap() = Some(sequence_data);
                    }
                }
                line_num += 1;
                if line_num == 5 {
                    line_num = 1
                }
            }
            println!("finished");
            *finished_clone.lock().unwrap() = true;
        });

        let seq_clone = Arc::clone(&seq);
        let finished_clone = Arc::clone(&finished);
        s.spawn(move |_| loop {
            while seq_clone.lock().unwrap().is_none() {
                if *finished_clone.lock().unwrap() {
                    break;
                }
            }
            if *finished_clone.lock().unwrap() {
                break;
            }
            println!("Sequence: {}", seq_clone.lock().unwrap().as_ref().unwrap());
            *seq_clone.lock().unwrap() = None;
        })
    });
}

/// Gets the command line arguments
pub fn arguments() -> Result<(String, String), Box<dyn std::error::Error>> {
    // parse arguments
    let args = App::new("DEL analysis")
        .version("0.1")
        .author("Rory Coffey <coffeyrt@gmail.com>")
        .about("Counts DEL hits from fastq files")
        .arg(
            Arg::with_name("fastq")
                .short("f")
                .long("fastq")
                .takes_value(true)
                .required(true)
                .help("FASTQ file unzipped"),
        )
        .arg(
            Arg::with_name("sequence_format")
                .short("s")
                .long("sequence_format")
                .takes_value(true)
                .help("Sequence format file"),
        )
        .get_matches();

    return Ok((
        args.value_of("fastq").unwrap().to_string(),
        args.value_of("sequence_format")
            .unwrap_or("nothing")
            .to_string(),
    ));
}
