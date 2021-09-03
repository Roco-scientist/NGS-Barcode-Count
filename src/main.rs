use clap::{App, Arg};
use num_cpus;
use rayon;
// use rayon::prelude::*;
use std::sync::{Arc, Mutex};

fn main() {
    let (fastq, format, _samples, threads) =
        arguments().unwrap_or_else(|err| panic!("Argument error: {}", err));

    let regex_string = del::del_info::regex_search(format).unwrap();

    rayon::scope(|s| {
        let seq = Arc::new(Mutex::new(Vec::new()));
        let finished = Arc::new(Mutex::new(false));

        let seq_clone = Arc::clone(&seq);
        let finished_clone = Arc::clone(&finished);
        s.spawn(move |_| {
            del::read_fastq(fastq, seq_clone).unwrap();
            println!("finished");
            *finished_clone.lock().unwrap() = true;
        });

        for thread in 1..threads {
            let seq_clone = Arc::clone(&seq);
            let finished_clone = Arc::clone(&finished);
            let regex_string_clone = regex_string.clone();
            s.spawn(move |_| {
                del::parse_sequences::parse(seq_clone, finished_clone, regex_string_clone, thread)
                    .unwrap();
            })
        }
    });
}

/// Gets the command line arguments
pub fn arguments() -> Result<(String, String, String, u8), Box<dyn std::error::Error>> {
    let total_cpus = num_cpus::get().to_string();
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
                .required(true)
                .help("Sequence format file"),
        )
        .arg(
            Arg::with_name("sample_barcodes")
                .short("b")
                .long("sample_barcodes")
                .takes_value(true)
                .required(true)
                .help("Sample barcodes file"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .default_value(&total_cpus)
                .help("Number of threads"),
        )
        .get_matches();

    return Ok((
        args.value_of("fastq").unwrap().to_string(),
        args.value_of("sequence_format").unwrap().to_string(),
        args.value_of("sample_barcodes").unwrap().to_string(),
        args.value_of("threads").unwrap().parse::<u8>().unwrap(),
    ));
}
