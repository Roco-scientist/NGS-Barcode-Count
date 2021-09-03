use clap::{App, Arg};
use num_cpus;
use rayon;
// use rayon::prelude::*;
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
};
use time::PreciseTime;

fn main() {
    let start = PreciseTime::now();
    let (fastq, format, samples, threads) =
        arguments().unwrap_or_else(|err| panic!("Argument error: {}", err));

    let regex_string = del::del_info::regex_search(format).unwrap();

    let results = Arc::new(Mutex::new(HashMap::new()));

    let samples_hashmap = del::del_info::sample_barcodes(samples).unwrap();

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
            let results_clone = Arc::clone(&results);
            let samples_clone = samples_hashmap.clone();
            s.spawn(move |_| {
                del::parse_sequences::parse(
                    seq_clone,
                    finished_clone,
                    regex_string_clone,
                    results_clone,
                    samples_clone,
                    thread,
                )
                .unwrap();
            })
        }
    });
    let end = PreciseTime::now();
    let seconds = start.to(end).num_milliseconds() / 1000;
    println!("Total time: {} seconds", seconds);
    // println!("results: {:?}", results.lock().unwrap());
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
