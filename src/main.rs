use chrono::Local;
use clap::{App, Arg};
use num_cpus;
use rayon;
// use rayon::prelude::*;
use std::{
    collections::HashMap,
    error::Error,
    sync::{Arc, Mutex},
    time::Instant,
};

fn main() {
    // Start a clock to measure how long the algorithm takes
    let start = Instant::now();

    // get the argument inputs
    let args = Args::new().unwrap_or_else(|err| panic!("Argument error: {}", err));

    let sequence_format = del::del_info::SequenceFormat::new(args.format.clone())
        .unwrap_or_else(|err| panic!("sequence format error: {}", err));
    sequence_format.display_format();

    // Create a results hashmap that will contain the counts.  This is passed between threads
    let results = Arc::new(Mutex::new(HashMap::new()));
    // Create a random_barcodes hashmap to keep track of the random barcodes.  This way if it already was found for the same sample and building blocks
    // it will not be counted
    let random_barcodes = Arc::new(Mutex::new(HashMap::new()));

    // Create a hashmap of the sample barcodes in order to convert sequence to sample ID
    let samples_hashmap_option;
    if let Some(ref samples) = args.sample_barcodes_option {
        let samples_hashmap = del::del_info::sample_barcode_file_conversion(samples).unwrap();
        samples_hashmap_option = Some(samples_hashmap);
    } else {
        samples_hashmap_option = None
    }

    // Create a hashmap of the building block barcodes in order to convert sequence to building block
    let bb_hashmap;
    if let Some(ref bb) = args.bb_barcodes_option {
        bb_hashmap =
            Some(del::del_info::bb_barcode_file_conversion(bb, sequence_format.bb_num).unwrap());
    } else {
        bb_hashmap = None
    }

    // Create a sequencing errors Struct to track errors.  This is passed between threads
    let sequence_errors = Arc::new(Mutex::new(del::del_info::SequenceErrors::new()));

    // Create a passed exit passed variable to stop reading when a thread has panicked
    let exit = Arc::new(Mutex::new(false));

    // Create a MaxSeqErrors struct which holds how many sequencing errors are allowed for each sequencing region
    let mut max_errors = del::del_info::MaxSeqErrors::new(
        args.sample_errors_option,
        args.bb_errors_option,
        args.constant_errors_option,
        &sequence_format.regex_string,
        &sequence_format.format_string,
    )
    .unwrap_or_else(|err| panic!("Max Sequencing Errors error: {}", err));
    // Display region sizes and errors allowed
    max_errors.display();

    // Start the multithreading scope
    rayon::scope(|s| {
        // Create a sequence vec which will have sequences entered by the reading thread, and sequences removed by the processing threads
        let seq = Arc::new(Mutex::new(Vec::new()));
        // Create a passed variable to let the processing threads know the reading thread is done
        let finished = Arc::new(Mutex::new(false));

        // Clone variables that are needed to be passed into the reading thread and create the reading thread
        let seq_clone = Arc::clone(&seq);
        let finished_clone = Arc::clone(&finished);
        let exit_clone = Arc::clone(&exit);
        let fastq = args.fastq.clone();
        s.spawn(move |_| {
            del::read_fastq(fastq, seq_clone, exit_clone).unwrap();
            *finished_clone.lock().unwrap() = true;
        });

        // Create processing threads.  One less than the total threads because of the single reading thread
        for _ in 1..args.threads {
            // Clone all variables needed to pass into each thread
            let seq_clone = Arc::clone(&seq);
            let finished_clone = Arc::clone(&finished);
            let sequence_format_clone = sequence_format.clone();
            let results_clone = Arc::clone(&results);
            let random_barcodes_clone = Arc::clone(&random_barcodes);
            let samples_clone = samples_hashmap_option.clone();
            let bb_clone = bb_hashmap.clone();
            let sequence_errors_clone = Arc::clone(&sequence_errors);
            let exit_clone = Arc::clone(&exit);
            let max_errors_clone = max_errors.clone();

            // Create a processing thread
            s.spawn(move |_| {
                del::parse_sequences::parse(
                    seq_clone,
                    finished_clone,
                    sequence_format_clone,
                    results_clone,
                    random_barcodes_clone,
                    samples_clone,
                    bb_clone,
                    sequence_errors_clone,
                    max_errors_clone,
                )
                .unwrap_or_else(|err| {
                    *exit_clone.lock().unwrap() = true;
                    panic!("Compute thread panic error: {}", err)
                });
            })
        }
    });

    // Print sequencing error counts to stdout
    sequence_errors.lock().unwrap().display();

    println!();
    // Get the end time and print compute time for the algorithm
    let elapsed_time = start.elapsed();
    if elapsed_time.as_secs() < 3 {
        println!("Compute time: {} milliseconds", elapsed_time.as_millis());
    } else {
        if elapsed_time.as_secs() > 600 {
            println!("Compute time: {} minutes", elapsed_time.as_secs() / 60)
        } else {
            println!("Compute time: {} seconds", elapsed_time.as_secs())
        }
    }

    println!();
    println!("Writing counts");
    del::output_counts(
        args.output_dir,
        results,
        sequence_format,
        bb_hashmap,
        args.prefix,
        args.merge_output,
    )
    .unwrap();
    // Get the end time and print total time for the algorithm
    let elapsed_time = start.elapsed();
    if elapsed_time.as_secs() < 3 {
        println!("Total time: {} milliseconds", elapsed_time.as_millis());
    } else {
        if elapsed_time.as_secs() > 600 {
            println!("Total time: {} minutes", elapsed_time.as_secs() / 60)
        } else {
            println!("Total time: {} seconds", elapsed_time.as_secs())
        }
    }
}

/// A struct that contains and initiates all input arguments
struct Args {
    fastq: String,                          // fastq file path
    format: String,                         // format scheme file path
    sample_barcodes_option: Option<String>, // sample barcode file path.  Optional
    bb_barcodes_option: Option<String>,     // building block barcode file path. Optional
    output_dir: String,                     // output directory.  Deafaults to './'
    threads: u8, // Number of threads to use.  Defaults to number of threads on the machine
    prefix: String, // Prefix string for the output files
    merge_output: bool, // Whether or not to create an additional output file that merges all samples
    bb_errors_option: Option<usize>, // Optional input of how many errors are allowed in each building block barcode.  Defaults to 20% of the length
    sample_errors_option: Option<usize>, // Optional input of how many errors are allowed in each sample barcode.  Defaults to 20% of the length
    constant_errors_option: Option<usize>, // Optional input of how many errors are allowed in each constant region barcode.  Defaults to 20% of the length
}

impl Args {
    fn new() -> Result<Args, Box<dyn Error>> {
        let total_cpus = num_cpus::get().to_string();
        let today = Local::today().format("%Y-%m-%d").to_string();
        // parse arguments
        let args = App::new("DEL analysis")
        .version("0.4.1")
        .author("Rory Coffey <coffeyrt@gmail.com>")
        .about("Counts DEL hits from fastq files and optional does conversions of sample IDs and building block IDs")
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
                .short("q")
                .long("sequence_format")
                .takes_value(true)
                .required(true)
                .help("Sequence format file"),
        )
        .arg(
            Arg::with_name("sample_barcodes")
                .short("s")
                .long("sample_barcodes")
                .takes_value(true)
                .help("Sample barcodes file"),
        )
        .arg(
            Arg::with_name("bb_barcodes")
                .short("b")
                .long("bb_barcodes")
                .takes_value(true)
                .help("Building block barcodes file"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .default_value(&total_cpus)
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("output_dir")
                .short("o")
                .long("output_dir")
                .takes_value(true)
                .default_value("./")
                .help("Directory to output the counts to"),
        )
        .arg(
            Arg::with_name("prefix")
                .short("p")
                .long("prefix")
                .takes_value(true)
                .default_value(&today)
                .help("File prefix name.  THe output will end with '_<sample_name>_counts.csv'"),
        )
        .arg(
            Arg::with_name("merge_output")
                .short("m")
                .long("merge_output")
                .takes_value(false)
                .help("Merge sample output counts into a single file.  Not necessary when there is only one sample"),
        )
        .arg(
            Arg::with_name("bb_errors")
                .long("bb_errors")
                .takes_value(true)
                .help("Maximimum number of sequence errors allowed within each building block barcode. Defaults to 20% of the total."),
        )
        .arg(
            Arg::with_name("sample_errors")
                .long("sample_errors")
                .takes_value(true)
                .help("Maximimum number of sequence errors allowed within sample barcode. Defaults to 20% of the total."),
        )
        .arg(
            Arg::with_name("contant_errors")
                .long("constant_errors")
                .takes_value(true)
                .help("Maximimum number of sequence errors allowed within constant region. Defaults to 20% of the total."),
        )
        .get_matches();

        let sample_barcodes_option;
        if let Some(sample) = args.value_of("sample_barcodes") {
            sample_barcodes_option = Some(sample.to_string())
        } else {
            sample_barcodes_option = None
        }

        let bb_barcodes_option;
        if let Some(bb) = args.value_of("bb_barcodes") {
            bb_barcodes_option = Some(bb.to_string())
        } else {
            bb_barcodes_option = None
        }

        let bb_errors_option;
        if let Some(bb) = args.value_of("bb_errors") {
            bb_errors_option = Some(bb.parse::<usize>()?)
        } else {
            bb_errors_option = None
        }

        let sample_errors_option;
        if let Some(sample) = args.value_of("sample_errors") {
            sample_errors_option = Some(sample.parse::<usize>()?)
        } else {
            sample_errors_option = None
        }

        let constant_errors_option;
        if let Some(constant) = args.value_of("constant_errors") {
            constant_errors_option = Some(constant.parse::<usize>()?)
        } else {
            constant_errors_option = None
        }

        let merge_output;
        if args.is_present("merge_output") {
            merge_output = true
        } else {
            merge_output = false
        }
        let fastq = args.value_of("fastq").unwrap().to_string();
        let format = args.value_of("sequence_format").unwrap().to_string();
        let output_dir = args.value_of("output_dir").unwrap().to_string();
        let threads = args.value_of("threads").unwrap().parse::<u8>().unwrap();
        let prefix = args.value_of("prefix").unwrap().to_string();

        Ok(Args {
            fastq,
            format,
            sample_barcodes_option,
            bb_barcodes_option,
            output_dir,
            threads,
            prefix,
            merge_output,
            bb_errors_option,
            sample_errors_option,
            constant_errors_option,
        })
    }
}
