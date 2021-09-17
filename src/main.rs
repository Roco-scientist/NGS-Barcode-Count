use chrono::Local;
use clap::{App, Arg};
use num_cpus;
use rayon;
// use rayon::prelude::*;
use std::{
    collections::HashMap,
    sync::{Arc, Mutex},
    time::Instant,
};

fn main() {
    // Start a clock to measure how long the algorithm takes
    let start = Instant::now();

    // get the argument inputs
    let (
        fastq,
        format,
        samples_barcodes,
        bb_barcodes,
        output_dir,
        threads,
        prefix,
        merge_output,
        bb_errors_option,
        sample_errors_option,
        constant_errors_option,
    ) = arguments().unwrap_or_else(|err| panic!("Argument error: {}", err));

    let sequence_format = del::del_info::SequenceFormat::new(format)
        .unwrap_or_else(|err| panic!("sequence format error: {}", err));
    sequence_format.display_format();

    // Create a results hashmap that will contain the counts.  This is passed between threads
    let results = Arc::new(Mutex::new(HashMap::new()));
    // Create a random_barcodes hashmap to keep track of the random barcodes.  This way if it already was found for the same sample and building blocks
    // it will not be counted
    let random_barcodes = Arc::new(Mutex::new(HashMap::new()));

    // Create a hashmap of the sample barcodes in order to convert sequence to sample ID
    let samples_hashmap_option;
    if let Some(samples) = samples_barcodes {
        let samples_hashmap = del::del_info::sample_barcode_file_conversion(samples).unwrap();
        samples_hashmap_option = Some(samples_hashmap);
    } else {
        samples_hashmap_option = None
    }

    // Create a hashmap of the building block barcodes in order to convert sequence to building block
    let bb_hashmap;
    if let Some(bb) = bb_barcodes {
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
        sample_errors_option,
        bb_errors_option,
        constant_errors_option,
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
        s.spawn(move |_| {
            del::read_fastq(fastq, seq_clone, exit_clone).unwrap();
            *finished_clone.lock().unwrap() = true;
        });

        // Create processing threads.  One less than the total threads because of the single reading thread
        for _ in 1..threads {
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
    // Get the end time and print total time for the algorithm
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
        output_dir,
        results,
        sequence_format,
        bb_hashmap,
        prefix,
        merge_output,
    )
    .unwrap();
    // Get the end time and print total time for the algorithm
    let elapsed_time = start.elapsed();
    if elapsed_time.as_secs() < 2 {
        println!("Total time: {} milliseconds", elapsed_time.as_millis());
    } else {
        if elapsed_time.as_secs() > 600 {
            println!("Total time: {} minutes", elapsed_time.as_secs() / 60)
        } else {
            println!("Total time: {} seconds", elapsed_time.as_secs())
        }
    }
}

/// Gets the command line arguments
pub fn arguments() -> Result<
    (
        String,
        String,
        Option<String>,
        Option<String>,
        String,
        u8,
        String,
        bool,
        Option<usize>,
        Option<usize>,
        Option<usize>,
    ),
    Box<dyn std::error::Error>,
> {
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

    let sample_barcodes;
    if let Some(sample) = args.value_of("sample_barcodes") {
        sample_barcodes = Some(sample.to_string())
    } else {
        sample_barcodes = None
    }

    let bb_barcodes;
    if let Some(bb) = args.value_of("bb_barcodes") {
        bb_barcodes = Some(bb.to_string())
    } else {
        bb_barcodes = None
    }

    let bb_errors;
    if let Some(bb) = args.value_of("bb_errors") {
        bb_errors = Some(bb.parse::<usize>()?)
    } else {
        bb_errors = None
    }

    let sample_errors;
    if let Some(sample) = args.value_of("sample_errors") {
        sample_errors = Some(sample.parse::<usize>()?)
    } else {
        sample_errors = None
    }

    let constant_errors;
    if let Some(constant) = args.value_of("constant_errors") {
        constant_errors = Some(constant.parse::<usize>()?)
    } else {
        constant_errors = None
    }

    let merge_output;
    if args.is_present("merge_output") {
        merge_output = true
    } else {
        merge_output = false
    }

    return Ok((
        args.value_of("fastq").unwrap().to_string(),
        args.value_of("sequence_format").unwrap().to_string(),
        sample_barcodes,
        bb_barcodes,
        args.value_of("output_dir").unwrap().to_string(),
        args.value_of("threads").unwrap().parse::<u8>().unwrap(),
        args.value_of("prefix").unwrap().to_string(),
        merge_output,
        bb_errors,
        sample_errors,
        constant_errors,
    ));
}
