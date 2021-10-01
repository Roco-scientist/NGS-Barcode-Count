use chrono::Local;
use std::sync::{
    atomic::{AtomicBool, AtomicU32, Ordering},
    Arc, Mutex,
};

fn main() {
    // Start a clock to measure how long the algorithm takes
    let start_time = Local::now();

    // get the argument inputs
    let args = barcode::Args::new().unwrap_or_else(|err| panic!("Argument error: {}", err));

    let sequence_format = barcode::barcode_info::SequenceFormat::new(args.format.clone())
        .unwrap_or_else(|err| panic!("sequence format error: {}", err));
    sequence_format.display_format();

    // Start getting the barcode conversion with the BarcodeConversions struct
    let mut barcode_conversions = barcode::barcode_info::BarcodeConversions::new();
    // Create a hashmap of the sample barcodes in order to convert sequence to sample ID
    if let Some(ref samples) = args.sample_barcodes_option {
        barcode_conversions
            .sample_barcode_file_conversion(samples)
            .unwrap();
        barcode_conversions.get_sample_seqs();
    }

    // Create a results struct that will contain the counts.  This is passed between threads
    let results = Arc::new(Mutex::new(barcode::barcode_info::Results::new(
        &barcode_conversions.samples_barcode_hash,
        sequence_format.random_barcode,
    )));

    // Create a hashmap of the building block barcodes in order to convert sequence to building block
    if let Some(ref barcodes) = args.barcodes_option {
        barcode_conversions
            .barcode_file_conversion(barcodes, sequence_format.barcode_num)
            .unwrap();
        barcode_conversions.get_barcode_seqs();
    }

    // Create a sequencing errors Struct to track errors.  This is passed between threads
    let mut sequence_errors = barcode::barcode_info::SequenceErrors::new();

    // Create a passed exit passed variable to stop reading when a thread has panicked
    let exit = Arc::new(AtomicBool::new(false));

    // Create a MaxSeqErrors struct which holds how many sequencing errors are allowed for each sequencing region
    let mut max_errors = barcode::barcode_info::MaxSeqErrors::new(
        args.sample_errors_option,
        sequence_format.sample_length_option().unwrap(),
        args.barcodes_errors_option,
        sequence_format.barcode_lengths().unwrap(),
        args.constant_errors_option,
        sequence_format.constant_region_length(),
        args.min_average_quality_score,
    )
    .unwrap_or_else(|err| panic!("Max Sequencing Errors error: {}", err));
    // Display region sizes and errors allowed
    max_errors.display();

    let total_reads_arc = Arc::new(AtomicU32::new(0));
    // Start the multithreading scope
    rayon::scope(|s| {
        // Create a sequence vec which will have sequences entered by the reading thread, and sequences removed by the processing threads
        let seq = Arc::new(Mutex::new(Vec::new()));
        // Create a passed variable to let the processing threads know the reading thread is done
        let finished = Arc::new(AtomicBool::new(false));

        // Clone variables that are needed to be passed into the reading thread and create the reading thread
        let seq_clone = Arc::clone(&seq);
        let finished_clone = Arc::clone(&finished);
        let exit_clone = Arc::clone(&exit);
        let fastq = args.fastq.clone();
        let total_reads_arc_clone = Arc::clone(&total_reads_arc);
        s.spawn(move |_| {
            barcode::io::read_fastq(fastq, seq_clone, exit_clone, total_reads_arc_clone)
                .unwrap_or_else(|err| {
                    finished_clone.store(true, Ordering::Relaxed);
                    panic!("Error: {}", err)
                });
            finished_clone.store(true, Ordering::Relaxed);
        });

        let shared_mut =
            barcode::parse_sequences::SharedMutData::new(seq, finished, Arc::clone(&results));
        // Create processing threads.  One less than the total threads because of the single reading thread
        for _ in 1..args.threads {
            // Clone all variables needed to pass into each thread
            let shared_mut_clone = shared_mut.arc_clone();
            let sequence_errors_clone = sequence_errors.arc_clone();
            let sequence_format_clone = sequence_format.clone();
            let exit_clone = &exit;
            let max_errors_clone = max_errors.clone();
            let sample_seqs_clone = barcode_conversions.sample_seqs.clone();
            let counted_barcode_seqs_clone = barcode_conversions.counted_barcode_seqs.clone();
            let min_quality_score = args.min_average_quality_score;

            // Create a processing thread
            s.spawn(move |_| {
                let mut parser = barcode::parse_sequences::SequenceParser::new(
                    shared_mut_clone,
                    sequence_errors_clone,
                    sequence_format_clone,
                    max_errors_clone,
                    sample_seqs_clone,
                    counted_barcode_seqs_clone,
                    min_quality_score,
                );
                parser.parse().unwrap_or_else(|err| {
                    exit_clone.store(true, Ordering::Relaxed);
                    panic!("Compute thread panic error: {}", err)
                });
            })
        }
    });

    // Print sequencing error counts to stdout
    sequence_errors.display();

    // Get the end time and print compute time for the algorithm
    let elapsed_time = Local::now() - start_time;
    println!(
        "Compute time: {} hours, {} minutes, {}.{} seconds",
        elapsed_time.num_hours(),
        elapsed_time.num_minutes() % 60,
        elapsed_time.num_seconds() % 60,
        elapsed_time.num_milliseconds() - (elapsed_time.num_seconds() * 1000)
    );
    println!();

    println!("Writing counts");
    println!();
    let mut output = barcode::io::Output::new(
        results,
        sequence_format,
        barcode_conversions.counted_barcodes_hash,
        barcode_conversions.samples_barcode_hash,
        args,
    )
    .unwrap_or_else(|err| panic!("Output error: {}", err));
    output
        .write_files()
        .unwrap_or_else(|err| panic!("Writing error: {}", err));
    // Get the end time and print total time for the algorithm
    output
        .write_stats(start_time, max_errors, sequence_errors, total_reads_arc)
        .unwrap_or_else(|err| panic!("Writing stats error: {}", err));

    // Get the end time and print total time for the algorithm
    let elapsed_time = Local::now() - start_time;
    println!(
        "Total time: {} hours, {} minutes, {}.{} seconds",
        elapsed_time.num_hours(),
        elapsed_time.num_minutes() % 60,
        elapsed_time.num_seconds() % 60,
        elapsed_time.num_milliseconds() - (elapsed_time.num_seconds() * 1000)
    );
}
