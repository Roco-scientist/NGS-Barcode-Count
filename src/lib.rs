pub mod del_info;
pub mod parse_sequences;

// use flate2::read::GzDecoder;
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    sync::{Arc, Mutex},
};

/// Reads in the FASTQ file line by line, then pushes every 2 out of 4 lines, which corresponds to the sequence line, into a Vec that is passed to other threads
///
/// FASTQ format:
/// Line 1: Sequence ID
/// Line 2: DNA sequence
/// Line 3: +
/// Line 4: Quality score
pub fn read_fastq(fastq: String, seq_clone: Arc<Mutex<Vec<String>>>) -> Result<(), Box<dyn Error>> {
    let fastq_file = File::open(fastq.clone())?; // open file

    let mut line_num = 1; // start line to know that each 2nd of 4 lines is pulled
                          // If the file is not zipped, proceed.  Still need to work on opening a zipped file
    if !fastq.ends_with("gz") {
        let mut total_reads = 0;

        // go line by line
        for line in BufReader::new(fastq_file).lines() {
            // if it is the sequence line which is line 2
            if line_num == 2 {
                if let Ok(sequence_data) = line {
                    // Pause if there are already 10000 sequences in the vec so memory is not overloaded
                    while seq_clone.lock().unwrap().len() >= 10000 {}
                    // Insert the sequence into the vec.  This will be popped out by other threads
                    seq_clone.lock().unwrap().insert(0, sequence_data);
                }
                // Add to read count to print numnber of sequences read by this thread
                total_reads += 1;
                if total_reads % 1000 == 0 {
                    print!("Total sequences: {}\r", total_reads);
                }
            }

            // increase line number and if it has passed line 4, reset to 1
            line_num += 1;
            if line_num == 5 {
                line_num = 1
            }
        }
    } else {
        // TODO make this work with gzip files
        // let mut buffer = [0; 1000000000];
        // fastq_file.read(&buffer)?;
        // let flate_line = GzDecoder::new(buffer);
        // seq_clone.lock().unwrap().insert(0, flate_line);
    }
    println!();
    Ok(())
}

/// Writes counts result to CSV file
pub fn output_counts(
    output_dir: String,
    results: Arc<Mutex<HashMap<String, u32>>>,
    bb_num: usize,
) -> Result<(), Box<dyn Error>> {
    let results_hasmap = results.lock().unwrap(); // get the results
    let mut output = File::create(format!("{}{}", output_dir, "Counts.csv"))?; // Create the output file

    // Create a comma separated header.  First column 'Sample_ID', the following columns 'BB_#'.  The last header is 'Count'
    let mut header = "Sample_ID".to_string();
    for num in 0..bb_num {
        header.push_str(&format!(",BB_{}", num + 1))
    }
    header.push_str(",Count\n");
    output.write_all(header.as_bytes())?; // Write the header to the file

    // Iterate through all results and write as comma separated.  The keys within the hashmap are already comma separated
    for (code, count) in results_hasmap.iter() {
        let row = format!("{},{}\n", code, count);
        output.write_all(row.as_bytes())?;
    }
    Ok(())
}
