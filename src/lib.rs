pub mod del_info;
pub mod parse_sequences;

// use flate2::read::GzDecoder;
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
    sync::{Arc, Mutex},
};

use itertools::Itertools;

/// Reads in the FASTQ file line by line, then pushes every 2 out of 4 lines, which corresponds to the sequence line, into a Vec that is passed to other threads
///
/// FASTQ format:
/// Line 1: Sequence ID
/// Line 2: DNA sequence
/// Line 3: +
/// Line 4: Quality score
pub fn read_fastq(
    fastq: String,
    seq_clone: Arc<Mutex<Vec<String>>>,
    exit_clone: Arc<Mutex<bool>>,
) -> Result<(), Box<dyn Error>> {
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
                    while seq_clone.lock().unwrap().len() >= 10000 {
                        // if threads have failed exit out of this thread
                        if *exit_clone.lock().unwrap() {
                            break;
                        }
                    }
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

            // if threads have failed exit out of this thread
            if *exit_clone.lock().unwrap() {
                break;
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
    results: Arc<Mutex<HashMap<String, HashMap<String, u32>>>>,
    bb_num: usize,
    bb_hashmap_option: Option<HashMap<usize, HashMap<String, String>>>,
    prefix: String,
) -> Result<(), Box<dyn Error>> {
    let results_hasmap = results.lock().unwrap(); // get the results

    // Create a comma separated header.  First columns are the building block bumbers, 'BB_#'.  The last header is 'Count'
    let mut header = "BB_1".to_string();
    for num in 1..bb_num {
        header.push_str(&format!(",BB_{}", num + 1))
    }
    header.push_str(",Count\n");

    // create the directory variable to join the file to
    let directory = Path::new(&output_dir);
    for sample_id in results_hasmap.keys() {
        // get the sample results
        let sample_counts_hash = results_hasmap.get(sample_id).unwrap();

        let file_name;
        // If no sample names are supplied, save as all counts, otherwise as sample name counts
        if results_hasmap.keys().count() == 1 && sample_id == "Unknown_sample_name" {
            file_name = format!("{}{}", prefix, "_all_counts.csv");
        } else {
            // create the filename as the sample_id_counts.csv
            file_name = format!("{}_{}{}", prefix, sample_id, "_counts.csv");
        }
        // join the filename with the directory to create the full path
        let output_path = directory.join(file_name);
        let mut output = File::create(output_path)?; // Create the output file

        output.write_all(header.as_bytes())?; // Write the header to the file

        // Iterate through all results and write as comma separated.  The keys within the hashmap are already comma separated
        // If there is an included building block barcode file, it is converted here
        if let Some(ref bb_hashmap) = bb_hashmap_option {
            for (code, count) in sample_counts_hash.iter() {
                // Convert the building block DNA barcodes and join them back to comma separated
                let converted = code
                    .split(",")
                    .enumerate()
                    .map(|(bb_index, bb_barcode)| {
                        let actual_bb_num = bb_index + 1;
                        let barcode_hash = bb_hashmap.get(&actual_bb_num).unwrap();
                        return barcode_hash.get(bb_barcode).unwrap().to_string();
                    })
                    .join(",");
                let row = format!("{},{}\n", converted, count);
                output.write_all(row.as_bytes())?;
            }
        } else {
            for (code, count) in sample_counts_hash.iter() {
                let row = format!("{},{}\n", code, count);
                output.write_all(row.as_bytes())?;
            }
        }
    }
    Ok(())
}
