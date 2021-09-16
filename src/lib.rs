pub mod del_info;
pub mod parse_sequences;

use flate2::read::GzDecoder;
use std::{
    collections::{HashMap, HashSet},
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
    let mut total_reads = 0u64;
    if !fastq.ends_with("fastq.gz") {
        if !fastq.ends_with("fastq") {
            panic!("This program only works with *.fastq files and *.fastq.gz files.  The latter is still experimental");
        }

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
                    print!("Total sequences:             {}\r", total_reads);
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
        println!("Warning: gzip files is still experimental.  The program may stop reading early. Best results come from using a decompressed fastq file");
        let mut reader = BufReader::new(GzDecoder::new(fastq_file));
        let mut total_reads = 0;

        // go line by line
        // reader.read_line(&mut line);
        let mut read_response = 10;
        while read_response != 0 {
            let mut line = String::new();
            read_response = reader.read_line(&mut line)?;
            //            if read_response == 0 {
            //               line = String::new();
            //              read_response = reader.read_line(&mut line)?;
            //         }
            // println!("Read response: {}\t{}", &read_response, &line);
            //print!("{}", &line);
            // if it is the sequence line which is line 2
            if line_num == 2 {
                // Pause if there are already 10000 sequences in the vec so memory is not overloaded
                while seq_clone.lock().unwrap().len() >= 10000 {}
                // Insert the sequence into the vec.  This will be popped out by other threads
                seq_clone.lock().unwrap().insert(0, line.clone());
                // Add to read count to print numnber of sequences read by this thread
                total_reads += 1;
                if total_reads % 1000 == 0 {
                    print!("Total sequences:             {}\r", total_reads);
                }
            }

            // increase line number and if it has passed line 4, reset to 1
            line_num += 1;
            if line_num == 5 {
                line_num = 1
            }
        }
    }
    print!("Total sequences:             {}\r", total_reads);
    println!();
    Ok(())
}

/// Writes counts result to CSV file
pub fn output_counts(
    output_dir: String,
    results: Arc<Mutex<HashMap<String, HashMap<String, u32>>>>,
    sequence_format: crate::del_info::SequenceFormat,
    bb_hashmap_option: Option<HashMap<usize, HashMap<String, String>>>,
    prefix: String,
    merge_output: bool,
) -> Result<(), Box<dyn Error>> {
    let results_hashmap = results.lock().unwrap(); // get the results

    let mut sample_ids = results_hashmap.keys().cloned().collect::<Vec<String>>();
    sample_ids.sort();

    // Create a comma separated header.  First columns are the building block bumbers, 'BB_#'.  The last header is 'Count'
    let mut header = "BB_1".to_string();
    for num in 1..sequence_format.bb_num {
        header.push_str(&format!(",BB_{}", num + 1))
    }
    // create the directory variable to join the file to
    let directory = Path::new(&output_dir);

    // Create the merge file and push the header, if merged called within arguments
    let merged_file_name = format!("{}{}", prefix, "_counts.all.csv");
    let merged_output_path = directory.join(merged_file_name);
    let mut merged_output_file = File::create(merged_output_path)?;
    // If merged called, create the header with the sample names as columns and write
    if merge_output {
        let mut merged_header = header.clone();
        for sample_id in &sample_ids {
            merged_header.push_str(",");
            merged_header.push_str(&sample_id);
        }
        merged_header.push_str("\n");
        merged_output_file.write_all(merged_header.as_bytes())?;
    }

    // Crate the header to be used with each sample file.  This is just BB_1..BB_n and Count
    header.push_str(",Count\n");

    // Create a HashSet for if there is merging to check what compounds have been written to the merged file
    let mut compounds_written = HashSet::new();
    for sample_id in &sample_ids {
        // get the sample results
        let sample_counts_hash = results_hashmap.get(sample_id).unwrap();

        let file_name;
        // If no sample names are supplied, save as all counts, otherwise as sample name counts
        if results_hashmap.keys().count() == 1 && sample_id == "Unknown_sample_name" {
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

                // If merge output argument is called, pull data for the compound and write to merged file
                if merge_output {
                    // If the compound has not already been written to the file proceed.  This will happen after the first sample is completed
                    if !compounds_written.contains(code) {
                        // Add the compound to the hashset so that it is not repeated later
                        compounds_written.insert(code);
                        // Start a new row with the converted building block barcodes
                        let mut merged_row = converted.clone();
                        // For every sample, retrieve the count and add to the row with a comma
                        for sample_id in &sample_ids {
                            merged_row.push_str(",");
                            merged_row.push_str(
                                &results_hashmap
                                    .get(sample_id)
                                    .unwrap()
                                    .get(code)
                                    .unwrap_or(&0)
                                    .to_string(),
                            )
                        }
                        merged_row.push_str("\n");
                        // write to the merged file
                        merged_output_file.write_all(merged_row.as_bytes())?;
                    }
                }
                // Create the row for the sample file and write
                let row = format!("{},{}\n", converted, count);
                output.write_all(row.as_bytes())?;
            }
        } else {
            // If there is no building block barcode conversion, write row with the DNA barcode instead
            for (code, count) in sample_counts_hash.iter() {
                let row = format!("{},{}\n", code, count);
                output.write_all(row.as_bytes())?;

                // If merge output argument is called, pull data for the compound and write to merged file
                if merge_output {
                    // If the compound has not already been written to the file proceed.  This will happen after the first sample is completed
                    if !compounds_written.contains(code) {
                        // Add the compound to the hashset so that it is not repeated later
                        compounds_written.insert(code);
                        // Start a new row
                        let mut merged_row = code.clone();
                        // For every sample, retrieve the count and add to the row with a comma
                        for sample_id in &sample_ids {
                            merged_row.push_str(",");
                            merged_row.push_str(
                                &results_hashmap
                                    .get(sample_id)
                                    .unwrap()
                                    .get(code)
                                    .unwrap_or(&0)
                                    .to_string(),
                            )
                        }
                        merged_row.push_str("\n");
                        // write to the merged file
                        merged_output_file.write_all(merged_row.as_bytes())?;
                    }
                }
            }
        }
    }
    Ok(())
}
