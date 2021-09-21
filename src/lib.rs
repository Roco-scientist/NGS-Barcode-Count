pub mod barcode_info;
pub mod parse_sequences;

use custom_error::custom_error;
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

custom_error! {FastqError
    NotFastq = "This program only works with *.fastq files and *.fastq.gz files.  The latter is still experimental",
    Line2NotSeq = "The second line within the FASTQ file is not a sequence. Check the FASTQ format",
    Line1Seq = "The first line within the FASTQ contains DNA sequences.  Check the FASTQ format",
}

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
    let mut test_fastq_format = true; // setup bool to test the first sequence read to make sure it is a DNA sequence
    let mut test_first_line = true; // setup bool to test the first read to make sure it is not a DNA sequence
    if !fastq.ends_with("fastq.gz") {
        if !fastq.ends_with("fastq") {
            return Err(Box::new(FastqError::NotFastq));
        }

        // go line by line
        for line_result in BufReader::new(fastq_file).lines() {
            let line = line_result?;
            if test_first_line {
                let linetype = test_sequence(&line);
                match linetype {
                    LineType::Sequence => return Err(Box::new(FastqError::Line1Seq)),
                    LineType::Metadata => (),
                }
                test_first_line = false
            }
            // if it is the sequence line which is line 2
            if line_num == 2 {
                // test the first sequence line for whether or not it is a sequence and therefor in the correct format
                if test_fastq_format {
                    let linetype = test_sequence(&line);
                    match linetype {
                        LineType::Sequence => (),
                        LineType::Metadata => return Err(Box::new(FastqError::Line2NotSeq)),
                    }
                    test_fastq_format = false
                }
                // Pause if there are already 10000 sequences in the vec so memory is not overloaded
                while seq_clone.lock().unwrap().len() >= 10000 {
                    // if threads have failed exit out of this thread
                    if *exit_clone.lock().unwrap() {
                        break;
                    }
                }
                // Insert the sequence into the vec.  This will be popped out by other threads
                seq_clone.lock().unwrap().insert(0, line);
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
        println!();
        println!("Warning: gzip files is still experimental.  The program may stop reading early. Best results come from using a decompressed fastq file\n");
        let mut reader = BufReader::new(GzDecoder::new(fastq_file));

        // go line by line
        // reader.read_line(&mut line);
        let mut read_response = 10;
        while read_response != 0 {
            let mut line = String::new();
            read_response = reader.read_line(&mut line)?;
            // Test the first line for whether or not it is sequence data.  It should be metadata for FASTQ formats
            if test_first_line {
                let linetype = test_sequence(&line);
                match linetype {
                    LineType::Sequence => return Err(Box::new(FastqError::Line1Seq)),
                    LineType::Metadata => (),
                }
                test_first_line = false
            }
            // if it is the sequence line which is line 2
            if line_num == 2 {
                // test the first sequence line for whether or not it is a sequence and therefor in the correct format
                if test_fastq_format {
                    let linetype = test_sequence(&line);
                    match linetype {
                        LineType::Sequence => (),
                        LineType::Metadata => return Err(Box::new(FastqError::Line2NotSeq)),
                    }
                    test_fastq_format = false
                }
                // Pause if there are already 10000 sequences in the vec so memory is not overloaded
                while seq_clone.lock().unwrap().len() >= 10000 {
                    // if threads have failed exit out of this thread
                    if *exit_clone.lock().unwrap() {
                        break;
                    }
                }
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
            if *exit_clone.lock().unwrap() {
                break;
            }
        }
    }
    print!("Total sequences:             {}\r", total_reads);
    println!();
    Ok(())
}

/// An enum of linetype to use with test_sequence.  Every line of a FASTQ should either be sequence or metadata for the sequence.
enum LineType {
    Sequence,
    Metadata,
}

/// Tests whether a line within the file String is a sequence by checking if over half of the line contains DNA neceotide letters
fn test_sequence(sequence: &String) -> LineType {
    let sequence_length = sequence.len(); // the the length of the line
    let adenines = sequence.matches("A").count(); // And the amount of each DNA nucleotide
    let guanines = sequence.matches("G").count();
    let cytosines = sequence.matches("C").count();
    let thymines = sequence.matches("T").count();
    let total_dna = adenines + guanines + cytosines + thymines;
    // Check if less than half of the line contains DNA nucleotides.  If so, return that the line is metadata.  Otherwise, a sequence
    if total_dna < sequence_length / 2 {
        return LineType::Metadata;
    }
    return LineType::Sequence;
}

/// Writes counts result to CSV file
pub fn output_counts(
    output_dir: String,
    results: Arc<Mutex<HashMap<String, HashMap<String, u32>>>>,
    sequence_format: crate::barcode_info::SequenceFormat,
    barcodes_hashmap_option: Option<HashMap<usize, HashMap<String, String>>>,
    prefix: String,
    merge_output: bool,
) -> Result<(), Box<dyn Error>> {
    let results_hashmap = results.lock().unwrap(); // get the results

    let mut sample_ids = results_hashmap.keys().cloned().collect::<Vec<String>>();
    sample_ids.sort();

    // Create a comma separated header.  First columns are the barcodes, 'Barcode_#'.  The last header is 'Count'
    let mut header = String::new();
    if sequence_format.barcode_num > 1 {
        let mut header = "Barcode_1".to_string();
        for num in 1..sequence_format.barcode_num {
            header.push_str(&format!(",Barcode_{}", num + 1))
        }
    } else {
        header.push_str("Barcode")
    }
    // create the directory variable to join the file to
    let directory = Path::new(&output_dir);

    // Create the merge file and push the header, if merged called within arguments
    let merged_file_name = format!("{}{}", prefix, "_counts.all.csv");
    let merged_output_path = directory.join(merged_file_name);
    let mut merged_output_file_option: Option<File> = None;
    // If merged called, create the header with the sample names as columns and write
    if merge_output {
        merged_output_file_option = Some(File::create(merged_output_path)?);
        let mut merged_header = header.clone();
        for sample_id in &sample_ids {
            merged_header.push_str(",");
            merged_header.push_str(&sample_id);
        }
        merged_header.push_str("\n");
        merged_output_file_option
            .as_ref()
            .unwrap()
            .write_all(merged_header.as_bytes())?;
    }

    // Crate the header to be used with each sample file.  This is just Barcode_1..Barcode_n and Count
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
        if let Some(ref barcodes_hashmap) = barcodes_hashmap_option {
            for (code, count) in sample_counts_hash.iter() {
                // Convert the building block DNA barcodes and join them back to comma separated
                let converted = code
                    .split(",")
                    .enumerate()
                    .map(|(barcode_index, barcode)| {
                        let actual_barcode_num = barcode_index + 1;
                        let barcode_hash = barcodes_hashmap.get(&actual_barcode_num).unwrap();
                        return barcode_hash.get(barcode).unwrap().to_string();
                    })
                    .join(",");

                // If merge output argument is called, pull data for the compound and write to merged file
                if let Some(ref mut merged_output_file) = merged_output_file_option {
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
                if let Some(ref mut merged_output_file) = merged_output_file_option {
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
