use custom_error::custom_error;
use flate2::read::GzDecoder;
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
    path::Path,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
};

use itertools::Itertools;

// Errors associated with checking the fastq format to make sure it is correct
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
    exit_clone: Arc<AtomicBool>,
) -> Result<(), Box<dyn Error>> {
    let fastq_file = File::open(fastq.clone())?; // open file

    // Create a fastq line reader which keeps track of line number, reads, and posts the sequence to the shared vector
    let mut fastq_line_reader = FastqLineReader::new(seq_clone, exit_clone);

    // If the file is not gzipped use BufReader to read in lines
    if !fastq.ends_with("fastq.gz") {
        // If the file does not end with fastq, return with an error
        if !fastq.ends_with("fastq") {
            return Err(Box::new(FastqError::NotFastq));
        }

        // go line by line
        for line_result in BufReader::new(fastq_file).lines() {
            let line = line_result?;
            // post the line to the shared vector and keep track of the number of sequences etc
            fastq_line_reader.read_and_post(line)?;
        }
    } else {
        println!("Warning: gzip files is still experimental.  The program may stop reading early. Best results come from using a decompressed fastq file\n");
        println!();
        // stream in first by decoding with GzDecoder, the reading into buffer
        let mut reader = BufReader::new(GzDecoder::new(fastq_file));

        // artificially set the read response to 10.  The first number does not matter
        let mut read_response = 10;
        // continue reading until there is a response of 0, which indicates the end of file.  This may be where some gzipped files abrupty end
        while read_response != 0 {
            let mut line = String::new();
            // move the read line to the line variable and get the response to check if it is 0 and therefore the file is done
            read_response = reader.read_line(&mut line)?;
            // post the line to the shared vector and keep track of the number of sequences etc
            fastq_line_reader.read_and_post(line)?;
        }
    }
    // Display the final total read count
    fastq_line_reader.display_total_reads();
    println!();
    Ok(())
}

/// A struct with functions for keeping track of read information and to post sequence lines to the shared vector
struct FastqLineReader {
    test_first_line: bool, // whether or not to keep testing line 1 as a sequence or metadata
    test_fastq_format: bool, // whether or not to test line 2, which should be a sequence
    line_num: u8,          // the current line number 1-4.  Resets back to 1
    total_reads: u64,      // total sequences read within the fastq file
    seq_clone: Arc<Mutex<Vec<String>>>, // the vector that is passed between threads which containst the sequences
    exit_clone: Arc<AtomicBool>, // a bool which is set to true when one of the other threads panic.  This is the prevent hanging and is used to exit this thread
}

impl FastqLineReader {
    /// Creates a new FastqLineReader struct
    pub fn new(seq_clone: Arc<Mutex<Vec<String>>>, exit_clone: Arc<AtomicBool>) -> FastqLineReader {
        FastqLineReader {
            test_first_line: true,
            test_fastq_format: true,
            line_num: 1,
            total_reads: 0,
            seq_clone,
            exit_clone,
        }
    }

    /// Reads in the line and either passes to the vec or discards it, depending if it is a sequence line.  Also increments on line count, sequence count etc.
    pub fn read_and_post(&mut self, line: String) -> Result<(), Box<dyn Error>> {
        // Test the first line for whether or not it is sequence data.  It should be metadata for FASTQ formats
        if self.test_first_line {
            let linetype = test_sequence(&line);
            match linetype {
                LineType::Sequence => return Err(Box::new(FastqError::Line1Seq)),
                LineType::Metadata => (),
            }
            self.test_first_line = false
        }
        // if it is the sequence line which is line 2
        if self.line_num == 2 {
            // test the first sequence line for whether or not it is a sequence and therefor in the correct format
            if self.test_fastq_format {
                let linetype = test_sequence(&line);
                match linetype {
                    LineType::Sequence => (),
                    LineType::Metadata => return Err(Box::new(FastqError::Line2NotSeq)),
                }
                self.test_fastq_format = false
            }
            // Pause if there are already 10000 sequences in the vec so memory is not overloaded
            while self.seq_clone.lock().unwrap().len() >= 10000 {
                // if threads have failed exit out of this thread
                if self.exit_clone.load(Ordering::Relaxed) {
                    break;
                }
            }
            // Insert the sequence into the vec.  This will be popped out by other threads
            self.seq_clone.lock().unwrap().insert(0, line);
            // Add to read count to print numnber of sequences read by this thread
            self.total_reads += 1;
            if self.total_reads % 1000 == 0 {
                self.display_total_reads();
            }
        }

        // increase line number and if it has passed line 4, reset to 1
        self.line_num += 1;
        if self.line_num == 5 {
            self.line_num = 1
        }
        Ok(())
    }

    /// Displays the total reads so far.  Used while reading to incrementally display, then used after finished reading the file to display total sequences that were read
    pub fn display_total_reads(&self) {
        print!("Total sequences:             {}\r", self.total_reads);
    }
}

/// An enum of linetype to use with test_sequence.  Every line of a FASTQ should either be sequence or metadata for the sequence.
enum LineType {
    Sequence,
    Metadata,
}

/// Tests whether a line within the file String is a sequence by checking if over half of the line contains DNA neceotide letters
fn test_sequence(sequence: &str) -> LineType {
    let sequence_length = sequence.len(); // the the length of the line
    let adenines = sequence.matches('A').count(); // And the amount of each DNA nucleotide
    let guanines = sequence.matches('G').count();
    let cytosines = sequence.matches('C').count();
    let thymines = sequence.matches('T').count();
    let any = sequence.matches('N').count();
    let total_dna = adenines + guanines + cytosines + thymines + any;
    // Check if less than half of the line contains DNA nucleotides.  If so, return that the line is metadata.  Otherwise, a sequence
    if total_dna < sequence_length / 2 {
        return LineType::Metadata;
    }
    LineType::Sequence
}

pub struct Output {
    results: crate::barcode_info::Results,
    sequence_format: crate::barcode_info::SequenceFormat,
    barcodes_hashmap_option: Option<HashMap<u8, HashMap<String, String>>>,
    prefix: String,
    merge_output: bool,
    merged_output_file_option: Option<File>,
    compounds_written: HashSet<String>,
}

impl Output {
    pub fn new(
        results_arc: Arc<Mutex<crate::barcode_info::Results>>,
        sequence_format: crate::barcode_info::SequenceFormat,
        barcodes_hashmap_option: Option<HashMap<u8, HashMap<String, String>>>,
        prefix: String,
        merge_output: bool,
    ) -> Result<Output, Box<dyn Error>> {
        let results = Arc::try_unwrap(results_arc).unwrap().into_inner().unwrap();
        Ok(Output {
            results,
            sequence_format,
            barcodes_hashmap_option,
            prefix,
            merge_output,
            merged_output_file_option: None,
            compounds_written: HashSet::new(),
        })
    }

    pub fn write_files(&mut self, output_dir: String) -> Result<(), Box<dyn Error>> {
        // Pull all sample IDs from either random hashmap or counts hashmap
        let mut sample_ids = match self.results.format_type {
            crate::barcode_info::FormatType::RandomBarcode => self
                .results
                .random_hashmap
                .keys()
                .cloned()
                .collect::<Vec<String>>(),
            crate::barcode_info::FormatType::NoRandomBarcode => self
                .results
                .count_hashmap
                .keys()
                .cloned()
                .collect::<Vec<String>>(),
        };
        sample_ids.sort();

        // create the directory variable to join the file to
        let directory = Path::new(&output_dir);

        let mut header = self.create_header();
        // If merged called, create the header with the sample names as columns and write
        if self.merge_output {
            // Create the merge file and push the header, if merged called within arguments
            let merged_file_name = format!("{}{}", self.prefix, "_counts.all.csv");
            let merged_output_path = directory.join(merged_file_name);

            self.merged_output_file_option = Some(File::create(merged_output_path)?);
            let mut merged_header = header.clone();
            for sample_id in &sample_ids {
                merged_header.push(',');
                merged_header.push_str(sample_id);
            }
            merged_header.push('\n');
            self.merged_output_file_option
                .as_ref()
                .unwrap()
                .write_all(merged_header.as_bytes())?;
        }

        // Crate the header to be used with each sample file.  This is just Barcode_1..Barcode_n and Count
        header.push_str(",Count\n");

        for sample_id in &sample_ids {
            let file_name;
            // If no sample names are supplied, save as all counts, otherwise as sample name counts
            if sample_ids.len() == 1 && sample_id == "Unknown_sample_name" {
                file_name = format!("{}{}", self.prefix, "_all_counts.csv");
            } else {
                // create the filename as the sample_id_counts.csv
                file_name = format!("{}_{}{}", self.prefix, sample_id, "_counts.csv");
            }
            // join the filename with the directory to create the full path
            let output_path = directory.join(file_name);
            let mut output = File::create(output_path)?; // Create the output file

            output.write_all(header.as_bytes())?; // Write the header to the file
            match self.results.format_type {
                crate::barcode_info::FormatType::RandomBarcode => {
                    self.write_random(sample_id, &sample_ids, &mut output)?
                }
                crate::barcode_info::FormatType::NoRandomBarcode => {
                    self.write_counts(sample_id, &sample_ids, &mut output)?
                }
            }
        }
        Ok(())
    }

    fn create_header(&self) -> String {
        // Create a comma separated header.  First columns are the barcodes, 'Barcode_#'.  The last header is 'Count'
        let mut header = String::new();
        if self.sequence_format.barcode_num > 1 {
            header = "Barcode_1".to_string();
            for num in 1..self.sequence_format.barcode_num {
                header.push_str(&format!(",Barcode_{}", num + 1))
            }
        } else {
            header.push_str("Barcode")
        }
        header
    }

    fn write_random(
        &mut self,
        sample_id: &str,
        sample_ids: &[String],
        output: &mut File,
    ) -> Result<(), Box<dyn Error>> {
        let sample_random_hash = self.results.random_hashmap.get(sample_id).unwrap();
        // Iterate through all results and write as comma separated.  The keys within the hashmap are already comma separated
        // If there is an included building block barcode file, it is converted here
        for (code, random_barcodes) in sample_random_hash.iter() {
            let written_barcodes;
            if let Some(ref barcodes_hashmap) = self.barcodes_hashmap_option {
                // Convert the building block DNA barcodes and join them back to comma separated
                written_barcodes = convert_code(code, barcodes_hashmap);
            } else {
                written_barcodes = code.clone();
            }

            // If merge output argument is called, pull data for the compound and write to merged file
            if let Some(ref mut merged_output_file) = self.merged_output_file_option {
                // If the compound has not already been written to the file proceed.  This will happen after the first sample is completed
                let new = self.compounds_written.insert(code.clone());
                if new {
                    // Start a new row with the converted building block barcodes
                    let mut merged_row = written_barcodes.clone();
                    // For every sample, retrieve the count and add to the row with a comma
                    for sample_id in sample_ids {
                        merged_row.push(',');
                        merged_row.push_str(
                            &self
                                .results
                                .random_hashmap
                                .get(sample_id)
                                .unwrap()
                                .get(code)
                                .unwrap_or(&HashSet::new())
                                .len()
                                .to_string(),
                        )
                    }
                    merged_row.push('\n');
                    // write to the merged file
                    merged_output_file.write_all(merged_row.as_bytes())?;
                }
            }
            // Create the row for the sample file and write
            let row = format!("{},{}\n", written_barcodes, random_barcodes.len());
            output.write_all(row.as_bytes())?;
        }
        Ok(())
    }

    fn write_counts(
        &mut self,
        sample_id: &str,
        sample_ids: &[String],
        output: &mut File,
    ) -> Result<(), Box<dyn Error>> {
        let sample_counts_hash = self.results.count_hashmap.get(sample_id).unwrap();
        for (code, count) in sample_counts_hash.iter() {
            let written_barcodes;
            if let Some(ref barcodes_hashmap) = self.barcodes_hashmap_option {
                // Convert the building block DNA barcodes and join them back to comma separated
                written_barcodes = convert_code(code, barcodes_hashmap);
            } else {
                written_barcodes = code.clone();
            }

            // If merge output argument is called, pull data for the compound and write to merged file
            if let Some(ref mut merged_output_file) = self.merged_output_file_option {
                // If the compound has not already been written to the file proceed.  This will happen after the first sample is completed
                let new = self.compounds_written.insert(code.clone());
                if new {
                    // Start a new row with the converted building block barcodes
                    let mut merged_row = written_barcodes.clone();
                    // For every sample, retrieve the count and add to the row with a comma
                    for sample_id in sample_ids {
                        merged_row.push(',');
                        merged_row.push_str(
                            &self
                                .results
                                .count_hashmap
                                .get(sample_id)
                                .unwrap()
                                .get(code)
                                .unwrap_or(&0)
                                .to_string(),
                        )
                    }
                    merged_row.push('\n');
                    // write to the merged file
                    merged_output_file.write_all(merged_row.as_bytes())?;
                }
            }
            // Create the row for the sample file and write
            let row = format!("{},{}\n", written_barcodes, count);
            output.write_all(row.as_bytes())?;
        }
        Ok(())
    }
}

fn convert_code(code: &str, barcodes_hashmap: &HashMap<u8, HashMap<String, String>>) -> String {
    code.split(',')
        .enumerate()
        .map(|(barcode_index, barcode)| {
            let actual_barcode_num = barcode_index as u8 + 1;
            let barcode_hash = barcodes_hashmap.get(&actual_barcode_num).unwrap();
            return barcode_hash.get(barcode).unwrap().to_string();
        })
        .join(",")
}
