use chrono::{DateTime, Local};
use custom_error::custom_error;
use flate2::read::GzDecoder;
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    fs::{File, OpenOptions},
    io::{BufRead, BufReader, Write},
    path::Path,
    sync::{
        atomic::{AtomicBool, AtomicU32, Ordering},
        Arc, Mutex,
    },
};

use itertools::Itertools;

// Errors associated with checking the fastq format to make sure it is correct
custom_error! {FastqError
    NotFastq = "This program only works with *.fastq files and *.fastq.gz files.  The latter is still experimental",
    LeftOver = "Sequence lines left over.  Check code",
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
    seq_clone: Arc<Mutex<Vec<crate::parse_sequences::RawSequenceRead>>>,
    exit_clone: Arc<AtomicBool>,
    total_reads_arc: Arc<AtomicU32>,
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
    fastq_line_reader.display_total_reads()?;
    total_reads_arc.store(fastq_line_reader.total_reads, Ordering::Relaxed);
    println!();
    Ok(())
}

/// A struct with functions for keeping track of read information and to post sequence lines to the shared vector
struct FastqLineReader {
    test: bool,   // whether or not to test the fastq format. Only does this for the first read
    line_num: u8, // the current line number 1-4.  Resets back to 1
    total_reads: u32, // total sequences read within the fastq file
    lines: Vec<String>,
    seq_clone: Arc<Mutex<Vec<crate::parse_sequences::RawSequenceRead>>>, // the vector that is passed between threads which containst the sequences
    exit_clone: Arc<AtomicBool>, // a bool which is set to true when one of the other threads panic.  This is the prevent hanging and is used to exit this thread
}

impl FastqLineReader {
    /// Creates a new FastqLineReader struct
    pub fn new(
        seq_clone: Arc<Mutex<Vec<crate::parse_sequences::RawSequenceRead>>>,
        exit_clone: Arc<AtomicBool>,
    ) -> FastqLineReader {
        FastqLineReader {
            test: true,
            line_num: 1,
            total_reads: 0,
            lines: Vec::new(),
            seq_clone,
            exit_clone,
        }
    }

    /// Reads in the line and either passes to the vec or discards it, depending if it is a sequence line.  Also increments on line count, sequence count etc.
    pub fn read_and_post(&mut self, line: String) -> Result<(), Box<dyn Error>> {
        // Pause if there are already 10000 sequences in the vec so memory is not overloaded
        while self.seq_clone.lock().unwrap().len() >= 10000 {
            // if threads have failed exit out of this thread
            if self.exit_clone.load(Ordering::Relaxed) {
                break;
            }
        }
        if self.line_num == 1 {
            if !self.lines.is_empty() {
                println!("Lines: {:?}", self.lines);
                self.lines.clear();
                return Err(Box::new(FastqError::LeftOver));
            }
            self.lines.push(line);
        } else if self.line_num == 2 {
            self.lines.push(line);
            // Add to read count to print numnber of sequences read by this thread
            self.total_reads += 1;
            if self.total_reads % 1000 == 0 {
                self.display_total_reads()?;
            }
        } else if self.line_num == 3 {
            self.lines.push(line);
        } else if self.line_num == 4 {
            let line_3 = self.lines.pop().unwrap_or(String::new());
            let line_2 = self.lines.pop().unwrap_or(String::new());
            let line_1 = self.lines.pop().unwrap_or(String::new());
            let raw_read =
                crate::parse_sequences::RawSequenceRead::new(line_1, line_2, line_3, line);
            if self.test {
                raw_read.check_fastq_format()?;
                self.test = false;
            }
            // Insert the sequence into the vec.  This will be popped out by other threads
            self.seq_clone.lock().unwrap().insert(0, raw_read);
        }
        // increase line number and if it has passed line 4, reset to 1
        self.line_num += 1;
        if self.line_num == 5 {
            self.line_num = 1
        }
        Ok(())
    }

    /// Displays the total reads so far.  Used while reading to incrementally display, then used after finished reading the file to display total sequences that were read
    pub fn display_total_reads(&self) -> Result<(), Box<dyn Error>> {
        print!("Total sequences:             {}\r", self.total_reads);
        std::io::stdout().flush()?;
        Ok(())
    }
}

/// A struct setup to output results and stat information into files
pub struct Output {
    results: crate::barcode_info::Results,
    sequence_format: crate::barcode_info::SequenceFormat,
    counted_barcodes_hash: Vec<HashMap<String, String>>,
    samples_barcode_hash: HashMap<String, String>,
    merged_output_file_option: Option<File>,
    compounds_written: HashSet<String>,
    args: crate::Args,
    output_files: Vec<String>,
}

impl Output {
    pub fn new(
        results_arc: Arc<Mutex<crate::barcode_info::Results>>,
        sequence_format: crate::barcode_info::SequenceFormat,
        counted_barcodes_hash: Vec<HashMap<String, String>>,
        samples_barcode_hash: HashMap<String, String>,
        args: crate::Args,
    ) -> Result<Output, Box<dyn Error>> {
        let results = Arc::try_unwrap(results_arc).unwrap().into_inner().unwrap();
        Ok(Output {
            results,
            sequence_format,
            counted_barcodes_hash,
            samples_barcode_hash,
            merged_output_file_option: None,
            compounds_written: HashSet::new(),
            args,
            output_files: Vec::new(),
        })
    }

    /// Sets up and writes the results file.  Works for either with or without a random barcode
    pub fn write_files(&mut self) -> Result<(), Box<dyn Error>> {
        let unknown_sample = "Unknown_sample_name".to_string();
        // Pull all sample IDs from either random hashmap or counts hashmap
        let mut sample_barcodes = match self.results.format_type {
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

        // If there was a sample conversion file, sort the barcodes by the sample IDs so that the columns for the merged file are in order
        if !self.samples_barcode_hash.is_empty() {
            sample_barcodes.sort_by_key(|barcode| {
                self.samples_barcode_hash
                    .get(barcode)
                    .unwrap_or(&unknown_sample)
            })
        }

        // create the directory variable to join the file to
        let output_dir = self.args.output_dir.clone();
        let directory = Path::new(&output_dir);

        let mut header = self.create_header();
        // If merged called, create the header with the sample names as columns and write
        if self.args.merge_output {
            if !self.samples_barcode_hash.is_empty() {
                // Create the merge file and push the header, if merged called within arguments
                let merged_file_name = format!("{}{}", self.args.prefix, "_counts.all.csv");
                self.output_files.push(merged_file_name.clone());
                let merged_output_path = directory.join(merged_file_name);

                self.merged_output_file_option = Some(File::create(merged_output_path)?);
                let mut merged_header = header.clone();
                for sample_barcode in &sample_barcodes {
                    // Get the sample name from the sample barcode
                    let sample_name = self
                        .samples_barcode_hash
                        .get(sample_barcode)
                        .unwrap_or(&unknown_sample);
                    merged_header.push(',');
                    merged_header.push_str(sample_name);
                }
                merged_header.push('\n');
                self.merged_output_file_option
                    .as_ref()
                    .unwrap()
                    .write_all(merged_header.as_bytes())?;
            } else {
                eprintln!("Merged file cannot be created without multiple sample barcodes");
                println!()
            }
        }

        // Crate the header to be used with each sample file.  This is just Barcode_1..Barcode_n and Count
        header.push_str(",Count\n");

        for sample_barcode in &sample_barcodes {
            let file_name;
            if !self.samples_barcode_hash.is_empty() {
                let sample_name = self
                    .samples_barcode_hash
                    .get(sample_barcode)
                    .unwrap_or(&unknown_sample);
                file_name = format!("{}_{}{}", self.args.prefix, sample_name, "_counts.csv");
            } else {
                file_name = format!("{}{}", self.args.prefix, "_all_counts.csv");
            }
            self.output_files.push(file_name.clone());
            // join the filename with the directory to create the full path
            let output_path = directory.join(file_name);
            let mut output = File::create(output_path)?; // Create the output file

            output.write_all(header.as_bytes())?; // Write the header to the file
            match self.results.format_type {
                crate::barcode_info::FormatType::RandomBarcode => {
                    self.write_random(sample_barcode, &sample_barcodes, &mut output)?
                }
                crate::barcode_info::FormatType::NoRandomBarcode => {
                    self.write_counts(sample_barcode, &sample_barcodes, &mut output)?
                }
            }
        }
        Ok(())
    }

    /// Creates the file header string for column headers
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

    /// Writes the files for when a random barcode is included
    fn write_random(
        &mut self,
        sample_barcode: &str,
        sample_barcodes: &[String],
        output: &mut File,
    ) -> Result<(), Box<dyn Error>> {
        let sample_random_hash = self.results.random_hashmap.get(sample_barcode).unwrap();
        // Iterate through all results and write as comma separated.  The keys within the hashmap are already comma separated
        // If there is an included building block barcode file, it is converted here
        for (code, random_barcodes) in sample_random_hash.iter() {
            let written_barcodes;
            if !self.counted_barcodes_hash.is_empty() {
                // Convert the building block DNA barcodes and join them back to comma separated
                written_barcodes = convert_code(code, &self.counted_barcodes_hash);
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
                    for sample_barcode in sample_barcodes {
                        merged_row.push(',');
                        merged_row.push_str(
                            &self
                                .results
                                .random_hashmap
                                .get(sample_barcode)
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

    /// Writes the files for when a random barcode is not included
    fn write_counts(
        &mut self,
        sample_barcode: &str,
        sample_barcodes: &[String],
        output: &mut File,
    ) -> Result<(), Box<dyn Error>> {
        let sample_counts_hash = self.results.count_hashmap.get(sample_barcode).unwrap();
        for (code, count) in sample_counts_hash.iter() {
            let written_barcodes;
            if !self.counted_barcodes_hash.is_empty() {
                // Convert the building block DNA barcodes and join them back to comma separated
                written_barcodes = convert_code(code, &self.counted_barcodes_hash);
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
                    for sample_barcode in sample_barcodes {
                        merged_row.push(',');
                        merged_row.push_str(
                            &self
                                .results
                                .count_hashmap
                                .get(sample_barcode)
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
    pub fn write_stats(
        &self,
        start_time: DateTime<Local>,
        max_sequence_errors: crate::barcode_info::MaxSeqErrors,
        mut seq_errors: crate::barcode_info::SequenceErrors,
        total_reads: Arc<AtomicU32>,
    ) -> Result<(), Box<dyn Error>> {
        // Create the stat file name
        let output_dir = self.args.output_dir.clone();
        let directory = Path::new(&output_dir);
        let today = Local::today().format("%Y-%m-%d").to_string();
        let stat_filename = directory.join(format!("{}_barcode_stats.txt", &today));
        // Make the stat file and make it an appending function
        let mut stat_file = OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open(stat_filename)?;

        // Get the total time the program took to run
        let now = Local::now();
        let elapsed_time = now - start_time;
        // Write the time information to the stat file
        stat_file.write_all(
            format!(
                "Start: {}\nFinish: {}\nTotal time: {} hours, {} minutes, {}.{} seconds\n\n",
                start_time.format("%Y-%m-%d %H:%M:%S").to_string(),
                now.format("%Y-%m-%d %H:%M:%S").to_string(),
                elapsed_time.num_hours(),
                elapsed_time.num_minutes() % 60,
                elapsed_time.num_seconds() % 60,
                elapsed_time.num_milliseconds() - (elapsed_time.num_seconds() * 1000)
            )
            .as_bytes(),
        )?;
        // Write the input file information
        stat_file.write_all(
            format!(
                "-Input files-\nFastq: {}\nFormat: {}\nSamples: {}\nBarcodes: {}\n\n",
                self.args.fastq,
                self.args.format,
                self.args
                    .sample_barcodes_option
                    .as_ref()
                    .unwrap_or(&"None".to_string()),
                self.args
                    .barcodes_option
                    .as_ref()
                    .unwrap_or(&"None".to_string())
            )
            .as_bytes(),
        )?;
        // Record the files that were created
        stat_file
            .write_all(format!("Output files: {}\n\n", self.output_files.join(", ")).as_bytes())?;
        // Record the barcode information
        stat_file.write_all(
            format!(
                "-Barcode info-\n{}\n\n",
                max_sequence_errors.display_string()
            )
            .as_bytes(),
        )?;
        // Record the total reads and errors
        stat_file.write_all(
            format!(
                "-Results-\nTotal sequences:             {}\n{}\n\n",
                total_reads.load(Ordering::Relaxed),
                seq_errors.display_string()
            )
            .as_bytes(),
        )?;
        // Close the writing with dashes so that it is separated from the next analysis if it is done on the same day
        stat_file.write_all("--------------------------------------------------------------------------------------------------\n\n\n".as_bytes())?;
        Ok(())
    }
}

/// Converst the DNA sequence from counted barcodes to the ID
fn convert_code(code: &str, barcodes_hashmap: &[HashMap<String, String>]) -> String {
    code.split(',')
        .enumerate()
        .map(|(barcode_index, barcode)| {
            let barcode_hash = &barcodes_hashmap[barcode_index];
            return barcode_hash.get(barcode).unwrap().to_string();
        })
        .join(",")
}

pub fn convert_sample_barcode(
    sample_barcode: &str,
    sample_barcodes_hash: &HashMap<String, String>,
) -> String {
    if let Some(sample_results) = sample_barcodes_hash.get(sample_barcode) {
        sample_results.to_string()
    } else {
        "Unknown_sample_name".to_string()
    }
}
