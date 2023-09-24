use anyhow::{anyhow, Result};
use chrono::{DateTime, Local};
use num_format::{Locale, ToFormattedString};
use std::{
    fs::{File, OpenOptions},
    io::{stdout, Write},
    path::Path,
    sync::{
        atomic::{AtomicU32, Ordering},
        Arc, Mutex,
    },
};

use ahash::{AHashSet, HashMap, HashMapExt};

use itertools::Itertools;

use crate::{
    arguments::Args,
    info::{
        MaxSeqErrors, Results, ResultsEnrichment, ResultsHashmap, SequenceErrors, SequenceFormat,
    },
};

#[derive(PartialEq, Clone)]
enum EnrichedType {
    Single,
    Double,
    Full,
}

/// A struct setup to output results and stat information into files
pub struct WriteFiles {
    results: Results,
    results_enriched: ResultsEnrichment,
    sequence_format: SequenceFormat,
    counted_barcodes_hash: Vec<HashMap<String, String>>,
    samples_barcode_hash: HashMap<String, String>,
    compounds_written: AHashSet<String>,
    args: Args,
    output_files: Vec<String>,
    output_counts: Vec<usize>,
    merged_count: usize,
    merge_text: String,
    sample_text: String,
}

impl WriteFiles {
    pub fn new(
        results_arc: Arc<Mutex<Results>>,
        sequence_format: SequenceFormat,
        counted_barcodes_hash: Vec<HashMap<String, String>>,
        samples_barcode_hash: HashMap<String, String>,
        args: Args,
    ) -> Result<Self> {
        let results = Arc::try_unwrap(results_arc).unwrap().into_inner().unwrap();
        Ok(WriteFiles {
            results,
            results_enriched: ResultsEnrichment::new(),
            sequence_format,
            counted_barcodes_hash,
            samples_barcode_hash,
            compounds_written: AHashSet::new(),
            args,
            output_files: Vec::new(),
            output_counts: Vec::new(),
            merged_count: 0,
            merge_text: String::new(),
            sample_text: String::new(),
        })
    }

    /// Sets up and writes the results file.  Works for either with or without a random barcode
    pub fn write_counts_files(&mut self) -> Result<()> {
        let unknown_sample = "barcode".to_string();
        // Pull all sample IDs from either random hashmap or counts hashmap
        let mut sample_barcodes = match &self.results.results_hashmap {
            ResultsHashmap::RandomBarcode(random_hashmap) => {
                random_hashmap.keys().cloned().collect::<Vec<String>>()
            }
            ResultsHashmap::NoRandomBarcode(count_hashmap) => {
                count_hashmap.keys().cloned().collect::<Vec<String>>()
            }
        };

        if self.args.enrich {
            self.results_enriched.add_sample_barcodes(&sample_barcodes);
        }

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
            if sample_barcodes.len() == 1 {
                eprintln!("Merged file cannot be created without multiple sample barcodes");
                println!();
                self.args.merge_output = false;
            } else {
                // Create the merge file and push the header
                let mut merged_header = header.clone();
                for sample_barcode in &sample_barcodes {
                    let sample_name = if self.samples_barcode_hash.is_empty() {
                        sample_barcode
                    } else {
                        // Get the sample name from the sample barcode
                        self
                            .samples_barcode_hash
                            .get(sample_barcode)
                            .unwrap_or(&unknown_sample)
                    };
                    merged_header.push(',');
                    merged_header.push_str(sample_name);
                }
                merged_header.push('\n');
                self.merge_text.push_str(&merged_header);
            }
        }

        // Crate the header to be used with each sample file.  This is just Barcode_1..Barcode_n and Count
        header.push_str(",Count\n");

        // For each sample, write the counts file
        for sample_barcode in &sample_barcodes {
            let sample_name = if !self.samples_barcode_hash.is_empty() {
                self
                    .samples_barcode_hash
                    .get(sample_barcode)
                    .unwrap_or(&unknown_sample)
            } else {
                sample_barcode
            };
            let file_name = format!("{}_{}_counts.csv", self.args.prefix, sample_name);
            println!("{}", file_name);
            self.output_files.push(file_name.clone());
            // join the filename with the directory to create the full path
            let output_path = directory.join(file_name);

            self.sample_text.push_str(&header);
            let count =
                self.add_counts_string(sample_barcode, &sample_barcodes, EnrichedType::Full)?;

            let mut output = File::create(output_path)?; // Create the output file
            output.write_all(self.sample_text.as_bytes())?;
            self.sample_text.clear();
            self.output_counts.push(count);
        }
        if self.args.merge_output {
            let merged_file_name = format!("{}{}", self.args.prefix, "_counts.all.csv");
            println!("{}", merged_file_name);
            println!(
                "Barcodes counted: {}",
                self.merged_count.to_formatted_string(&Locale::en)
            );
            self.output_files.push(merged_file_name.clone());
            let merged_output_path = directory.join(merged_file_name);
            let mut merged_output_file = File::create(merged_output_path)?;
            merged_output_file.write_all(self.merge_text.as_bytes())?;
            self.merge_text.clear();
            self.output_counts.insert(0, self.merged_count);
            self.merged_count = 0;
        }
        if self.args.enrich {
            self.write_enriched_files(EnrichedType::Single)?;
            if self.sequence_format.barcode_num > 2 {
                self.write_enriched_files(EnrichedType::Double)?;
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

    /// Writes the files for when a random barcode is not included
    fn add_counts_string(
        &mut self,
        sample_barcode: &str,
        sample_barcodes: &[String],
        enrichment: EnrichedType, // In order to make this non redundant with writing single and double barcodes, this enum determines some aspects
    ) -> Result<usize> {
        let mut hash_holder: HashMap<String, HashMap<String, usize>> = HashMap::new(); // a hodler hash to hold the hashmap from sample_counts_hash for a longer lifetime.  Also used later
                                                                                       // Select from the hashmap connected the the EnrichedType
        let codes = match enrichment {
            EnrichedType::Single => {
                hash_holder = self.results_enriched.single_hashmap.clone();
                hash_holder
                    .get(sample_barcode)
                    .unwrap()
                    .keys()
                    .cloned()
                    .collect::<Vec<String>>()
            }
            EnrichedType::Double => {
                hash_holder = self.results_enriched.double_hashmap.clone();
                hash_holder
                    .get(sample_barcode)
                    .unwrap()
                    .keys()
                    .cloned()
                    .collect::<Vec<String>>()
            }
            EnrichedType::Full => match &self.results.results_hashmap {
                ResultsHashmap::NoRandomBarcode(count_hashmap) => count_hashmap
                    .get(sample_barcode)
                    .unwrap()
                    .keys()
                    .cloned()
                    .collect::<Vec<String>>(),
                ResultsHashmap::RandomBarcode(random_hashmap) => random_hashmap
                    .get(sample_barcode)
                    .unwrap()
                    .keys()
                    .cloned()
                    .collect::<Vec<String>>(),
            },
        };

        let mut barcode_num = 0;
        for (line_num, code) in codes.iter().enumerate() {
            let count = match enrichment {
                EnrichedType::Single => *self
                    .results_enriched
                    .single_hashmap
                    .get(sample_barcode)
                    .unwrap()
                    .get(code)
                    .unwrap(),
                EnrichedType::Double => *self
                    .results_enriched
                    .double_hashmap
                    .get(sample_barcode)
                    .unwrap()
                    .get(code)
                    .unwrap(),
                EnrichedType::Full => match &self.results.results_hashmap {
                    ResultsHashmap::NoRandomBarcode(count_hashmap) => *count_hashmap
                        .get(sample_barcode)
                        .unwrap()
                        .get(code)
                        .unwrap(),
                    ResultsHashmap::RandomBarcode(random_hashmap) => random_hashmap
                        .get(sample_barcode)
                        .unwrap()
                        .get(code)
                        .unwrap()
                        .len(),
                },
            };
            barcode_num = line_num + 1;
            // Print the number counted so far ever 50,000 writes
            if barcode_num % 50000 == 0 {
                print!(
                    "Barcodes counted: {}\r",
                    barcode_num.to_formatted_string(&Locale::en)
                );
                stdout().flush()?;
            }
            let written_barcodes = if enrichment == EnrichedType::Full && !self.counted_barcodes_hash.is_empty() {
                // Convert the building block DNA barcodes and join them back to comma separated
                convert_code(code, &self.counted_barcodes_hash)
            } else {
                code.to_string()
            };

            // If merge output argument is called, pull data for the compound and write to merged file
            if self.args.merge_output {
                // If the compound has not already been written to the file proceed.  This will happen after the first sample is completed
                let new = self.compounds_written.insert(code.to_string());
                if new {
                    self.merged_count += 1;
                    // Start a new row with the converted building block barcodes
                    let mut merged_row = written_barcodes.clone();
                    // For every sample, retrieve the count and add to the row with a comma
                    for sample_barcode in sample_barcodes {
                        merged_row.push(',');
                        // Get teh sample count from the hashmap that corresponds to the EnrichedType.  For single and double, it is the holding hashmap created earlier
                        let sample_count = match enrichment {
                            EnrichedType::Single => hash_holder
                                .get(sample_barcode)
                                .unwrap()
                                .get(code)
                                .unwrap_or(&0)
                                .to_string(),

                            EnrichedType::Double => hash_holder
                                .get(sample_barcode)
                                .unwrap()
                                .get(code)
                                .unwrap_or(&0)
                                .to_string(),

                            EnrichedType::Full => match &self.results.results_hashmap {
                                ResultsHashmap::RandomBarcode(random_hashmap) => random_hashmap
                                    .get(sample_barcode)
                                    .unwrap()
                                    .get(code)
                                    .unwrap_or(&AHashSet::new())
                                    .len()
                                    .to_string(),
                                ResultsHashmap::NoRandomBarcode(count_hashmap) => count_hashmap
                                    .get(sample_barcode)
                                    .unwrap()
                                    .get(code)
                                    .unwrap_or(&0)
                                    .to_string(),
                            },
                        };
                        merged_row.push_str(&sample_count);
                    }
                    merged_row.push('\n');
                    // write to the merged file
                    self.merge_text.push_str(&merged_row);
                }
            }
            // Create the row for the sample file and write
            let row = format!("{},{}\n", written_barcodes, count);
            self.sample_text.push_str(&row);
            // If enrichment type is Full, which is neither single nor double for adding string,
            // and enrich is called.  Add 1 and 2 synthon enrichment.  This is becuase this smae
            // method is called to create the 1 and 2 synthon strings, and therefore should only
            // run when Full is used
            if enrichment == EnrichedType::Full && self.args.enrich {
                self.results_enriched
                    .add_single(sample_barcode, &written_barcodes, count);
                if self.sequence_format.barcode_num > 2 {
                    self.results_enriched
                        .add_double(sample_barcode, &written_barcodes, count);
                }
            }
        }
        print!(
            "Barcodes counted: {}\r",
            barcode_num.to_formatted_string(&Locale::en)
        );
        println!();
        Ok(barcode_num)
    }

    /// Write enriched files for either single or double barcodes if either flag is called
    fn write_enriched_files(&mut self, enrichment: EnrichedType) -> Result<()> {
        let unknown_sample = "barcode".to_string();
        // Pull all sample IDs from either single or double hashmap, which was added to in either random or counts write
        let mut sample_barcodes = match enrichment {
            EnrichedType::Single => self
                .results_enriched
                .single_hashmap
                .keys()
                .cloned()
                .collect::<Vec<String>>(),
            EnrichedType::Double => self
                .results_enriched
                .double_hashmap
                .keys()
                .cloned()
                .collect::<Vec<String>>(),
            EnrichedType::Full => {
                return Err(anyhow!(
                    "Does not work with Full enrichment type.  Only Single and Double"
                ))
            }
        };

        // If there was a sample conversion file, sort the barcodes by the sample IDs so that the columns for the merged file are in order
        if !self.samples_barcode_hash.is_empty() {
            sample_barcodes.sort_by_key(|barcode| {
                self.samples_barcode_hash
                    .get(barcode)
                    .unwrap_or(&unknown_sample)
            })
        }

        // Create a descriptor for output file names
        let descriptor = match enrichment {
            EnrichedType::Single => "Single",
            EnrichedType::Double => "Double",
            EnrichedType::Full => {
                return Err(anyhow!(
                    "Does not work with Full enrichment type.  Only Single and Double"
                ))
            }
        };

        // create the directory variable to join the file to
        let output_dir = self.args.output_dir.clone();
        let directory = Path::new(&output_dir);

        let mut header = self.create_header();
        // If merged called, create the header with the sample names as columns and write
        if self.args.merge_output {
            let mut merged_header = header.clone();
            for sample_barcode in &sample_barcodes {
                let sample_name = if self.samples_barcode_hash.is_empty() {
                    sample_barcode
                } else {
                    // Get the sample name from the sample barcode
                    self
                        .samples_barcode_hash
                        .get(sample_barcode)
                        .unwrap_or(&unknown_sample)
                };
                merged_header.push(',');
                merged_header.push_str(sample_name);
            }
            merged_header.push('\n');
            self.merge_text.push_str(&merged_header);
        }

        // Crate the header to be used with each sample file.  This is just Barcode_1..Barcode_n and Count
        header.push_str(",Count\n");

        // For each sample, write the enriched file
        for sample_barcode in &sample_barcodes {
            // Create the file_name with the single or double descriptor
            let sample_name = if !self.samples_barcode_hash.is_empty() {
                self
                    .samples_barcode_hash
                    .get(sample_barcode)
                    .unwrap_or(&unknown_sample)
            } else {
                sample_barcode
            };
            let file_name = format!(
                "{}_{}_counts.{}.csv",
                self.args.prefix, sample_name, descriptor
            );
            println!("{}", file_name);
            self.output_files.push(file_name.clone());
            // join the filename with the directory to create the full path
            let output_path = directory.join(file_name);

            self.sample_text.push_str(&header);
            let count =
                self.add_counts_string(sample_barcode, &sample_barcodes, enrichment.clone())?;
            let mut output = File::create(output_path)?; // Create the output file
            output.write_all(self.sample_text.as_bytes())?;
            self.sample_text.clear();
            // add the counts to output to stats later
            self.output_counts.push(count);
        }
        // Add the count of merged barcodes if the flag is called
        if self.args.merge_output {
            // Create the merge file and push the header, if merged called within arguments
            let merged_file_name = format!("{}_counts.all.{}.csv", self.args.prefix, descriptor);
            println!("{}", merged_file_name);
            self.output_files.push(merged_file_name.clone());
            let merged_output_path = directory.join(merged_file_name);
            let mut merged_output_file = File::create(merged_output_path)?;
            merged_output_file.write_all(self.merge_text.as_bytes())?;
            println!(
                "Barcodes counted: {}",
                self.merged_count.to_formatted_string(&Locale::en)
            );
            self.merge_text.clear();
            self.output_counts.insert(
                self.output_counts.len() - sample_barcodes.len(),
                self.merged_count,
            );
            self.merged_count = 0;
        }
        Ok(())
    }

    /// Appends the stats information for record keeping
    pub fn write_stats_file(
        &self,
        start_time: DateTime<Local>,
        max_sequence_errors: MaxSeqErrors,
        seq_errors: SequenceErrors,
        total_reads: Arc<AtomicU32>,
        sequence_format: SequenceFormat,
    ) -> Result<()> {
        // Create the stat file name
        let output_dir = self.args.output_dir.clone();
        let directory = Path::new(&output_dir);
        let stat_filename = directory.join(format!("{}_barcode_stats.txt", self.args.prefix));
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
                "-TIME INFORMATION-\nStart: {}\nFinish: {}\nTotal time: {} hours, {} minutes, {}.{} seconds\n\n",
                start_time.format("%Y-%m-%d %H:%M:%S"),
                now.format("%Y-%m-%d %H:%M:%S"),
                elapsed_time.num_hours(),
                elapsed_time.num_minutes() % 60,
                elapsed_time.num_seconds() % 60,
                millisecond_decimal(elapsed_time)
            )
            .as_bytes(),
        )?;
        // Write the input file information
        stat_file.write_all(
            format!(
                "-INPUT FILES-\nFastq: {}\nFormat: {}\nSamples: {}\nBarcodes: {}\n\n",
                self.args.fastq,
                self.args.format,
                self.args
                    .sample_barcodes_option
                    .as_ref()
                    .unwrap_or(&"None".to_string()),
                self.args
                    .counted_barcodes_option
                    .as_ref()
                    .unwrap_or(&"None".to_string())
            )
            .as_bytes(),
        )?;
        // Record the sequence_format
        stat_file.write_all(format!("{}\n\n", sequence_format).as_bytes())?;
        // Record the barcode information
        stat_file.write_all(format!("{}\n", max_sequence_errors).as_bytes())?;
        // Record the total reads and errors
        stat_file.write_all(
            format!(
                "-RESULTS-\nTotal sequences:             {}\n{}\n\n",
                total_reads
                    .load(Ordering::Relaxed)
                    .to_formatted_string(&Locale::en),
                seq_errors
            )
            .as_bytes(),
        )?;
        // Record the files that were created
        stat_file.write_all("-OUTPUT FILES-\n".as_bytes())?;
        for (file_name, counts) in self.output_files.iter().zip(self.output_counts.iter()) {
            stat_file.write_all(
                format!(
                    "File & barcodes counted: {}\t{}\n",
                    file_name,
                    counts.to_formatted_string(&Locale::en)
                )
                .as_bytes(),
            )?;
        }
        stat_file.write_all("\n".as_bytes())?;
        if self.args.fastq.ends_with("gz") && total_reads.load(Ordering::Relaxed) < 1_000_000 {
            let warning = "WARNING: The program may have stopped early with the gzipped file.  Unzip the fastq.gz and rerun the algorithm on the unzipped fastq file if the number of reads is expected to be above 1,000,000 ";
            println!("\n{}\n", warning);
            stat_file.write_all(format!("\n{}\n", warning).as_bytes())?;
        }
        // Close the writing with dashes so that it is separated from the next analysis if it is done on the same day
        stat_file.write_all("--------------------------------------------------------------------------------------------------\n\n\n".as_bytes())?;
        Ok(())
    }
}

pub fn millisecond_decimal(elapsed_time: chrono::Duration) -> String {
    let milliseconds =
        (elapsed_time.num_milliseconds() - (elapsed_time.num_seconds() * 1000)).to_string();
    let mut final_string = String::new();
    for _ in milliseconds.chars().count()..3 {
        final_string.push('0');
    }
    final_string.push_str(&milliseconds);
    final_string
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
        "barcode".to_string()
    }
}
