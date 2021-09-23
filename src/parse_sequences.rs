use std::{
    collections::HashMap,
    error::Error,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
};

type CountedBarcode = String;
type BarcodeID = String;
type BarcodeBarcodeID = HashMap<CountedBarcode, BarcodeID>;
type BarcodeNum = usize;
type BarcodeNumBarcode = HashMap<BarcodeNum, BarcodeBarcodeID>;
type RandomBarcodes = Vec<String>;
type BarcodeRandomBarcode = HashMap<CountedBarcode, RandomBarcodes>;
type SampleName = String;
type RandomBarcodeHolder = HashMap<SampleName, BarcodeRandomBarcode>;

pub struct SequenceParser {
    seq_clone: Arc<Mutex<Vec<String>>>,
    finished_clone: Arc<AtomicBool>,
    sequence_format_clone: crate::barcode_info::SequenceFormat,
    results_clone: Arc<Mutex<HashMap<String, HashMap<String, u32>>>>,
    random_barcodes_clone: Arc<Mutex<RandomBarcodeHolder>>,
    samples_clone: Option<HashMap<String, String>>,
    barcodes_clone: Option<BarcodeNumBarcode>,
    sequence_errors_clone: Arc<Mutex<crate::barcode_info::SequenceErrors>>,
    max_errors_clone: crate::barcode_info::MaxSeqErrors,
    sample_seqs: Option<Vec<String>>,
    barcodes_seqs_option: Option<Vec<Vec<String>>>,
    raw_sequence: RawSequence,
    seq_match_result_option: Option<SequenceMatchResult>,
}

impl SequenceParser {
    pub fn new(
        seq_clone: Arc<Mutex<Vec<String>>>,
        finished_clone: Arc<AtomicBool>,
        sequence_format_clone: crate::barcode_info::SequenceFormat,
        results_clone: Arc<Mutex<HashMap<String, HashMap<String, u32>>>>,
        random_barcodes_clone: Arc<Mutex<RandomBarcodeHolder>>,
        samples_clone: Option<HashMap<String, String>>,
        barcodes_clone: Option<BarcodeNumBarcode>,
        sequence_errors_clone: Arc<Mutex<crate::barcode_info::SequenceErrors>>,
        max_errors_clone: crate::barcode_info::MaxSeqErrors,
    ) -> SequenceParser {
        SequenceParser {
            seq_clone,
            finished_clone,
            sequence_format_clone,
            results_clone,
            random_barcodes_clone,
            samples_clone,
            barcodes_clone,
            sequence_errors_clone,
            max_errors_clone,
            sample_seqs: None,
            barcodes_seqs_option: None,
            raw_sequence: RawSequence::new("".to_string()),
            seq_match_result_option: None,
        }
    }
    pub fn parse(&mut self) -> Result<(), Box<dyn Error>> {
        self.get_sample_seqs();
        // Get a vec of all possible building block barcodes for error correction
        self.get_barcode_seqs();

        // Loop until there are no sequences left to parse.  These are fed into seq_clone vec by the reader thread
        loop {
            // If there are no sequences in seq_clone, pause, unless the reader thread is finished
            while self.seq_clone.lock().unwrap().is_empty() {
                if self.finished_clone.load(Ordering::Relaxed) {
                    break;
                }
            }
            // If thre are no sequences and the reader thread is finished, break out of the loop
            if self.seq_clone.lock().unwrap().is_empty()
                && self.finished_clone.load(Ordering::Relaxed)
            {
                break;
            }

            self.get_seqeunce();
            // If this returns a Some value proceed
            if !self.raw_sequence.sequence.is_empty() {
                // get the result string match from the sequence.  This is a String with commas between sample_ID and builbding block sequences.
                // This is a convenience for later writing to a file since it will already be comma separated
                self.match_seq()?;
                // If the constant region matched proceed to add to results, otherwise try and fix the constant region
                if let Some(seq_match_result) = &self.seq_match_result_option {
                    // Alwasy add value unless random barcode is included and it has already been found for the sample and building blocks
                    let mut add_value = true;
                    // If there is a random barcode included
                    if let Some(random_barcode) = &seq_match_result.random_barcode_option {
                        // Unlock the random_hashmap
                        let mut random_hashmap = self.random_barcodes_clone.lock().unwrap();

                        // If it does not already have the sample name, insert the sample name -> building_block -> random_barcodes
                        if !random_hashmap.contains_key(&seq_match_result.sample_name) {
                            let mut intermediate_hashmap = HashMap::new();
                            let intermediate_vec = vec![random_barcode.to_string()];
                            intermediate_hashmap
                                .insert(seq_match_result.barcode_string.clone(), intermediate_vec);
                            random_hashmap
                                .insert(seq_match_result.sample_name.clone(), intermediate_hashmap);
                        } else {
                            let barcodes_hashmap = random_hashmap
                                .get_mut(&seq_match_result.sample_name)
                                .unwrap();
                            // If the random hashmap does not have the building blocks yet, insert building_block -> random_barcodes
                            if !barcodes_hashmap.contains_key(&seq_match_result.barcode_string) {
                                let intermediate_vec = vec![random_barcode.to_string()];
                                barcodes_hashmap.insert(
                                    seq_match_result.barcode_string.clone(),
                                    intermediate_vec,
                                );
                            } else {
                                let random_vec = barcodes_hashmap
                                    .get_mut(&seq_match_result.barcode_string)
                                    .unwrap();
                                // else check if the random barcode already used for the sample_name and building_blocks
                                // if the random barcode is already in the vector, change add_value to false
                                // otherqise add the random barcode to the random_barcodes vector
                                if random_vec.contains(&random_barcode) {
                                    add_value = false;
                                    self.sequence_errors_clone.lock().unwrap().duplicated()
                                } else {
                                    random_vec.push(random_barcode.to_string())
                                }
                            }
                        }
                    }
                    // Add 1 count to the results hashmap
                    if add_value {
                        let mut results_hashmap = self.results_clone.lock().unwrap();

                        // If results hashmap does not already contain the sample_name, insert the sanmle_name -> barcodes -> 0
                        if !results_hashmap.contains_key(&seq_match_result.sample_name) {
                            let mut intermediate_hashmap = HashMap::new();
                            intermediate_hashmap.insert(seq_match_result.barcode_string.clone(), 0);
                            results_hashmap
                                .insert(seq_match_result.sample_name.clone(), intermediate_hashmap);
                        }

                        // Insert 0 if the barcodes is not within the sample_name -> barcodes
                        // Then add one regardless
                        *results_hashmap
                            .get_mut(&seq_match_result.sample_name)
                            .unwrap()
                            .entry(seq_match_result.barcode_string.clone())
                            .or_insert(0) += 1;
                        self.sequence_errors_clone.lock().unwrap().correct_match()
                    }
                }
            }
        }
        Ok(())
    }
    fn get_seqeunce(&mut self) {
        // Pop off the last sequence from the seq_clone vec
        if let Some(new_sequence) = self.seq_clone.lock().unwrap().pop() {
            self.raw_sequence = RawSequence::new(new_sequence)
        } else {
            self.raw_sequence = RawSequence::new("".to_string())
        }
    }
    fn get_sample_seqs(&mut self) {
        // Get a vec of all possible sample barcodes for error correction
        if let Some(ref samples) = self.samples_clone {
            self.sample_seqs = Some(samples.keys().map(|key| key.to_string()).collect());
        }
    }
    fn get_barcode_seqs(&mut self) {
        if let Some(ref barcodes) = self.barcodes_clone {
            let mut barcodes_vec = Vec::new();
            let mut barcodes_keys = barcodes.keys().collect::<Vec<&usize>>();
            barcodes_keys.sort();
            for key in barcodes_keys {
                let barcodes_data = barcodes.get(key).unwrap();
                let barcodes = barcodes_data
                    .keys()
                    .map(|key| key.to_string())
                    .collect::<Vec<String>>();
                barcodes_vec.push(barcodes);
            }
            self.barcodes_seqs_option = Some(barcodes_vec);
        }
    }
    /// Does a regex search and captures the barcodes.  Converts the sample barcode to ID.  Returns a String with commas between Sample_ID and
    /// building block barcodes.  This is used as a key within the results vector, where the value can be used as the count
    fn match_seq(&mut self) -> Result<(), Box<dyn Error>> {
        if !self
            .sequence_format_clone
            .format_regex
            .is_match(&self.raw_sequence.sequence)
        {
            self.raw_sequence.fix_constant_region(
                &self.sequence_format_clone.format_string,
                self.max_errors_clone.max_constant_errors(),
            )?;
        }
        // find the barcodes with the reges search
        let barcode_search = self
            .sequence_format_clone
            .format_regex
            .captures(&self.raw_sequence.sequence);

        // if the barcodes are found continue, else return None and record a constant region error
        if let Some(barcodes) = barcode_search {
            // Look for sample conversion
            let sample_barcode_match_option = barcodes.name("sample");
            let sample_barcode_option;
            // If there is a sample barcode included, save it to the sample barcode option variable
            // Otherwise record None and it will be ignored
            if let Some(sample_barcode) = sample_barcode_match_option {
                sample_barcode_option = Some(sample_barcode.as_str().to_string())
            } else {
                sample_barcode_option = None
            }

            // look for random barcode
            let random_barcode_match_option = barcodes.name("random");
            let random_barcode_option;

            // If there is a random barcode included, save it to the random barcode option variable
            // Otherwise record None and it will be ignored
            if let Some(random_barcode) = random_barcode_match_option {
                random_barcode_option = Some(random_barcode.as_str().to_string())
            } else {
                random_barcode_option = None
            }

            // If sample barcode is in the sample conversion file, convert. Otherwise try and fix the error
            let sample_name_option;
            if let Some(sample_barcode) = sample_barcode_option {
                if let Some(sample) = self.samples_clone.as_ref().unwrap().get(&sample_barcode) {
                    sample_name_option = Some(sample.to_string())
                } else {
                    let sample_seq_option = fix_error(
                        &sample_barcode,
                        self.sample_seqs.as_ref().unwrap(),
                        self.max_errors_clone.max_sample_errors(),
                    )?;
                    if let Some(sample_seq_new) = sample_seq_option {
                        sample_name_option = Some(
                            self.samples_clone
                                .as_ref()
                                .unwrap()
                                .get(&sample_seq_new)
                                .unwrap()
                                .to_string(),
                        );
                    } else {
                        // if sample barcode cannot be fixed, record the error
                        self.sequence_errors_clone
                            .lock()
                            .unwrap()
                            .sample_barcode_error();
                        self.seq_match_result_option = None;
                        return Ok(());
                    }
                }
            } else {
                // If there is not a sample barcode, record the sample name as Unknown sample
                sample_name_option = Some("Unknown_sample_name".to_string())
            }
            // If the sample barcode -> ID is found, add the building block suquences to the result string.  Otherwise, return None
            if let Some(sample_name) = sample_name_option {
                // if there is a building block conversion file used, convert
                if let Some(barcodes_seqs) = &self.barcodes_seqs_option {
                    let mut barcodes_string = String::new();
                    // fore each building block, convert and add as comma separated to a key text for results hashmap
                    for x in 0..self.sequence_format_clone.barcode_num {
                        let mut barcodes_seq =
                            barcodes[format!("barcode{}", x + 1).as_str()].to_string();
                        // If the building block sequence does not exists within the conversion file, try and fix
                        // If it cannnot fix the sequence, add to barcode_error
                        if !barcodes_seqs[x].contains(&barcodes_seq) {
                            let barcodes_seq_fix_option = fix_error(
                                &barcodes_seq,
                                &barcodes_seqs[x],
                                self.max_errors_clone.max_barcode_errors(),
                            )?;
                            if let Some(barcodes_seq_fix) = barcodes_seq_fix_option {
                                barcodes_seq = barcodes_seq_fix
                            } else {
                                self.sequence_errors_clone.lock().unwrap().barcode_error();
                                self.seq_match_result_option = None;
                                return Ok(());
                            }
                        }
                        // If it is the start just push the barcodes DNA barcode, otherwise a comma and the barcodes DNA barcode
                        // This is converted while writing to disk in case the memory size of the conversion would be too large
                        // ie, 6 DNA nucleotides takes up less memory than a long SMIILES string
                        if x != 0 {
                            barcodes_string.push(',');
                        }
                        barcodes_string.push_str(&barcodes_seq);
                    }
                    // if all goes well with barcodes_conversion and sample conversion, return Some
                    self.seq_match_result_option = Some(SequenceMatchResult::new(
                        sample_name,
                        barcodes_string,
                        random_barcode_option,
                    ));
                    Ok(())
                } else {
                    // If there is not a building block conversion file, do not try and fix the barcode errors. Push the raw DNA barcode seqeunces
                    let mut barcodes_string = barcodes["barcode1"].to_string();
                    for x in 1..self.sequence_format_clone.barcode_num {
                        let barcodes_num = format!("barcode{}", x + 1);
                        barcodes_string.push(',');
                        barcodes_string.push_str(&barcodes[barcodes_num.as_str()]);
                    }
                    self.seq_match_result_option = Some(SequenceMatchResult::new(
                        sample_name,
                        barcodes_string,
                        random_barcode_option,
                    ));
                    Ok(())
                }
            } else {
                // If the sample barcode was not found record the error and return None
                self.sequence_errors_clone
                    .lock()
                    .unwrap()
                    .sample_barcode_error();
                self.seq_match_result_option = None;
                Ok(())
            }
        } else {
            // If the constant region was not found, record the error and return None
            self.sequence_errors_clone
                .lock()
                .unwrap()
                .constant_region_error();
            self.seq_match_result_option = None;
            Ok(())
        }
    }
}

struct RawSequence {
    sequence: String,
}

impl RawSequence {
    pub fn new(sequence: String) -> RawSequence {
        RawSequence { sequence }
    }

    pub fn insert_barcodes_constant_region(&mut self, format_string: &str, best_sequence: String) {
        // Start a new string to push to
        let mut fixed_sequence = String::new();
        // Push the correct constant region nucleotides.  If the constant string has an N, push the nucleotides from the original
        // sequence corresponding to the barcodes
        for (old_char, new_char) in best_sequence.chars().zip(format_string.chars()) {
            if new_char == 'N' {
                fixed_sequence.push_str(&old_char.to_string());
            } else {
                fixed_sequence.push_str(&new_char.to_string());
            }
        }
        self.sequence = fixed_sequence
    }

    pub fn fix_constant_region(
        &mut self,
        format_string: &str,
        max_constant_errors: usize,
    ) -> Result<(), Box<dyn Error>> {
        // Find the region of the sequence that best matches the constant region.  This is doen by iterating through the sequence
        // Get the length difference between what was sequenced and the barcode region with constant regions
        // This is to stop the iteration in the next step
        let length_diff = self.sequence.len() - format_string.len();

        // Create a vector of sequences the length of the constant region + barcodes to check for where the best match is located
        let mut possible_seqs = Vec::new();
        for index in 0..length_diff {
            let possible_seq = self
                .sequence
                .chars()
                .skip(index) // skip to where the current index is and take the next amount equal to the length of the constant region + barcodes
                .take(format_string.len())
                .collect::<String>();
            // Add the new sequence to the vector of all possible matches
            possible_seqs.push(possible_seq);
        }
        // Find the closest match within what was sequenced to the constant region
        let best_sequence_option = fix_error(format_string, &possible_seqs, max_constant_errors)?;

        if let Some(best_sequence) = best_sequence_option {
            self.insert_barcodes_constant_region(format_string, best_sequence);
            Ok(())
        } else {
            self.sequence = "".to_string();
            Ok(())
        }
    }
}

struct SequenceMatchResult {
    sample_name: String,
    barcode_string: String,
    random_barcode_option: Option<String>,
}

impl SequenceMatchResult {
    pub fn new(
        sample_name: String,
        barcode_string: String,
        random_barcode_option: Option<String>,
    ) -> SequenceMatchResult {
        SequenceMatchResult {
            sample_name,
            barcode_string,
            random_barcode_option,
        }
    }
}

/// Fix an error in a sequence by comparing it to all possible sequences.  If no sequence matches with fewer or equal to the number of mismatches 'None' is returned.
/// 'None' is also returned if two or more sequences are best matches,
///
/// # Example
///
/// ```
/// use barcode::parse_sequences::fix_error;
///
/// let barcode = "AGTAG".to_string();
///
/// let possible_barcodes_one_match = vec!["AGCAG".to_string(), "ACAAG".to_string(), "AGCAA".to_string()]; // only the first has a single mismatch
/// let possible_barcodes_two_match = vec!["AGCAG".to_string(), "AGAAG".to_string(), "AGCAA".to_string()]; // first and second have a single mismatch
///
/// let max_mismatches = barcode.chars().count() / 5; // allow up to 20% mismatches
///
/// let fixed_error_one = fix_error(barcode, &possible_barcodes_one_match, max_mismatches).unwrap();
/// let fixed_error_two = fix_error(barcode, &possible_barcodes_two_match, max_mismatches).unwrap();
///
/// assert_eq!(fixed_error_one, Some("AGCAG".to_string()));
/// assert_eq!(fixed_error_two, None);
/// ```
pub fn fix_error(
    mismatch_seq: &str,
    possible_seqs: &[String],
    mismatches: usize,
) -> Result<Option<String>, Box<dyn Error>> {
    let mut best_match = None; // start the best match with None
    let mut best_mismatch_count = mismatches + 1; // Add 1 and start the best.  This allows a match with the same mismatches as required
    let mut keep = true; // An initiated variable to check if there is more than one best match

    // Iterate through possible matches
    for true_seq in possible_seqs {
        // Initiate the number of mismatches for the current iterated possible sequece match
        let mut mismatches = 0;

        // Iterate through the nucleotides of the possible match and the sequence to be fixed finding how many mismatches
        // If the mismatches exceed the current best mismatched, end this early
        for (possible_char, current_char) in true_seq.chars().zip(mismatch_seq.chars()) {
            if possible_char != current_char && current_char != 'N' && possible_char != 'N' {
                mismatches += 1;
            }
            if mismatches > best_mismatch_count {
                break;
            }
        }
        // If there are more than one best match, don't keep
        if mismatches == best_mismatch_count {
            keep = false
        }
        // If this is the best match, keep and reset best mismatches to this value
        if mismatches < best_mismatch_count {
            keep = true;
            best_mismatch_count = mismatches;
            best_match = Some(true_seq.to_string());
        }
    }
    // If there is one best match and it is some, return it.  Otherwise return None
    if keep && best_match.is_some() {
        Ok(best_match)
    } else {
        Ok(None)
    }
}
