use regex::Regex;
use std::{
    collections::HashMap,
    error::Error,
    sync::{Arc, Mutex},
};

/// Parses the sequence into its barcodes and converts the sample barcode to ID.
/// Setup with thread safe variables to allow multithreading.
/// There is some error correction.  For the constant region, up to 20% error is allowed.
/// It finds the best match among possible candidates.  If there are two best matches, then it leads to an error recorded and no count
pub fn parse(
    seq_clone: Arc<Mutex<Vec<String>>>,
    finished_clone: Arc<Mutex<bool>>,
    regex_string: String,
    constant_clone: String,
    results_clone: Arc<Mutex<HashMap<String, HashMap<String, u32>>>>,
    random_barcodes_clone: Arc<Mutex<HashMap<String, HashMap<String, Vec<String>>>>>,
    samples_clone: Option<HashMap<String, String>>,
    bb_clone: Option<HashMap<usize, HashMap<String, String>>>,
    sequence_errors_clone: Arc<Mutex<super::del_info::SequenceErrors>>,
) -> Result<(), Box<dyn Error>> {
    // Create a new regex search that has captures for each barcode
    let format_search = Regex::new(&regex_string)?;

    // Get the number of barcodes for later use
    let building_block_num = regex_string.matches("bb").count();

    // Get a vec of all possible sample barcodes for error correction
    let sample_seqs: Option<Vec<String>>;
    if let Some(ref samples) = samples_clone {
        sample_seqs = Some(samples.keys().map(|key| key.to_string()).collect());
    } else {
        sample_seqs = None
    }

    // Get a vec of all possible sample barcodes for error correction
    let bb_seqs_option: Option<Vec<Vec<String>>>;
    if let Some(ref bb) = bb_clone {
        let mut bb_vec = Vec::new();
        let mut bb_keys = bb.keys().collect::<Vec<&usize>>();
        bb_keys.sort();
        for key in bb_keys {
            let bb_data = bb.get(key).unwrap();
            let bb_barcodes = bb_data
                .keys()
                .map(|key| key.to_string())
                .collect::<Vec<String>>();
            bb_vec.push(bb_barcodes);
        }
        bb_seqs_option = Some(bb_vec);
    } else {
        bb_seqs_option = None
    }

    // Loop until there are no sequences left to parse.  These are fed into seq_clone vec by the reader thread
    loop {
        // If there are no sequences in seq_clone, pause, unless the reader thread is finished
        while seq_clone.lock().unwrap().is_empty() {
            if *finished_clone.lock().unwrap() {
                break;
            }
        }
        // If thre are no sequences and the reader thread is finished, break out of the loop
        if seq_clone.lock().unwrap().is_empty() && *finished_clone.lock().unwrap() {
            break;
        }

        // Pop off the last sequence from the seq_clone vec
        let seq = seq_clone.lock().unwrap().pop();
        // If this returns a Some value proceed
        if let Some(sequence) = seq {
            // get the result string match from the sequence.  This is a String with commas between sample_ID and builbding block sequences.
            // This is a convenience for later writing to a file since it will already be comma separated
            let result_tuple_option = match_seq(
                &sequence,
                &format_search,
                &samples_clone,
                &sequence_errors_clone,
                building_block_num,
                &sample_seqs,
                &bb_seqs_option,
                &constant_clone,
            )?;
            // If the constant region matched proceed to add to results, otherwise try and fix the constant region
            if let Some((sample_name, bb_string, random_barcode_option)) = result_tuple_option {
                // Alwasy add value unless random barcode is included and it has already been found for the sample and building blocks
                let mut add_value = true;
                // If there is a random barcode included
                if let Some(random_barcode) = random_barcode_option {
                    // Unlock the random_hashmap
                    let mut random_hashmap = random_barcodes_clone.lock().unwrap();

                    // If it does not already have the sample name, insert the sample name -> building_block -> random_barcodes
                    if !random_hashmap.contains_key(&sample_name) {
                        let mut intermediate_hashmap = HashMap::new();
                        let intermediate_vec = vec![random_barcode];
                        intermediate_hashmap.insert(bb_string.clone(), intermediate_vec);
                        random_hashmap.insert(sample_name.clone(), intermediate_hashmap);
                    } else {
                        let bb_hashmap = random_hashmap.get_mut(&sample_name).unwrap();
                        // If the random hashmap does not have the building blocks yet, insert building_block -> random_barcodes
                        if !bb_hashmap.contains_key(&bb_string) {
                            let intermediate_vec = vec![random_barcode];
                            bb_hashmap.insert(bb_string.clone(), intermediate_vec);
                        } else {
                            let random_vec = bb_hashmap.get_mut(&bb_string).unwrap();
                            // else check if the random barcode already used for the sample_name and building_blocks
                            // if the random barcode is already in the vector, change add_value to false
                            // otherqise add the random barcode to the random_barcodes vector
                            if random_vec.contains(&random_barcode) {
                                add_value = false
                            } else {
                                random_vec.push(random_barcode)
                            }
                        }
                    }
                }
                // Add 1 count to the results hashmap
                if add_value {
                    let mut results_hashmap = results_clone.lock().unwrap();

                    // If results hashmap does not already contain the sample_name, insert the sanmle_name -> barcodes -> 0
                    if !results_hashmap.contains_key(&sample_name) {
                        let mut intermediate_hashmap = HashMap::new();
                        intermediate_hashmap.insert(bb_string.clone(), 0);
                        results_hashmap.insert(sample_name.clone(), intermediate_hashmap);
                    }

                    // Insert 0 if the barcodes is not within the sample_name -> barcodes
                    // Then add one regardless
                    *results_hashmap
                        .get_mut(&sample_name)
                        .unwrap()
                        .entry(bb_string)
                        .or_insert(0) += 1;
                }
            }
        }
    }
    Ok(())
}

/// Does a regex search and captures the barcodes.  Converts the sample barcode to ID.  Returns a String with commas between Sample_ID and
/// building block barcodes.  This is used as a key within the results vector, where the value can be used as the count
fn match_seq(
    sequence: &String,
    format_search: &Regex,
    samples_clone: &Option<HashMap<String, String>>,
    sequence_errors_clone: &Arc<Mutex<super::del_info::SequenceErrors>>,
    building_block_num: usize,
    sample_seqs: &Option<Vec<String>>,
    bb_seqs_option: &Option<Vec<Vec<String>>>,
    constant_clone: &String,
) -> Result<Option<(String, String, Option<String>)>, Box<dyn Error>> {
    // find the barcodes with the reges search
    let mut barcode_search = format_search.captures(&sequence);

    // If the barcode search results in None, try and fix the constant regions
    let fixed_sequence;
    if barcode_search.is_none() {
        let fixed_sequence_option = fix_constant_region(&sequence, constant_clone)?;
        // If a suitable fix for the constant region was found, recreate the barcode search,
        // otherwise record the constant region error and return None
        if fixed_sequence_option.is_some() {
            fixed_sequence = fixed_sequence_option.unwrap();
            barcode_search = format_search.captures(&fixed_sequence);
        } else {
            sequence_errors_clone
                .lock()
                .unwrap()
                .constant_region_error();
            return Ok(None);
        }
    }

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
            if let Some(sample) = samples_clone.as_ref().unwrap().get(&sample_barcode) {
                sample_name_option = Some(sample.to_string())
            } else {
                let sample_seq_option =
                    fix_error(&sample_barcode, &sample_seqs.as_ref().unwrap(), 3)?;
                if let Some(sample_seq_new) = sample_seq_option {
                    sample_name_option = Some(
                        samples_clone
                            .as_ref()
                            .unwrap()
                            .get(&sample_seq_new)
                            .unwrap()
                            .to_string(),
                    );
                } else {
                    // if sample barcode cannot be fixed, record the error
                    sequence_errors_clone.lock().unwrap().sample_barcode_error();
                    return Ok(None);
                }
            }
        } else {
            // If there is not a sample barcode, record the sample name as Unknown sample
            sample_name_option = Some("Unknown_sample_name".to_string())
        }
        // If the sample barcode -> ID is found, add the building block suquences to the result string.  Otherwise, return None
        if let Some(sample_name) = sample_name_option {
            // if there is a building block conversion file used, convert
            if let Some(bb_seqs) = bb_seqs_option {
                let mut bb_string = String::new();
                // fore each building block, convert and add as comma separated to a key text for results hashmap
                for x in 0..building_block_num {
                    let mut bb_seq = barcodes[format!("bb{}", x + 1).as_str()].to_string();
                    // If the building block sequence does not exists within the conversion file, try and fix
                    // If it cannnot fix the sequence, add to bb_barcode_error
                    if !bb_seqs[x].contains(&bb_seq) {
                        let bb_seq_fix_option =
                            fix_error(&bb_seq, &bb_seqs[x], bb_seq.chars().count() / 5)?;
                        if let Some(bb_seq_fix) = bb_seq_fix_option {
                            bb_seq = bb_seq_fix
                        } else {
                            sequence_errors_clone.lock().unwrap().bb_barcode_error();
                            return Ok(None);
                        }
                    }
                    // If it is the start just push the bb DNA barcode, otherwise a comma and the bb DNA barcode
                    // This is converted while writing to disk in case the memory size of the conversion would be too large
                    // ie, 6 DNA nucleotides takes up less memory than a long SMIILES string
                    if x == 0 {
                        bb_string.push_str(&bb_seq);
                    } else {
                        bb_string.push_str(",");
                        bb_string.push_str(&bb_seq);
                    }
                }
                // if all goes well with bb_conversion and sample conversion, return Some
                return Ok(Some((sample_name, bb_string, random_barcode_option)));
            } else {
                // If there is not a building block conversion file, do not try and fix the barcode errors. Push the raw DNA barcode seqeunces
                let mut bb_string = barcodes["bb1"].to_string();
                for x in 1..building_block_num {
                    let bb_num = format!("bb{}", x + 1);
                    bb_string.push_str(",");
                    bb_string.push_str(&barcodes[bb_num.as_str()]);
                }
                return Ok(Some((sample_name, bb_string, random_barcode_option)));
            }
        } else {
            // If the sample barcode was not found record the error and return None
            sequence_errors_clone.lock().unwrap().sample_barcode_error();
            Ok(None)
        }
    } else {
        // If the constant region was not found, record the error and return None
        sequence_errors_clone
            .lock()
            .unwrap()
            .constant_region_error();
        Ok(None)
    }
}

/// Fix an error in a sequence by comparing it to all possible sequences.  If no sequence matches with fewer or equal to the number of mismatches 'None' is returned.
/// 'None' is also returned if two or more sequences are best matches,
///
/// # Example
///
/// ```
/// use del::parse_sequences::fix_error;
///
/// let barcode = "AGTAG".to_string();
///
/// let possible_barcodes_one_match = vec!["AGCAG".to_string(), "ACAAG".to_string(), "AGCAA".to_string()]; // only the first has a single mismatch
/// let possible_barcodes_two_match = vec!["AGCAG".to_string(), "AGAAG".to_string(), "AGCAA".to_string()]; // first and second have a single mismatch
///
/// let max_mismatches = barcode.chars().count() / 5; // allow up to 20% mismatches
///
/// let fixed_error_one = fix_error(&barcode, &possible_barcodes_one_match, max_mismatches).unwrap();
/// let fixed_error_two = fix_error(&barcode, &possible_barcodes_two_match, max_mismatches).unwrap();
///
/// assert_eq!(fixed_error_one, Some("AGCAG".to_string()));
/// assert_eq!(fixed_error_two, None);
/// ```
pub fn fix_error(
    mismatch_seq: &String,
    possible_seqs: &Vec<String>,
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

fn fix_constant_region(
    sequence: &String,
    constant_clone: &String,
) -> Result<Option<String>, Box<dyn Error>> {
    // Find the region of the sequence that best matches the constant region.  This is doen by iterating through the sequence
    let errors_allowed = (constant_clone.len() - constant_clone.matches("N").count()) / 5; // errors allowed is the length of the constant region - the Ns / 5 or 20%
                                                                                           // Get the length difference between what was sequenced and the barcode region with constant regions
                                                                                           // This is to stop the iteration in the next step
    let length_diff = sequence.len() - constant_clone.len();

    // Create a vector of sequences the length of the constant region + barcodes to check for where the best match is located
    let mut possible_seqs = Vec::new();
    for index in 0..length_diff {
        let possible_seq = sequence
            .chars()
            .skip(index) // skip to where the current index is and take the next amount equal to the length of the constant region + barcodes
            .take(constant_clone.len())
            .collect::<String>();
        // Add the new sequence to the vector of all possible matches
        possible_seqs.push(possible_seq);
    }
    // Find the closest match within what was sequenced to the constant region
    fix_error(&constant_clone, &possible_seqs, errors_allowed)
}

/// Fixes the constant region of the sequence by flipping the barcodes into the constant region format string within the locations of the 'N's
/// and returning the resulted string
///
/// # Example
///
/// ```
/// use del::parse_sequences::insert_barcodes_constant_region;
///
/// let sequence = "AGTAGATCTGAGATAGACAGC".to_string();// A 'CG' sequencing error is in the middle when compared to the format
/// let sequence_format = "AGTAGNNNTGACGTANNNAGC".to_string();
/// let fixed_sequence = insert_barcodes_constant_region(sequence, &sequence_format); // flips in the 'CG' sequencing error
/// assert_eq!(fixed_sequence, "AGTAGATCTGACGTAGACAGC".to_string())
/// ```
pub fn insert_barcodes_constant_region(old_sequence: String, sequence_fix: &String) -> String {
    // Start a new string to push to
    let mut fixed_sequence = String::new();
    // Push the correct constant region nucleotides.  If the constant string has an N, push the nucleotides from the original
    // sequence corresponding to the barcodes
    for (old_char, new_char) in old_sequence.chars().zip(sequence_fix.chars()) {
        if new_char == 'N' {
            fixed_sequence.push_str(&old_char.to_string());
        } else {
            fixed_sequence.push_str(&new_char.to_string());
        }
    }
    return fixed_sequence;
}
