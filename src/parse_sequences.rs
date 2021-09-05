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
    results_clone: Arc<Mutex<HashMap<String, u32>>>,
    samples_clone: HashMap<String, String>,
    sequence_errors_clone: Arc<Mutex<super::del_info::SequenceErrors>>,
) -> Result<(), Box<dyn Error>> {
    // Create a new regex search that has captures for each barcode
    let format_search = Regex::new(&regex_string)?;

    // Get the number of barcodes for later use
    let building_block_num = regex_string.matches("bb").count();

    // Get a vec of all possible sample barcodes for error correction
    let sample_seqs: Vec<String> = samples_clone.keys().map(|key| key.to_string()).collect();

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
            let result_string_option = match_seq(
                &sequence,
                &format_search,
                &samples_clone,
                &sequence_errors_clone,
                building_block_num,
                &sample_seqs,
            )?;
            // If the constant region matched proceed to add to results, otherwise try and fix the constant region
            if let Some(result_string) = result_string_option {
                let mut results_hashmap = results_clone.lock().unwrap();

                let count = results_hashmap.entry(result_string).or_insert(0);
                *count += 1;
            } else {
                // Find the region of the sequence that best matches the constant region.  This is doen by iterating through the sequence
                let errors_allowed =
                    (constant_clone.len() - constant_clone.matches("N").count()) / 5; // errors allowed is the length of the constant region - the Ns / 5 or 20%
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
                let fixed_constant_option =
                    fix_error(&constant_clone, &possible_seqs, errors_allowed as u8)?;

                // If there is a match proceed and count, otherwise record that there was a constant region mismatch
                if let Some(new_sequence) = fixed_constant_option {
                    // Flip all barcodes into the constant's 'N's for the fixed sequence
                    let fixed_sequence = fix_constant_region(new_sequence, &constant_clone);
                    let result_string_option = match_seq(
                        &fixed_sequence,
                        &format_search,
                        &samples_clone,
                        &sequence_errors_clone,
                        building_block_num,
                        &sample_seqs,
                    )?;
                    if let Some(result_string) = result_string_option {
                        let mut results_hashmap = results_clone.lock().unwrap();

                        let count = results_hashmap.entry(result_string).or_insert(0);
                        *count += 1;
                    }
                } else {
                    sequence_errors_clone
                        .lock()
                        .unwrap()
                        .constant_region_error();
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
    samples_clone: &HashMap<String, String>,
    sequence_errors_clone: &Arc<Mutex<super::del_info::SequenceErrors>>,
    building_block_num: usize,
    sample_seqs: &Vec<String>,
) -> Result<Option<String>, Box<dyn Error>> {
    // find the barcodes with the reges search
    let barcode_search = format_search.captures(&sequence);

    // if the barcodes are found continue, else return None
    if let Some(barcodes) = barcode_search {
        // Look for sample conversion
        let sample_seq = barcodes["sample"].to_string();
        // If sample barcode is in the sample conversion file, convert. Otherwise try and fix the error
        let sample_name_option;
        if let Some(sample) = samples_clone.get(&sample_seq) {
            sample_name_option = Some(sample.to_string())
        } else {
            let sample_seq_option = fix_error(&sample_seq, &sample_seqs, 3)?;
            if let Some(sample_seq_new) = sample_seq_option {
                sample_name_option = Some(samples_clone.get(&sample_seq_new).unwrap().to_string());
            } else {
                sequence_errors_clone.lock().unwrap().sample_barcode_error();
                return Ok(None);
            }
        }
        // If the sample barcode -> ID is found, add the building block suquences to the result string.  Otherwise, return None
        if let Some(mut result_string) = sample_name_option {
            for x in 0..building_block_num {
                let bb_num = format!("bb{}", x + 1);
                result_string.push_str(",");
                result_string.push_str(&barcodes[bb_num.as_str()]);
            }
            return Ok(Some(result_string));
        }
    }
    Ok(None)
}

/// Fix an error in a sequence by comparing it to all possible sequences.  If no sequence matches with fewer or equal to the number of mismatches 'None' is returned.
/// 'None' is also returned if two or more sequences are best matches,
fn fix_error(
    mismatch_seq: &String,
    possible_seqs: &Vec<String>,
    mismatches: u8,
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
            if possible_char != current_char && current_char != 'N' {
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

/// Fixes the constant region of the sequence by flipping the barcodes into the constant region format string within the locations of the 'N's
/// and returning the resulted string
fn fix_constant_region(old_sequence: String, sequence_fix: &String) -> String {
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
