use regex::Regex;
use std::{
    collections::HashMap,
    error::Error,
    sync::{Arc, Mutex},
};

pub fn parse(
    seq_clone: Arc<Mutex<Vec<String>>>,
    finished_clone: Arc<Mutex<bool>>,
    regex_string: String,
    constant_clone: String,
    results_clone: Arc<Mutex<HashMap<String, u32>>>,
    samples_clone: HashMap<String, String>,
    sequence_errors_clone: Arc<Mutex<super::del_info::SequenceErrors>>,
) -> Result<(), Box<dyn Error>> {
    let format_search = Regex::new(&regex_string)?;
    let building_block_num = regex_string.matches("bb").count();
    let sample_seqs: Vec<String> = samples_clone.keys().map(|key| key.to_string()).collect();
    loop {
        while seq_clone.lock().unwrap().is_empty() {
            if *finished_clone.lock().unwrap() {
                break;
            }
        }
        if seq_clone.lock().unwrap().is_empty() && *finished_clone.lock().unwrap() {
            break;
        }
        let seq = seq_clone.lock().unwrap().pop();
        if let Some(sequence) = seq {
            let result_string_option = match_seq(
                &sequence,
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
            } else {
                let barcode_length = constant_clone.len();
                let ns = constant_clone.matches("N").count();
                let errors_allowed = (barcode_length - ns) / 5;
                let sequence_length = sequence.len();
                let length_diff = sequence_length - barcode_length;
                let mut possible_seqs = Vec::new();
                for index in 0..length_diff {
                    let possible_seq = sequence
                        .chars()
                        .skip(index)
                        .take(barcode_length)
                        .collect::<String>();
                    possible_seqs.push(possible_seq);
                }
                let fixed_constant_option =
                    fix_error(&constant_clone, &possible_seqs, errors_allowed as u8)?;
                if let Some(new_sequence) = fixed_constant_option {
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

fn fix_error(
    mismatch_seq: &String,
    possible_seqs: &Vec<String>,
    mismatches: u8,
) -> Result<Option<String>, Box<dyn Error>> {
    let mut best_match = None;
    let mut best_mismatch_count = mismatches + 1;
    let mut keep = true;
    for true_seq in possible_seqs {
        let mut mismatches = 0;
        for (possible_char, current_char) in true_seq.chars().zip(mismatch_seq.chars()) {
            if possible_char != current_char && current_char != 'N' {
                mismatches += 1;
            }
            if mismatches > best_mismatch_count {
                break;
            }
        }
        if mismatches == best_mismatch_count {
            keep = false
        }
        if mismatches < best_mismatch_count {
            keep = true;
            best_mismatch_count = mismatches;
            best_match = Some(true_seq.to_string());
        }
    }
    if keep && best_match.is_some() {
        Ok(best_match)
    } else {
        Ok(None)
    }
}

fn fix_constant_region(old_sequence: String, sequence_fix: &String) -> String {
    let mut fixed_sequence = String::new();
    for (old_char, new_char) in old_sequence.chars().zip(sequence_fix.chars()) {
        if new_char == 'N' {
            fixed_sequence.push_str(&old_char.to_string());
        } else {
            fixed_sequence.push_str(&new_char.to_string());
        }
    }
    return fixed_sequence;
}

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
