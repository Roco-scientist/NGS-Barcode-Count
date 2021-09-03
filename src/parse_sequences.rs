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
    results_clone: Arc<Mutex<HashMap<String, u32>>>,
    samples_clone: HashMap<String, String>,
    thread: u8,
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
            let barcode_search = format_search.captures(&sequence);
            if let Some(barcodes) = barcode_search {
                let sample_seq = barcodes["sample"].to_string();
                let sample_name_option;
                if let Some(sample) = samples_clone.get(&sample_seq) {
                    sample_name_option = Some(sample.to_string())
                } else {
                    let sample_seq_option = fix_error(&sample_seq, &sample_seqs, 3)?;
                    if let Some(sample_seq_new) = sample_seq_option {
                        sample_name_option =
                            Some(samples_clone.get(&sample_seq_new).unwrap().to_string());
                    } else {
                        sample_name_option = None;
                    }
                }
                if let Some(mut result_string) = sample_name_option {
                    for x in 0..building_block_num {
                        let bb_num = format!("bb{}", x + 1);
                        result_string.push_str(",");
                        result_string.push_str(&barcodes[bb_num.as_str()]);
                    }
                    // println!("Thread {} - result: {}", thread, &result_string);

                    let mut results_hashmap = results_clone.lock().unwrap();

                    let count = results_hashmap.entry(result_string).or_insert(0);
                    *count += 1;
                }
            } else {
                // println!("Thread {} - barcodes not found: {}", thread, sequence);
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
            if possible_char != current_char {
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
            best_match = Some(true_seq);
        }
    }
    if keep && best_match.is_some() {
        Ok(Some(best_match.unwrap().to_string()))
    } else {
        Ok(None)
    }
}
