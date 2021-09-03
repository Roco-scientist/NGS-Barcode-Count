use regex::Regex;
use std::{
    error::Error,
    sync::{Arc, Mutex},
};

pub fn parse(
    seq_clone: Arc<Mutex<Vec<String>>>,
    finished_clone: Arc<Mutex<bool>>,
    regex_string: String,
    thread: u8,
) -> Result<(), Box<dyn Error>> {
    let format_search = Regex::new(&regex_string)?;
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
                println!(
                    "Thread {}\tSample: {}\tBB1: {}\tBB2: {}\tBB3: {}",
                    &thread,
                    &barcodes["sample"],
                    &barcodes["bb1"],
                    &barcodes["bb2"],
                    &barcodes["bb3"],
                )
            } else {
                println!("Barcodes not found")
            }
        }
    }
    Ok(())
}
