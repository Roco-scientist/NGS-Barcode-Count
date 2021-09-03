use std::{
    error::Error,
    sync::{Arc, Mutex},
};

pub fn parse(
    seq_clone: Arc<Mutex<Vec<String>>>,
    finished_clone: Arc<Mutex<bool>>,
    thread: u8,
) -> Result<(), Box<dyn Error>> {
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
            // println!("Thread - {} - Sequence: {}", thread, sequence);
        }
    }
    Ok(())
}
