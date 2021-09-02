use std::{
    error::Error,
    sync::{Arc, Mutex},
};

pub fn parse(
    seq_clone: Arc<Mutex<Vec<String>>>,
    finished_clone: Arc<Mutex<bool>>,
) -> Result<(), Box<dyn Error>> {
    loop {
        while seq_clone.lock().unwrap().len() == 0 {
            if *finished_clone.lock().unwrap() {
                break;
            }
        }
        if *finished_clone.lock().unwrap() {
            break;
        }
        let seq = seq_clone.lock().unwrap().pop();
        if let Some(sequence) = seq {
            println!("Sequence: {}", sequence);
        }
    }
    Ok(())
}
