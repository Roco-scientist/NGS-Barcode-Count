pub mod parse_sequences;

use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader},
    sync::{Arc, Mutex},
};

pub fn read_fastq(fastq: String, seq_clone: Arc<Mutex<Vec<String>>>) -> Result<(), Box<dyn Error>> {
    let fastq_file =
        File::open(fastq).unwrap_or_else(|err| panic!("Error reading fastq file: {}", err));
    let mut line_num = 1;
    for line in BufReader::new(fastq_file).lines() {
        if line_num == 2 {
            if let Ok(sequence_data) = line {
                while seq_clone.lock().unwrap().len() >= 10000 {}
                seq_clone.lock().unwrap().insert(0, sequence_data);
            }
        }
        line_num += 1;
        if line_num == 5 {
            line_num = 1
        }
    }
    Ok(())
}
