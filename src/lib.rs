pub mod del_info;
pub mod parse_sequences;

// use flate2::read::GzDecoder;
use std::{
    error::Error,
    fs::File,
    io::{BufRead, BufReader},
    sync::{Arc, Mutex},
};

pub fn read_fastq(fastq: String, seq_clone: Arc<Mutex<Vec<String>>>) -> Result<(), Box<dyn Error>> {
    let fastq_file = File::open(fastq.clone())?;
    let mut line_num = 1;
    if !fastq.ends_with("gz") {
        let mut total_reads = 0;
        for line in BufReader::new(fastq_file).lines() {
            if line_num == 2 {
                if let Ok(sequence_data) = line {
                    while seq_clone.lock().unwrap().len() >= 10000 {}
                    seq_clone.lock().unwrap().insert(0, sequence_data);
                }
                total_reads += 1;
                print!("Total sequences: {}\r", total_reads);
            }
            line_num += 1;
            if line_num == 5 {
                line_num = 1
            }
        }
    } else {
        // TODO make this work with gzip files
        // let mut buffer = [0; 1000000000];
        // fastq_file.read(&buffer)?;
        // let flate_line = GzDecoder::new(buffer);
        // seq_clone.lock().unwrap().insert(0, flate_line);
    }
    println!();
    Ok(())
}
