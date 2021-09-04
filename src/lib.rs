pub mod del_info;
pub mod parse_sequences;

// use flate2::read::GzDecoder;
use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{BufRead, BufReader, Write},
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
                if total_reads % 1000 == 0 {
                    print!("Total sequences: {}\r", total_reads);
                }
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

pub fn output_counts(
    output_dir: String,
    results: Arc<Mutex<HashMap<String, u32>>>,
    bb_num: usize,
) -> Result<(), Box<dyn Error>> {
    let results_hasmap = results.lock().unwrap();
    let mut output = File::create(format!("{}{}", output_dir, "Counts.csv"))?;
    let mut header = "Sample_ID".to_string();
    for num in 0..bb_num {
        header.push_str(&format!(",BB_{}", num + 1))
    }
    header.push_str(",Count\n");
    output.write_all(header.as_bytes())?;
    for (code, count) in results_hasmap.iter() {
        let row = format!("{},{}\n", code, count);
        output.write_all(row.as_bytes())?;
    }
    Ok(())
}
