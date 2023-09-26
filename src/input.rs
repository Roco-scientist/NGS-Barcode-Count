use anyhow::{bail, Context, Result};
use num_format::{Locale, ToFormattedString};
use std::{
    collections::VecDeque,
    fmt,
    fs::File,
    io::{BufRead, BufReader, Write},
    sync::{
        atomic::{AtomicBool, AtomicU32, Ordering},
        Arc, Mutex,
    },
};
use bgzip::BGZFReader;

use crate::parse::RawSequenceRead;

/// Reads in the FASTQ file line by line, then pushes every 2 out of 4 lines, which corresponds to the sequence line, into a Vec that is passed to other threads
///
/// FASTQ format:
/// Line 1: Sequence ID
/// Line 2: DNA sequence
/// Line 3: +
/// Line 4: Quality score
pub fn read_fastq(
    fastq: String,
    seq_clone: Arc<Mutex<VecDeque<String>>>,
    exit_clone: Arc<AtomicBool>,
    total_reads_arc: Arc<AtomicU32>,
) -> Result<()> {

    // Create a fastq line reader which keeps track of line number, reads, and posts the sequence to the shared vector
    let mut fastq_line_reader = FastqLineReader::new(seq_clone, exit_clone);

    let fastq_file = File::open(&fastq).context(format!("Failed to open file: {}", fastq))?; // open file
    // If the file is not gzipped use BufReader to read in lines
    if !fastq.ends_with("fastq.gz") {
        // If the file does not end with fastq, return with an error
        if !fastq.ends_with("fastq") {
            bail!("This program only works with *.fastq files and *.fastq.gz files.  The latter is still experimental")
        }

        // go line by line
        let mut stdout = std::io::stdout();
        let mut lock = stdout.lock();
        for line_result in BufReader::new(fastq_file).lines() {
            let mut line =
                line_result.context(format!("Bufread could not read line for file: {}", fastq))?;
            line.push('\n');
            // post the line to the shared vector and keep track of the number of sequences etc
            fastq_line_reader.read(line);
            if fastq_line_reader.line_num == 4 {
                fastq_line_reader.post()?;
            }
            // Add to read count to print numnber of sequences read by this thread
            if fastq_line_reader.total_reads % 10000 == 0 {
                write!(lock, "{}", fastq_line_reader)?;
                stdout.flush()?;
            }
        }
    } else {
        println!("If this program stops reading before the expected number of sequencing reads, unzip the gzipped fastq and rerun.");
        println!();
        // stream in first by decoding with GzDecoder, the reading into buffer
        let mut reader = BGZFReader::new(fastq_file)?;

        let mut stdout = std::io::stdout();
        let mut lock = stdout.lock();
        let mut read_response = 10;
        // continue reading until there is a response of 0, which indicates the end of file.  This may be where some gzipped files abrupty end
        while read_response != 0 {
            let mut line = String::new();
            read_response = reader.read_line(&mut line)?;
            // post the line to the shared vector and keep track of the number of sequences etc
            fastq_line_reader.read(line);
            if fastq_line_reader.line_num == 4 {
                fastq_line_reader.post()?;
            }
            // Add to read count to print numnber of sequences read by this thread
            if fastq_line_reader.total_reads % 10000 == 0 {
                write!(lock, "{}", fastq_line_reader)?;
                stdout.flush()?;
            }
        }
    }
    // Display the final total read count
    print!("{}", fastq_line_reader);
    total_reads_arc.store(fastq_line_reader.total_reads, Ordering::Relaxed);
    println!();
    Ok(())
}

/// A struct with functions for keeping track of read information and to post sequence lines to the shared vector
struct FastqLineReader {
    test: bool,   // whether or not to test the fastq format. Only does this for the first read
    line_num: u8, // the current line number 1-4.  Resets back to 1
    total_reads: u32, // total sequences read within the fastq file
    raw_sequence_read_string: String,
    seq_clone: Arc<Mutex<VecDeque<String>>>, // the vector that is passed between threads which containst the sequences
    exit_clone: Arc<AtomicBool>, // a bool which is set to true when one of the other threads panic.  This is the prevent hanging and is used to exit this thread
}

impl FastqLineReader {
    /// Creates a new FastqLineReader struct
    pub fn new(seq_clone: Arc<Mutex<VecDeque<String>>>, exit_clone: Arc<AtomicBool>) -> Self {
        FastqLineReader {
            test: true,
            line_num: 0,
            total_reads: 0,
            raw_sequence_read_string: String::new(),
            seq_clone,
            exit_clone,
        }
    }

    /// Reads in the line and either passes to the vec or discards it, depending if it is a sequence line.  Also increments on line count, sequence count etc.
    pub fn read(&mut self, line: String) {
        // Pause if there are already 10000 sequences in the vec so memory is not overloaded
        while self.seq_clone.lock().unwrap().len() >= 10000 {
            // if threads have failed exit out of this thread
            if self.exit_clone.load(Ordering::Relaxed) {
                break;
            }
        }
        // increase line number and if it has passed line 4, reset to 1
        self.line_num += 1;
        if self.line_num == 5 {
            self.line_num = 1
        }
        if self.line_num == 1 {
            self.total_reads += 1;
            self.raw_sequence_read_string = line;
        } else {
            self.raw_sequence_read_string.push_str(&line);
        }
    }

    pub fn post(&mut self) -> Result<()> {
        self.raw_sequence_read_string.pop(); // removes the last \n
                                             // Insert the sequence into the vec.  This will be popped out by other threads
        if self.test {
            RawSequenceRead::unpack(self.raw_sequence_read_string.clone())?.check_fastq_format()?;
            self.test = false;
        }
        self.seq_clone
            .lock()
            .unwrap()
            .push_front(self.raw_sequence_read_string.clone());
        Ok(())
    }
}

impl fmt::Display for FastqLineReader {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Total sequences:             {}\r",
            self.total_reads.to_formatted_string(&Locale::en)
        )
    }
}
