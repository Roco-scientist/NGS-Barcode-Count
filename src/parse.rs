use anyhow::{anyhow, Result};
use regex::Captures;
use std::{
    collections::VecDeque,
    fmt,
    sync::{
        atomic::{AtomicBool, Ordering},
        Arc, Mutex,
    },
};

use crate::info::{MaxSeqErrors, Results, SequenceErrors, SequenceFormat};
use ahash::AHashSet;

pub struct SequenceParser {
    shared_mut_clone: SharedMutData,
    sequence_errors_clone: SequenceErrors,
    sequence_format_clone: SequenceFormat,
    max_errors_clone: MaxSeqErrors,
    sample_seqs: AHashSet<String>,
    counted_barcode_seqs: Vec<AHashSet<String>>,
    raw_sequence: RawSequenceRead,
    barcode_groups: Vec<String>,
    min_quality_score: f32,
}

impl SequenceParser {
    pub fn new(
        shared_mut_clone: SharedMutData,
        sequence_errors_clone: SequenceErrors,
        sequence_format_clone: SequenceFormat,
        max_errors_clone: MaxSeqErrors,
        sample_seqs: AHashSet<String>,
        counted_barcode_seqs: Vec<AHashSet<String>>,
        min_quality_score: f32,
    ) -> Self {
        let mut barcode_groups = Vec::new();
        for x in 0..sequence_format_clone.barcode_num {
            barcode_groups.push(format!("barcode{}", x + 1))
        }
        SequenceParser {
            shared_mut_clone,
            sequence_errors_clone,
            sequence_format_clone,
            max_errors_clone,
            sample_seqs,
            counted_barcode_seqs,
            raw_sequence: RawSequenceRead::new(),
            barcode_groups,
            min_quality_score,
        }
    }
    pub fn parse(&mut self) -> Result<()> {
        // Loop until there are no sequences left to parse.  These are fed into seq vec by the reader thread
        loop {
            if self.get_seqeunce()? {
                if let Some(seq_match_result) = self.match_seq()? {
                    let barcode_string = seq_match_result.barcode_string();
                    // If there is a random barcode included
                    let added = self.shared_mut_clone.results.lock().unwrap().add_count(
                        &seq_match_result.sample_barcode,
                        seq_match_result.random_barcode.as_ref(),
                        barcode_string,
                    );
                    if added {
                        self.sequence_errors_clone.correct_match()
                    } else {
                        self.sequence_errors_clone.duplicated();
                    }
                }
            } else if self.shared_mut_clone.finished.load(Ordering::Relaxed) {
                break;
            }
        }
        Ok(())
    }

    fn get_seqeunce(&mut self) -> Result<bool> {
        // Pop off the last sequence from the seq vec
        if let Some(new_raw_sequence) = self.shared_mut_clone.seq.lock().unwrap().pop_back() {
            self.raw_sequence = RawSequenceRead::unpack(new_raw_sequence)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Does a regex search and captures the barcodes.  Returns a struct of the results.  
    fn match_seq(&mut self) -> Result<Option<SequenceMatchResult>> {
        self.check_and_fix_consant_region();
        // if the barcodes are found continue, else return None and record a constant region error
        if let Some(barcodes) = self
            .sequence_format_clone
            .format_regex
            .captures(&self.raw_sequence.sequence)
        {
            // If there was a minimum set for quality, check each barcode's quality
            if self.min_quality_score > 0.0 {
                if let Some(format_match) = self
                    .sequence_format_clone
                    .format_regex
                    .find(&self.raw_sequence.sequence)
                {
                    let start = format_match.start();
                    if self.raw_sequence.low_quality(
                        self.min_quality_score,
                        &self.sequence_format_clone.regions_string,
                        start,
                    ) {
                        // If any are low qualty, add to the low quality count and return
                        self.sequence_errors_clone.low_quality_barcode();
                        return Ok(None);
                    }
                } else {
                    return Err(anyhow!(
                        "Regex find failed after regex captures was successful"
                    ));
                }
            }

            // Create a match results struct which tests the regex regions
            let match_results = SequenceMatchResult::new(
                barcodes,
                &self.barcode_groups,
                &self.counted_barcode_seqs,
                self.max_errors_clone.max_barcode_errors(),
                &self.sample_seqs,
                self.max_errors_clone.max_sample_errors(),
            );

            // If the sample barcode was not found, record the error and return none so that the algorithm stops for this sequence
            if match_results.sample_barcode_error {
                self.sequence_errors_clone.sample_barcode_error();
                return Ok(None);
            }
            // If any of the counted barcodes were not found, even with error handling, record the error and return none so that the algorithm stops for this sequence
            if match_results.counted_barcode_error {
                self.sequence_errors_clone.barcode_error();
                return Ok(None);
            }
            // If all went well, return the match results struct
            Ok(Some(match_results))
        } else {
            // If the constant region was not found, record the error and return None
            self.sequence_errors_clone.constant_region_error();
            Ok(None)
        }
    }

    /// Checks the constant region of the sequence then finds the best fix if it is not found.  Basically whether or not the regex search worked
    fn check_and_fix_consant_region(&mut self) {
        // If the regex search does not work, try to fix the constant region
        if !self
            .sequence_format_clone
            .format_regex
            .is_match(&self.raw_sequence.sequence)
        {
            self.raw_sequence.fix_constant_region(
                &self.sequence_format_clone.format_string,
                self.max_errors_clone.max_constant_errors(),
            );
        }
    }
}

pub struct SharedMutData {
    pub seq: Arc<Mutex<VecDeque<String>>>,
    pub finished: Arc<AtomicBool>,
    pub results: Arc<Mutex<Results>>,
}

impl SharedMutData {
    pub fn new(
        seq: Arc<Mutex<VecDeque<String>>>,
        finished: Arc<AtomicBool>,
        results: Arc<Mutex<Results>>,
    ) -> Self {
        SharedMutData {
            seq,
            finished,
            results,
        }
    }

    pub fn arc_clone(&self) -> SharedMutData {
        let seq = Arc::clone(&self.seq);
        let finished = Arc::clone(&self.finished);
        let results = Arc::clone(&self.results);
        SharedMutData {
            seq,
            finished,
            results,
        }
    }
}

/// A struct to hold the raw sequencing information and transform it if there are sequencing errors
#[derive(Clone)]
pub struct RawSequenceRead {
    description: String,     // line 1 of fastq
    pub sequence: String,    // line 2 of fastq
    add_description: String, // line 3 of fastq
    quality_values: String,  // line 4 of fastq
}

impl Default for RawSequenceRead {
    fn default() -> Self {
        Self::new()
    }
}

impl RawSequenceRead {
    pub fn new() -> Self {
        RawSequenceRead {
            description: String::new(),
            sequence: String::new(),
            add_description: String::new(),
            quality_values: String::new(),
        }
    }

    pub fn new_fill(
        line_1: String,
        line_2: String,
        line_3: String,
        line_4: String,
    ) -> RawSequenceRead {
        RawSequenceRead {
            description: line_1,
            sequence: line_2,
            add_description: line_3,
            quality_values: line_4,
        }
    }

    pub fn add_line(&mut self, line_num: u16, line: String) -> Result<()> {
        match line_num {
            1 => self.description = line,
            2 => self.sequence = line,
            3 => self.add_description = line,
            4 => self.quality_values = line,
            _ => {
                return Err(anyhow!(
                "Too many new lines found within fastq read\nCurrent read\n{}\nCurrent line: {}",
                self,
                line
            ))
            }
        }
        Ok(())
    }

    pub fn pack(&self) -> String {
        format!(
            "{}\n{}\n{}\n{}",
            self.description, self.sequence, self.add_description, self.quality_values
        )
    }

    pub fn unpack(raw_string: String) -> Result<Self> {
        let mut raw_sequence_read = RawSequenceRead::new();
        for (line_num, line) in raw_string.split('\n').enumerate() {
            let true_line = line_num as u16 + 1;
            raw_sequence_read.add_line(true_line, line.to_string())?
        }
        Ok(raw_sequence_read)
    }

    /// Replaces the 'N's in the sequencing format with the barcodes to fix any sequencing errrors that would cause the regex search not to work
    pub fn insert_barcodes_constant_region(&mut self, format_string: &str, best_sequence: String) {
        // Start a new string to push to
        let mut fixed_sequence = String::new();
        // Push the correct constant region nucleotides.  If the constant string has an N, push the nucleotides from the original
        // sequence corresponding to the barcodes
        for (old_char, new_char) in best_sequence.chars().zip(format_string.chars()) {
            if new_char == 'N' {
                fixed_sequence.push_str(&old_char.to_string());
            } else {
                fixed_sequence.push_str(&new_char.to_string());
            }
        }
        self.sequence = fixed_sequence
    }

    /// Fixes the constant region by finding the closest match within the full seqeuence that has fewer than the max errors allowed,
    /// then uses the format string to flip the barcodes into the 'N's and have a fixed constant region string
    pub fn fix_constant_region(&mut self, format_string: &str, max_constant_errors: u16) {
        // Find the region of the sequence that best matches the constant region.  This is doen by iterating through the sequence
        // Get the length difference between what was sequenced and the barcode region with constant regions
        // This is to stop the iteration in the next step
        let length_diff = self.sequence.len() - format_string.len();

        // Create a vector of sequences the length of the constant region + barcodes to check for where the best match is located
        let mut possible_seqs = Vec::new();
        for index in 0..length_diff {
            let possible_seq = self
                .sequence
                .chars()
                .skip(index) // skip to where the current index is and take the next amount equal to the length of the constant region + barcodes
                .take(format_string.len())
                .collect::<String>();
            // Add the new sequence to the vector of all possible matches
            possible_seqs.push(possible_seq);
        }
        // Find the closest match within what was sequenced to the constant region
        let best_sequence_option = fix_error(format_string, &possible_seqs, max_constant_errors);

        if let Some(best_sequence) = best_sequence_option {
            self.insert_barcodes_constant_region(format_string, best_sequence);
        } else {
            self.sequence = "".to_string();
        }
    }

    /// Each DNA base read score within FASTQ is the ascii number - 33.
    /// This returns the number scores associated with the ascii values
    ///
    /// Score    Error Probability
    /// 40       0.0001
    /// 30       0.001
    /// 20       0.01
    /// 10       0.1
    pub fn quality_scores(&self) -> Vec<u8> {
        self.quality_values
            .chars()
            .map(|ch| ch as u8 - 33)
            .collect::<Vec<u8>>()
    }

    /// Test for if any of the barcode average quality score falls below the min_average cutoff
    pub fn low_quality(
        &self,
        min_average: f32,
        barcode_indicator_string: &str,
        start: usize,
    ) -> bool {
        let mut scores = Vec::new(); // vec to hold the quality scores for each barcode
        let mut previous_type = '\0'; // setup previoius barcode inidator type for the first comparison

        for (score, seq_type) in self
            .quality_scores()
            .iter()
            .skip(start)
            .zip(barcode_indicator_string.chars())
        // Zip score and barcode indicator
        {
            // Check for change in barcode/sequence type to know to calculate average score and create a new start to scores
            if seq_type != previous_type {
                // If scores is empty.  This avoids if there is a transition from consant region to barcode.  Constant region is not calculated
                if !scores.is_empty() {
                    // Get the average quality score for the barcode and if it is below the max average return true for low_quality
                    let sum: f32 = scores.iter().sum();
                    let average_score: f32 = sum / scores.len() as f32;
                    if average_score < min_average {
                        return true;
                    }
                    // Start a new vec for the next barcode
                    scores = Vec::new();
                }
                // set previous indicator type to the current type to check next
                previous_type = seq_type;
                // if indicator is not for a constant region, create a new vec with the new score
                if seq_type != 'C' {
                    scores = vec![*score as f32];
                }
            } else {
                // If indicator type is not for a constant region, add the score to the vec
                if seq_type != 'C' {
                    scores.push(*score as f32);
                }
            }
        }
        // If no average scores cause a true return, then return low_quality as false
        false
    }

    pub fn check_fastq_format(&self) -> Result<()> {
        // Test to see if the first line is not a sequence and the second is a sequence, which is typical fastq format
        match test_sequence(&self.description) {
            LineType::Sequence => {
                println!("{}", self);
                return Err(anyhow!("The first line within the FASTQ contains DNA sequences.  Check the FASTQ format"));
            }
            LineType::Metadata => (),
        }
        match test_sequence(&self.sequence) {
            LineType::Sequence => (),
            LineType::Metadata => {
                println!("{}", self);
                return Err(anyhow!("The second line within the FASTQ file is not a sequence. Check the FASTQ format"));
            }
        }
        Ok(())
    }
}

impl fmt::Display for RawSequenceRead {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Line 1: {}\nLine 2: {}\nLine 3: {}\nLine 4: {}",
            self.description, self.sequence, self.add_description, self.quality_values
        )
    }
}

/// An enum of linetype to use with test_sequence.  Every line of a FASTQ should either be sequence or metadata for the sequence.
enum LineType {
    Sequence,
    Metadata,
}

/// Tests whether a line within the file String is a sequence by checking if over half of the line contains DNA neceotide letters
fn test_sequence(sequence: &str) -> LineType {
    let sequence_length = sequence.len(); // the the length of the line
    let adenines = sequence.matches('A').count(); // And the amount of each DNA nucleotide
    let guanines = sequence.matches('G').count();
    let cytosines = sequence.matches('C').count();
    let thymines = sequence.matches('T').count();
    let any = sequence.matches('N').count();
    let total_dna = adenines + guanines + cytosines + thymines + any;
    // Check if less than half of the line contains DNA nucleotides.  If so, return that the line is metadata.  Otherwise, a sequence
    if total_dna < sequence_length / 2 {
        return LineType::Metadata;
    }
    LineType::Sequence
}

/// A struct to hold the results of the regex search on the sequence along with perform the functions to fix and find
pub struct SequenceMatchResult {
    pub sample_barcode: String,
    pub counted_barcodes: Vec<String>,
    pub counted_barcode_error: bool,
    pub sample_barcode_error: bool,
    pub random_barcode: Option<String>,
}

impl SequenceMatchResult {
    pub fn new(
        barcodes: Captures, // The regex result on the sequence
        barcode_groups: &[String],
        counted_barcode_seqs: &[AHashSet<String>], // The vec of known counted barcode sequences in order to fix sequencing errors.  Will be empty if none are known or included
        counted_barcode_max_errors: &[u16], // The maximum errors allowed for each counted barcode
        sample_seqs: &AHashSet<String>, // A hashset of all known sample barcodes. Will be empty if none are known or included
        sample_seqs_max_errors: u16,    // Maximum allowed sample barcode sequencing errors
    ) -> SequenceMatchResult {
        // Check for sample barcode and start with setting error to false
        let mut sample_barcode_error = false;
        let sample_barcode;
        // If 'sample' is within the regex returned search continue with checking and fixing
        if let Some(sample_barcode_match) = barcodes.name("sample") {
            let sample_barcode_str = sample_barcode_match.as_str();
            if sample_seqs.is_empty() {
                sample_barcode = sample_barcode_str.to_string();
            } else {
                // If the sample barcode is known save it
                if sample_seqs.contains(sample_barcode_str) {
                    sample_barcode = sample_barcode_str.to_string();
                } else {
                    // Otherwise try and fix it.  If the fix returns none, then save the error and an empty string
                    let sample_barcode_fix_option =
                        fix_error(sample_barcode_str, sample_seqs, sample_seqs_max_errors);
                    if let Some(fixed_barcode) = sample_barcode_fix_option {
                        sample_barcode = fixed_barcode;
                    } else {
                        sample_barcode = String::new();
                        sample_barcode_error = true;
                    }
                }
            }
        } else {
            // If there was no sample, save an empty string which should not have any allocation
            sample_barcode = "barcode".to_string();
        }

        // Check the counted barcodes and start with setting the error to false
        let mut counted_barcode_error = false;
        // Create an empty vec to hold the barcodes
        let mut counted_barcodes = Vec::new();
        // Only continue if the sample barcode was found
        if !sample_barcode_error {
            // Iterate through the counted barcocdes.  Fix if they are not within the known barcodes
            for (index, barcode_group) in barcode_groups.iter().enumerate() {
                let mut counted_barcode =
                    barcodes.name(barcode_group).unwrap().as_str().to_string();
                // If a barcode conversion file was included and there are known barcodes, check for sequencing errors
                if !counted_barcode_seqs.is_empty() {
                    // If the barcode is not known, try and fix
                    if !counted_barcode_seqs[index].contains(&counted_barcode) {
                        let barcode_seq_fix_option = fix_error(
                            &counted_barcode,
                            &counted_barcode_seqs[index],
                            counted_barcode_max_errors[index],
                        );
                        if let Some(fixed_barcode) = barcode_seq_fix_option {
                            counted_barcode = fixed_barcode;
                        } else {
                            // If a fix was not found, return the error and stop going through more barcodes
                            counted_barcode_error = true;
                            break;
                        }
                    }
                }
                // If all is well, add the counted barcode to the vec
                counted_barcodes.push(counted_barcode);
            }
        }

        // Chceck for a random barcode
        let random_barcode;
        // If a random barcode exists, add it.  Otherwise set it to an empty string
        if let Some(random_barcode_match) = barcodes.name("random") {
            random_barcode = Some(random_barcode_match.as_str().to_string())
        } else {
            random_barcode = None
        }
        SequenceMatchResult {
            sample_barcode,
            counted_barcodes,
            counted_barcode_error,
            sample_barcode_error,
            random_barcode,
        }
    }

    /// Returns a comma separated counted barcodes string.  Perfect for CSV file writing
    pub fn barcode_string(&self) -> String {
        self.counted_barcodes.join(",")
    }
}

/// Fix an error in a sequence by comparing it to all possible sequences.  If no sequence matches with fewer or equal to the number of mismatches 'None' is returned.
/// 'None' is also returned if two or more sequences are best matches.  Will work with vec and hashset
///
/// # Example
///
/// ```
/// use barcode_count::parse::fix_error;
///
/// let barcode = "AGTAG";
///
/// let possible_barcodes_one_match: std::collections::AHashSet<String> = ["AGCAG".to_string(), "ACAAG".to_string(), "AGCAA".to_string()].iter().cloned().collect(); // only the first has a single mismatch
/// let possible_barcodes_two_match: std::collections::AHashSet<String> = ["AGCAG".to_string(), "AGAAG".to_string(), "AGCAA".to_string()].iter().cloned().collect(); // first and second have a single mismatch
///
/// let max_mismatches = barcode.chars().count() as u16 / 5; // allow up to 20% mismatches
///
/// let fixed_error_one = fix_error(barcode, &possible_barcodes_one_match, max_mismatches);
/// let fixed_error_two = fix_error(barcode, &possible_barcodes_two_match, max_mismatches);
///
/// assert_eq!(fixed_error_one, Some("AGCAG".to_string()));
/// assert_eq!(fixed_error_two, None);
/// ```
pub fn fix_error<'a, I>(mismatch_seq: &str, possible_seqs: I, mismatches: u16) -> Option<String>
where
    I: IntoIterator<Item = &'a String>,
{
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
            if possible_char != current_char && current_char != 'N' && possible_char != 'N' {
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
        best_match
    } else {
        None
    }
}
