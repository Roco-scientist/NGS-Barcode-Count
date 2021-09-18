use itertools::Itertools;
use regex::Regex;
use std::{collections::HashMap, error::Error, fs};

// Struct to keep track of sequencing errors and correct matches.  This is displayed at the end of the algorithm for QC measures
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SequenceErrors {
    // errors within the constant region
    constant_region: u64,
    // errors within the sample barcode
    sample_barcode: u64,
    // erors within the building block barcode
    bb_barcode: u64,
    // total matched
    matched: u64,
    // total random barcode duplicates
    duplicates: u64,
}

impl SequenceErrors {
    /// Create a new sequence error struct.  Starts with 0 errors in all regions, then is added to later.
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// ```
    pub fn new() -> SequenceErrors {
        SequenceErrors {
            constant_region: 0,
            sample_barcode: 0,
            bb_barcode: 0,
            matched: 0,
            duplicates: 0,
        }
    }

    /// Add one to constant region error
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.constant_region_error();
    /// ```
    pub fn constant_region_error(&mut self) {
        self.constant_region += 1;
    }

    /// Add one to sample barcode error
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.sample_barcode_error();
    /// ```
    pub fn sample_barcode_error(&mut self) {
        self.sample_barcode += 1;
    }

    /// Add one to building block barcode error
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.bb_barcode_error();
    /// ```
    pub fn bb_barcode_error(&mut self) {
        self.bb_barcode += 1;
    }

    /// Add one to correct match
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.correct_match();
    /// ```
    pub fn correct_match(&mut self) {
        self.matched += 1;
    }

    /// Add one to duplicates
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.duplicated();
    /// ```
    pub fn duplicated(&mut self) {
        self.duplicates += 1;
    }

    /// Print to stdout all sequencing error counts
    ///
    /// # Example
    /// ```
    /// use del::del_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.bb_barcode_error();
    /// sequence_errors.display();
    /// ```
    pub fn display(&mut self) {
        println!(
            "Correctly matched sequences: {}\nConstant Region Mismatches:  {}\nSample Barcode Mismatches:   {}\nBuilding Block Mismatches:   {}\nDuplicates:                  {}",
            self.matched, self.constant_region, self.sample_barcode, self.bb_barcode, self.duplicates
        )
    }
}

#[derive(Debug, Clone)]
pub struct SequenceFormat {
    pub format_string: String,
    // Not implemented yet
    pub format_string_multiple: Option<Vec<String>>,
    pub format_regex: Regex,
    pub regex_string: String,
    pub bb_num: usize,
    // Not implemented yet
    pub multiple: bool,
    format_data: String,
}

impl SequenceFormat {
    /// Creates a new SequenceFormat struct which holds the sequencing format information, such as, where the barcodes are located within the sequence
    pub fn new(format: String) -> Result<SequenceFormat, Box<dyn Error>> {
        // Read sequenc format file to string
        let format_data = fs::read_to_string(format)?
            .lines() // split into lines
            .filter(|line| !line.starts_with("#")) // remove any line that starts with '#'
            .collect::<String>(); // collect into a String

        let regex_string = build_regex_captures(&format_data)?; // Build the regex string from the input format file information
        let format_regex = Regex::new(&regex_string)?; // Convert the regex string to a Regex
        let format_string = build_format_string(&format_data)?; // Create the format string replacing 'N's where there is a barcode
        let bb_num = regex_string.matches("bb").count(); // Count the number of barcodes.  This is used later for retrieving barcodes etc.

        // Create and return the SequenceFormat struct
        Ok(SequenceFormat {
            format_string,
            format_string_multiple: None,
            format_regex,
            regex_string,
            bb_num,
            multiple: false,
            format_data,
        })
    }

    /// Displays the sequence format information with 'N's replacing all barcodes
    pub fn display_format(&self) -> () {
        println!("Format: {}", self.format_string)
    }

    /// Returns a Vec of the size of all building block barcodes within the DEL scheme
    pub fn bb_lengths(&self) -> Result<Vec<usize>, Box<dyn Error>> {
        let bb_search = Regex::new(r"(\{\d+\})")?; // Create a search that finds the '{#}'
        let digit_search = Regex::new(r"\d+")?; // Create a search that pulls out the number
        let mut bb_lengths = Vec::new(); // Create a Vec that will contain on the building block barcode lengths

        // For each building block barcode found in the format file string
        for group in bb_search.find_iter(&self.format_data) {
            let group_str = group.as_str();
            // Pull out the numeric value
            let digits = digit_search
                .find(group_str)
                .unwrap()
                .as_str()
                .parse::<usize>()?;
            // And add to the vector
            bb_lengths.push(digits)
        }
        Ok(bb_lengths)
    }

    /// Returns the sample barcode length found in the format file string
    pub fn sample_length_option(&self) -> Result<Option<usize>, Box<dyn Error>> {
        let sample_search = Regex::new(r"(\[\d+\])")?; // Create a search that finds the '[#]'
        let digit_search = Regex::new(r"\d+")?; // Create a search that pulls out the numeric value

        // If there is a sample barcode inluded in the format file, find the size.  If not, return None
        if let Some(sample_match) = sample_search.find(&self.format_data) {
            let sample_str = sample_match.as_str();
            // Get the numeric value of the sample barcode size
            let digits = digit_search
                .find(sample_str)
                .unwrap()
                .as_str()
                .parse::<usize>()?;
            return Ok(Some(digits));
        } else {
            return Ok(None);
        }
    }

    /// Returns the amount of nucleotides within the constant regions from the format file
    pub fn constant_region_length(&self) -> usize {
        // Get the full length of the format_string and subtract the amount of 'N's found to get the constant nucleotide count
        self.format_string.len() - self.format_string.matches("N").count()
    }
}

/// Builds the catpure groups from the file format
///
/// # Example
///
/// ```
/// use del::del_info::build_format_string;
/// let format_data = "[8]AGCTAGATC{6}TGGA{6}TGGA{6}TGATTGCGC(6)NNNNAT".to_string();
///
/// assert_eq!(build_format_string(&format_data).unwrap(),  "NNNNNNNNAGCTAGATCNNNNNNTGGANNNNNNTGGANNNNNNTGATTGCGCNNNNNNNNNNAT".to_string())
/// ```
pub fn build_format_string(format_data: &String) -> Result<String, Box<dyn Error>> {
    let digit_search = Regex::new(r"\d+")?;
    let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|N+|[ATGC]+")?;
    let mut final_format = String::new();
    for group in barcode_search.find_iter(format_data) {
        let group_str = group.as_str();
        let digits_option = digit_search.find(group_str);
        if let Some(digit) = digits_option {
            let digit_value = digit.as_str().parse::<usize>()?;
            for _ in 0..digit_value {
                final_format.push_str("N")
            }
        } else {
            final_format.push_str(group_str)
        }
    }
    Ok(final_format)
}

/// Builds the catpure groups from the file format
///
/// # Example
///
/// ```
/// use del::del_info::build_regex_captures;
/// let format_data = "[8]AGCTAGATC{6}TGGA{6}TGGA{6}TGATTGCGC(6)NNNNAT".to_string();
///
/// assert_eq!(build_regex_captures(&format_data).unwrap(),  "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6}).{4}AT".to_string())
/// ```
pub fn build_regex_captures(format_data: &String) -> Result<String, Box<dyn Error>> {
    let digit_search = Regex::new(r"\d+")?;
    let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|N+|[ATGC]+")?;
    // the previous does not bumber each barcode but names each caputre with bb#
    // The '#' needs to bre replaced with teh sequential number
    let mut final_format = String::new();
    let mut bb_num = 0;
    // Fore each character, if the character is #, replace with the sequential building block number
    for group in barcode_search.find_iter(format_data) {
        let group_str = group.as_str();
        let mut group_name_option = None;
        if group_str.contains("[") {
            group_name_option = Some("sample".to_string())
        } else {
            if group_str.contains("{") {
                bb_num += 1;
                group_name_option = Some(format!("bb{}", bb_num));
            } else {
                if group_str.contains("(") {
                    group_name_option = Some("random".to_string());
                }
            }
        }
        if let Some(group_name) = group_name_option {
            let digits = digit_search
                .captures(&group_str)
                .unwrap()
                .get(0)
                .unwrap()
                .as_str()
                .parse::<usize>()
                .unwrap();
            let mut capture_group = format!("(?P<{}>.", group_name);
            capture_group.push('{');
            capture_group.push_str(&digits.to_string());
            capture_group.push_str("})");
            final_format.push_str(&capture_group);
        } else {
            if group_str.contains("N") {
                let num_of_ns = group_str.matches("N").count();
                let mut n_group = ".{".to_string();
                n_group.push_str(&num_of_ns.to_string());
                n_group.push('}');
                final_format.push_str(&n_group);
            } else {
                final_format.push_str(&group_str.to_uppercase())
            }
        }
    }
    return Ok(final_format);
}

/// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
/// and the second needs to be the ID
pub fn sample_barcode_file_conversion(
    barcode_path: &String,
) -> Result<HashMap<String, String>, Box<dyn Error>> {
    // read in the sample barcode file
    let barcode_data: HashMap<String, String> = fs::read_to_string(barcode_path)?
        .lines() // split the lines
        .skip(1) // skip the first line which should be the header
        .map(|line| {
            line.split(",")
                .take(2) // take only the first two values, or columns
                .map(|value| value.to_string())
                .collect_tuple()
                .unwrap_or(("".to_string(), "".to_string()))
        }) // comma split the line into a tuple with the first being the key and the last the value
        .collect::<Vec<(String, String)>>()
        .iter()
        .cloned()
        .collect(); // collect into a hashmap
    Ok(barcode_data)
}

/// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
/// the second needs to be the ID, and the third needs to be the building block number
///
/// # Panics
///
/// This panics if the third column of the barcode conversion file does not contain integers.  Also
/// panics if not all integers for barcode numbers is within this columns
pub fn bb_barcode_file_conversion(
    barcode_path: &String,
    bb_num: usize,
) -> Result<HashMap<usize, HashMap<String, String>>, Box<dyn Error>> {
    // read in the sample barcode file
    let barcode_vecs = fs::read_to_string(barcode_path)?
        .lines() // split the lines
        .skip(1) // skip the first line which should be the header
        .map(|line| {
            line.split(",")
                .take(3) // take only the first three values, or columns
                .map(|value| value.to_string())
                .collect_tuple()
                .unwrap_or(("".to_string(), "".to_string(), "".to_string()))
        }) // comma split the line into a tuple with the first being the key and the last the value
        .collect::<Vec<(String, String, String)>>();
    let mut barcode_data = HashMap::new();
    let mut barcode_num_contained = Vec::new();
    for (barcode, id, bb_num) in barcode_vecs {
        let bb_num_usize = bb_num.parse::<usize>().unwrap_or_else(|err| {
            panic!("Third column of building block barcode file contains something other than an integer: {}\nError: {}", bb_num, err)
        });
        if !barcode_num_contained.contains(&bb_num_usize) {
            barcode_num_contained.push(bb_num_usize)
        };
        if !barcode_data.contains_key(&bb_num_usize) {
            let mut intermediate_hash = HashMap::new();
            intermediate_hash.insert(barcode, id);
            barcode_data.insert(bb_num_usize, intermediate_hash);
        } else {
            barcode_data
                .get_mut(&bb_num_usize)
                .unwrap()
                .insert(barcode, id);
        }
    }
    let mut missing_barcode_num = Vec::new();
    for x in 0..bb_num {
        let actual_num = x + 1;
        if !barcode_num_contained.contains(&actual_num) {
            missing_barcode_num.push(actual_num)
        }
    }
    if !missing_barcode_num.is_empty() {
        panic!("Building block barcode conversion file missing barcode numers {:?} in the third column", missing_barcode_num)
    }
    Ok(barcode_data)
}

// Struct of how many sequencing errrors are allowed
#[derive(Debug, Clone, PartialEq)]
pub struct MaxSeqErrors {
    // errors within the constant region
    constant_region: usize,
    constant_region_size: usize,
    // errors within the sample barcode
    sample_barcode: usize,
    sample_size: usize,
    // erors within the building block barcode
    bb_barcode: usize,
    bb_sizes: Vec<usize>,
}

impl MaxSeqErrors {
    /// Create a new sequence error struct
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let bb_errors_option = None;
    /// let bb_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// ```
    pub fn new(
        sample_errors_option: Option<usize>,
        sample_barcode_size_option: Option<usize>,
        bb_errors_option: Option<usize>,
        bb_sizes: Vec<usize>,
        constant_errors_option: Option<usize>,
        constant_region_size: usize,
    ) -> Result<MaxSeqErrors, Box<dyn Error>> {
        let max_sample_errors;
        // start with a sample size of 0 in case there is no sample barcode.  If there is then mutate
        let mut sample_size = 0;
        // If sample barcode was included, calculate the maximum error, otherwise set error to 0
        if let Some(sample_size_actual) = sample_barcode_size_option {
            sample_size = sample_size_actual;
            // if there was sample errors input from arguments, use that, otherwise calculate 20% for max errors
            if let Some(sample_errors) = sample_errors_option {
                max_sample_errors = sample_errors
            } else {
                max_sample_errors = sample_size_actual / 5;
            }
        } else {
            max_sample_errors = 0;
        }

        let max_bb_errors;
        // If max error was set by input arguments, use that value, otherwise calculate 20% of barcode size for max error
        if let Some(bb_errors) = bb_errors_option {
            max_bb_errors = bb_errors
        } else {
            max_bb_errors = bb_sizes.iter().max().unwrap() / 5;
        }

        let max_constant_errors;
        // If max error was set by input arguments, use that value, otherwise calculate 20% of barcode size for max error
        if let Some(constant_errors) = constant_errors_option {
            max_constant_errors = constant_errors
        } else {
            max_constant_errors = constant_region_size / 5;
            // errors allowed is the length of the constant region - the Ns / 5 or 20%
        }

        Ok(MaxSeqErrors {
            constant_region: max_constant_errors,
            constant_region_size,
            sample_barcode: max_sample_errors,
            sample_size,
            bb_barcode: max_bb_errors,
            bb_sizes,
        })
    }

    /// Returns the maximum allowed constant region errors
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let bb_errors_option = None;
    /// let bb_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 6);
    /// let bb_sizes = vec![8,8,8];
    /// let constant_errors_option = Some(3);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 3);
    /// ```
    pub fn max_constant_errors(&mut self) -> usize {
        self.constant_region
    }

    /// Returns the maximum allowed sample barcode errors
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let bb_errors_option = None;
    /// let bb_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 2);
    /// let bb_sizes = vec![8,8,8];
    /// let sample_errors_option = Some(3);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 3);
    /// ```
    pub fn max_sample_errors(&mut self) -> usize {
        self.sample_barcode
    }

    /// Returns the maximum allowed building block barcode errors
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let bb_errors_option = None;
    /// let bb_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_bb_errors(), 1);
    /// let bb_sizes = vec![8,8,8];
    /// let bb_errors_option = Some(2);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_bb_errors(), 2);
    /// ```
    pub fn max_bb_errors(&mut self) -> usize {
        self.bb_barcode
    }

    /// Print to stdout all maximum sequencing errors
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let bb_errors_option = None;
    /// let bb_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, bb_errors_option, bb_sizes, constant_errors_option, constant_region_size).unwrap();
    /// max_sequence_errors.display();
    /// ```
    pub fn display(&mut self) {
        println!(
            "
            \n########## Barcode Info ###################################\n\
            Constant region size: {}\n\
            Maximum constant region mismatches allowed per sequence: {}\n\
            Sample barcode size: {}\n\
            Maximum sample barcode mismatches allowed per sequence: {}\n\
            Building block sizes: {:?}\n\
            Maximum building block mismatches allowed per barcode: {}\n",
            self.constant_region_size,
            self.constant_region,
            self.sample_size,
            self.sample_barcode,
            self.bb_sizes,
            self.bb_barcode
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_sequence_errors_test() {
        let sample_errors_option = None;
        let sample_barcode_size_option = Some(10);
        let bb_errors_option = None;
        let bb_sizes = vec![8, 8, 8];
        let constant_errors_option = None;
        let constant_region_size = 30;
        let max_sequence_errors = MaxSeqErrors::new(
            sample_errors_option,
            sample_barcode_size_option,
            bb_errors_option,
            bb_sizes,
            constant_errors_option,
            constant_region_size,
        )
        .unwrap();
        assert_eq!(
            max_sequence_errors,
            MaxSeqErrors {
                constant_region: 6,
                constant_region_size: 30,
                sample_barcode: 2,
                sample_size: 10,
                bb_barcode: 1,
                bb_sizes: vec![8, 8, 8],
            }
        )
    }

    #[test]
    fn seq_errors_test() {
        let mut sequence_errors = SequenceErrors::new();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 0,
                sample_barcode: 0,
                bb_barcode: 0,
                matched: 0,
                duplicates: 0,
            }
        );
        sequence_errors.correct_match();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 0,
                sample_barcode: 0,
                bb_barcode: 0,
                matched: 1,
                duplicates: 0,
            }
        );
        sequence_errors.constant_region_error();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 1,
                sample_barcode: 0,
                bb_barcode: 0,
                matched: 1,
                duplicates: 0,
            }
        );
        sequence_errors.sample_barcode_error();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 1,
                sample_barcode: 1,
                bb_barcode: 0,
                matched: 1,
                duplicates: 0,
            }
        );
        sequence_errors.bb_barcode_error();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 1,
                sample_barcode: 1,
                bb_barcode: 1,
                matched: 1,
                duplicates: 0,
            }
        );
        sequence_errors.duplicated();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 1,
                sample_barcode: 1,
                bb_barcode: 1,
                matched: 1,
                duplicates: 1,
            }
        );
    }

    #[test]
    fn bb_barcode_file_conversion_test() {
        let bb_barcodes =
            bb_barcode_file_conversion(&"building_block_barcode.example.csv".to_string(), 3)
                .unwrap();
        let mut bb_comparison = HashMap::new();
        for barcode_num in [1usize, 2, 3] {
            if barcode_num == 1 {
                let start_hash: HashMap<String, String> = [
                    ("CAGAGAC".to_string(), "BB_name_1".to_string()),
                    ("TGATTGC".to_string(), "BB_name_2".to_string()),
                ]
                .iter()
                .cloned()
                .collect();
                bb_comparison.insert(barcode_num, start_hash);
            }
            if barcode_num == 2 {
                let start_hash: HashMap<String, String> = [
                    ("ATGAAAT".to_string(), "BB_name_3".to_string()),
                    ("GCGCCAT".to_string(), "BB_name_4".to_string()),
                ]
                .iter()
                .cloned()
                .collect();
                bb_comparison.insert(barcode_num, start_hash);
            }
            if barcode_num == 3 {
                let start_hash: HashMap<String, String> = [
                    ("GATAGCT".to_string(), "BB_name_5".to_string()),
                    ("TTAGCTA".to_string(), "BB_name_6".to_string()),
                ]
                .iter()
                .cloned()
                .collect();
                bb_comparison.insert(barcode_num, start_hash);
            }
        }
        assert_eq!(bb_barcodes, bb_comparison);
    }

    #[test]
    fn sample_barcode_file_conversion_test() {
        let sample_barcodes =
            sample_barcode_file_conversion(&"sample_barcode.example.csv".to_string()).unwrap();
        let sample_barcodes_comparison: HashMap<String, String> = [
            ("AGCATAC".to_string(), "Sample_name_1".to_string()),
            ("AACTTAC".to_string(), "Sample_name_2".to_string()),
        ]
        .iter()
        .cloned()
        .collect();
        assert_eq!(sample_barcodes, sample_barcodes_comparison);
    }
}
