use itertools::Itertools;
use regex::Regex;
use std::{collections::HashMap, error::Error, fs};

// Struct to keep track of sequencing errors
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SequenceErrors {
    // errors within the constant region
    constant_region: u64,
    // errors within the sample barcode
    sample_barcode: u64,
    // erors within the building block barcode
    bb_barcode: u64,
}

impl SequenceErrors {
    /// Create a new sequence error struct
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
            "Constant Region Mismatches: {}\nSample Barcode Mismatches: {}\nBuilding Block Mismatches: {}",
            self.constant_region, self.sample_barcode, self.bb_barcode
        )
    }
}

/// Reads in the sequencing format file and outputs a regex string with captures
pub fn regex_search(format: String) -> Result<String, Box<dyn Error>> {
    // Read sequenc format file to string
    let format_data = fs::read_to_string(format)?
        .lines() // split into lines
        .filter(|line| !line.starts_with("#")) // remove any line that starts with '#'
        .collect::<String>(); // collect into a String

    let final_format = build_regex_captures(format_data)?;
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
/// assert_eq!(build_regex_captures(format_data).unwrap(),  "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6}).{4}AT".to_string())
/// ```
pub fn build_regex_captures(format_data: String) -> Result<String, Box<dyn Error>> {
    let digit_search = Regex::new(r"\d+")?;
    let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|N+|[ATGC]+")?;
    // the previous does not bumber each barcode but names each caputre with bb#
    // The '#' needs to bre replaced with teh sequential number
    let mut final_format = String::new();
    let mut bb_num = 0;
    // Fore each character, if the character is #, replace with the sequential building block number
    for group in barcode_search.find_iter(&format_data) {
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

/// Replaces the capture groups in the regex string with 'N's hte lenght of the barcode
/// returns the reformatted string
///
/// # Example
/// ```
/// use del::del_info;
/// let regex_string = "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6})";
/// let sequence_format = del_info::replace_group(regex_string);
/// assert_eq!(sequence_format.unwrap(), "NNNNNNNNAGCTAGATCNNNNNNTGGANNNNNNTGGANNNNNNTGATTGCGCNNNNNN")
/// ```
pub fn replace_group(regex_string: &str) -> Result<String, Box<dyn Error>> {
    // Create a search for all of the captures in order to remove
    let remove_search = Regex::new(r"(\(\?P<sample>.)|(\(\?P<random>.)|(\(\?P<bb\d+>.)|(\}\))|\{")?;
    // Create a capture for digits or sequences
    let captures_search = Regex::new(r"(\d+)|([ATGC]+)")?;

    // Clean the regex string by removing the captures only leaving the numbers
    let regex_string_cleaned = remove_search.replace_all(regex_string, "");

    // setup a new string to push into for the final return
    let mut new_format = String::new();

    // For each group of numbers or nucleotids, if it is a number, repeat 'N' for the amount of number.
    // If it is a group of nucleotides, just push those to the final string
    for capture in captures_search.find_iter(&regex_string_cleaned) {
        let found = capture.as_str().to_string();
        // If it is a number, repeat Ns, otherwise just add the nucleotides
        let digit_check = found.parse::<u32>();
        if let Ok(digit) = digit_check {
            for _ in 0..digit {
                new_format.push_str("N");
            }
        } else {
            new_format.push_str(&found);
        }
    }
    Ok(new_format)
}

/// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
/// and the second needs to be the ID
pub fn sample_barcode_file_conversion(
    barcode_path: String,
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
    barcode_path: String,
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
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct MaxSeqErrors {
    // errors within the constant region
    constant_region: usize,
    constant_size: usize,
    // errors within the sample barcode
    sample_barcode: usize,
    sample_size: usize,
    // erors within the building block barcode
    bb_barcode: usize,
    bb_size: usize,
}

impl MaxSeqErrors {
    /// Create a new sequence error struct
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let bb_errors_option = None;
    /// let constant_errors_option = None;
    /// let regex_string = "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6})".to_string();
    /// let constant_region_string = del::del_info::replace_group(&regex_string).unwrap();
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, bb_errors_option, constant_errors_option, &regex_string, &constant_region_string).unwrap();
    /// ```
    pub fn new(
        sample_errors_option: Option<usize>,
        bb_errors_option: Option<usize>,
        constant_errors_option: Option<usize>,
        regex_string: &String,
        constant_region_string: &String,
    ) -> Result<MaxSeqErrors, Box<dyn Error>> {
        // Find the length of the sample barcodes and set maximum error to either 20% of that length, or the input value from args
        // Create a regex to pull the number from the sample group
        let sample_search = Regex::new(r"sample>\.\{(?P<sample_size>\d+)")?;
        let sample_size_search_option = sample_search.captures(regex_string);
        let max_sample_errors;
        // start with a sample size of 0 in case there is no sample barcode.  If there is then mutate
        let mut sample_size = 0;
        // If sample barcode was included, calculate the maximum error, otherwise set error to 0
        if let Some(sample_size_search) = sample_size_search_option {
            // get sample size from the regex string
            sample_size = sample_size_search
                .name("sample_size")
                .unwrap()
                .as_str()
                .parse::<usize>()?;
            // if there was sample errors input from arguments, use that, otherwise calculate 20% for max errors
            if let Some(sample_errors) = sample_errors_option {
                max_sample_errors = sample_errors
            } else {
                max_sample_errors = sample_size / 5;
            }
        } else {
            max_sample_errors = 0;
        }

        // Find the length of the building block barcodes and set maximum error to either 20% of that length, or the input value from args
        let bb_search = Regex::new(r"bb\d+>\.\{(?P<bb_size>\d+)")?;
        let max_bb_errors;
        // Get the building block size from the first building block.  This assumes they are all the same size.  Can come back and make this better
        let bb_size = bb_search
            .captures(regex_string)
            .unwrap()
            .name("bb_size")
            .unwrap()
            .as_str()
            .parse::<usize>()?;
        // If max error was set by input arguments, use that value, otherwise calculate 20% of barcode size for max error
        if let Some(bb_errors) = bb_errors_option {
            max_bb_errors = bb_errors
        } else {
            max_bb_errors = bb_size / 5;
        }

        // Find the length of the constant region and set maximum error to either 20% of that length, or the input value from args
        let max_constant_errors;
        // Get the constant region size by subtracting the 'N's from the format string
        let constant_size =
            constant_region_string.len() - constant_region_string.matches("N").count();
        // If max error was set by input arguments, use that value, otherwise calculate 20% of barcode size for max error
        if let Some(constant_errors) = constant_errors_option {
            max_constant_errors = constant_errors
        } else {
            max_constant_errors = constant_size / 5;
            // errors allowed is the length of the constant region - the Ns / 5 or 20%
        }

        Ok(MaxSeqErrors {
            constant_region: max_constant_errors,
            constant_size,
            sample_barcode: max_sample_errors,
            sample_size,
            bb_barcode: max_bb_errors,
            bb_size,
        })
    }

    /// Returns the maximum allowed constant region errors
    ///
    /// # Example
    /// ```
    /// use del::del_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let bb_errors_option = None;
    /// let constant_errors_option = None;
    /// let regex_string = "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6})".to_string();
    /// let constant_region_string = del::del_info::replace_group(&regex_string).unwrap();
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, bb_errors_option, constant_errors_option, &regex_string, &constant_region_string).unwrap();
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 5);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, bb_errors_option, Some(3), &regex_string, &constant_region_string).unwrap();
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
    /// let bb_errors_option = None;
    /// let constant_errors_option = None;
    /// let regex_string = "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6})".to_string();
    /// let constant_region_string = del::del_info::replace_group(&regex_string).unwrap();
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, bb_errors_option, constant_errors_option, &regex_string, &constant_region_string).unwrap();
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 1);
    /// let mut max_sequence_errors = MaxSeqErrors::new(Some(2), bb_errors_option, constant_errors_option, &regex_string, &constant_region_string).unwrap();
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 2);
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
    /// let bb_errors_option = None;
    /// let constant_errors_option = None;
    /// let regex_string = "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6})".to_string();
    /// let constant_region_string = del::del_info::replace_group(&regex_string).unwrap();
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, bb_errors_option, constant_errors_option, &regex_string, &constant_region_string).unwrap();
    /// assert_eq!(max_sequence_errors.max_bb_errors(), 1);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, Some(2), constant_errors_option, &regex_string, &constant_region_string).unwrap();
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
    /// let bb_errors_option = None;
    /// let constant_errors_option = None;
    /// let regex_string = "(?P<sample>.{8})AGCTAGATC(?P<bb1>.{6})TGGA(?P<bb2>.{6})TGGA(?P<bb3>.{6})TGATTGCGC(?P<random>.{6})".to_string();
    /// let constant_region_string = del::del_info::replace_group(&regex_string).unwrap();
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, bb_errors_option, constant_errors_option, &regex_string, &constant_region_string).unwrap();
    /// max_sequence_errors.display();
    /// ```
    pub fn display(&mut self) {
        println!(
            "Constant region size: {}\n\
            Maximum constant region mismatches allowed per sequence: {}\n\
            Sample barcode size: {}\n\
            Maximum sample barcode mismatches allowed per sequence: {}\n\
            Building block size: {}\n\
            Maximum building block mismatches allowed per barcode: {}",
            self.constant_size,
            self.constant_region,
            self.sample_size,
            self.sample_barcode,
            self.bb_size,
            self.bb_barcode
        )
    }
}
