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
    // erors within the counted barcode
    barcode: u64,
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
    /// use barcode::barcode_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// ```
    pub fn new() -> SequenceErrors {
        SequenceErrors {
            constant_region: 0,
            sample_barcode: 0,
            barcode: 0,
            matched: 0,
            duplicates: 0,
        }
    }

    /// Add one to constant region error
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::SequenceErrors;
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
    /// use barcode::barcode_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.sample_barcode_error();
    /// ```
    pub fn sample_barcode_error(&mut self) {
        self.sample_barcode += 1;
    }

    /// Add one to barcode error
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.barcode_error();
    /// ```
    pub fn barcode_error(&mut self) {
        self.barcode += 1;
    }

    /// Add one to correct match
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::SequenceErrors;
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
    /// use barcode::barcode_info::SequenceErrors;
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
    /// use barcode::barcode_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.barcode_error();
    /// sequence_errors.display();
    /// ```
    pub fn display(&mut self) {
        println!(
            "Correctly matched sequences: {}\nConstant region mismatches:  {}\nSample barcode mismatches:   {}\nBarcode mismatches:          {}\nDuplicates:                  {}",
            self.matched, self.constant_region, self.sample_barcode, self.barcode, self.duplicates
        );
        println!()
    }
}

#[derive(Debug, Clone)]
pub struct SequenceFormat {
    pub format_string: String,
    // Not implemented yet
    pub format_string_multiple: Option<Vec<String>>,
    pub format_regex: Regex,
    pub regex_string: String,
    pub barcode_num: usize,
    // Not implemented yet
    pub multiple: bool,
    format_data: String,
    pub random_barcode: bool,
}

impl SequenceFormat {
    /// Creates a new SequenceFormat struct which holds the sequencing format information, such as, where the barcodes are located within the sequence
    pub fn new(format: String) -> Result<SequenceFormat, Box<dyn Error>> {
        // Read sequenc format file to string
        let format_data = fs::read_to_string(format)?
            .lines() // split into lines
            .filter(|line| !line.starts_with('#')) // remove any line that starts with '#'
            .collect::<String>(); // collect into a String

        let regex_string = build_regex_captures(&format_data)?; // Build the regex string from the input format file information
        let random_barcode = regex_string.contains("random");
        let format_regex = Regex::new(&regex_string)?; // Convert the regex string to a Regex
        let format_string = build_format_string(&format_data)?; // Create the format string replacing 'N's where there is a barcode
        let barcode_num = regex_string.matches("barcode").count(); // Count the number of barcodes.  This is used later for retrieving barcodes etc.

        // Create and return the SequenceFormat struct
        Ok(SequenceFormat {
            format_string,
            format_string_multiple: None,
            format_regex,
            regex_string,
            barcode_num,
            multiple: false,
            format_data,
            random_barcode,
        })
    }

    /// Displays the sequence format information with 'N's replacing all barcodes
    pub fn display_format(&self) {
        println!("Format: {}", self.format_string);
        // println!();
    }

    /// Returns a Vec of the size of all counted barcodes within the seqeunce format
    pub fn barcode_lengths(&self) -> Result<Vec<usize>, Box<dyn Error>> {
        let barcode_search = Regex::new(r"(\{\d+\})")?; // Create a search that finds the '{#}'
        let digit_search = Regex::new(r"\d+")?; // Create a search that pulls out the number
        let mut barcode_lengths = Vec::new(); // Create a Vec that will contain the counted barcode lengths

        // For each counted barcode found in the format file string.  This allows multiple counted barcodes, as found in DEL
        for group in barcode_search.find_iter(&self.format_data) {
            let group_str = group.as_str();
            // Pull out the numeric value
            let digits = digit_search
                .find(group_str)
                .unwrap()
                .as_str()
                .parse::<usize>()?;
            // And add to the vector
            barcode_lengths.push(digits)
        }
        Ok(barcode_lengths)
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
            Ok(Some(digits))
        } else {
            Ok(None)
        }
    }

    /// Returns the amount of nucleotides within the constant regions from the format file
    pub fn constant_region_length(&self) -> usize {
        // Get the full length of the format_string and subtract the amount of 'N's found to get the constant nucleotide count
        self.format_string.len() - self.format_string.matches('N').count()
    }
}

/// Builds the catpure groups from the file format
///
/// # Example
///
/// ```
/// use barcode::barcode_info::build_format_string;
/// let format_data = "[8]AGCTAGATC{6}TGGA{6}TGGA{6}TGATTGCGC(6)NNNNAT".to_string();
///
/// assert_eq!(build_format_string(&format_data).unwrap(),  "NNNNNNNNAGCTAGATCNNNNNNTGGANNNNNNTGGANNNNNNTGATTGCGCNNNNNNNNNNAT".to_string())
/// ```
pub fn build_format_string(format_data: &str) -> Result<String, Box<dyn Error>> {
    let digit_search = Regex::new(r"\d+")?;
    let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|N+|[ATGC]+")?;
    let mut final_format = String::new();
    for group in barcode_search.find_iter(format_data) {
        let group_str = group.as_str();
        let digits_option = digit_search.find(group_str);
        if let Some(digit) = digits_option {
            let digit_value = digit.as_str().parse::<usize>()?;
            for _ in 0..digit_value {
                final_format.push('N')
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
/// use barcode::barcode_info::build_regex_captures;
/// let format_data = "[8]AGCTAGATC{6}TGGA{6}TGGA{6}TGATTGCGC(6)NNNNAT".to_string();
///
/// assert_eq!(build_regex_captures(&format_data).unwrap(),  "(?P<sample>.{8})AGCTAGATC(?P<barcode1>.{6})TGGA(?P<barcode2>.{6})TGGA(?P<barcode3>.{6})TGATTGCGC(?P<random>.{6}).{4}AT".to_string())
/// ```
pub fn build_regex_captures(format_data: &str) -> Result<String, Box<dyn Error>> {
    let digit_search = Regex::new(r"\d+")?;
    let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|N+|[ATGC]+")?;
    // the previous does not bumber each barcode but names each caputre with barcode#
    // The '#' needs to bre replaced with teh sequential number
    let mut final_format = String::new();
    let mut barcode_num = 0;
    // Fore each character, if the character is #, replace with the sequential barcode number
    for group in barcode_search.find_iter(format_data) {
        let group_str = group.as_str();
        let mut group_name_option = None;
        if group_str.contains('[') {
            group_name_option = Some("sample".to_string())
        } else if group_str.contains('{') {
            barcode_num += 1;
            group_name_option = Some(format!("barcode{}", barcode_num));
        } else if group_str.contains('(') {
            group_name_option = Some("random".to_string());
        }

        if let Some(group_name) = group_name_option {
            let digits = digit_search
                .captures(group_str)
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
        } else if group_str.contains('N') {
            let num_of_ns = group_str.matches('N').count();
            let mut n_group = ".{".to_string();
            n_group.push_str(&num_of_ns.to_string());
            n_group.push('}');
            final_format.push_str(&n_group);
        } else {
            final_format.push_str(&group_str.to_uppercase())
        }
    }
    Ok(final_format)
}

/// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
/// and the second needs to be the ID
pub fn sample_barcode_file_conversion(
    barcode_path: &str,
) -> Result<HashMap<String, String>, Box<dyn Error>> {
    // read in the sample barcode file
    let barcode_data: HashMap<String, String> = fs::read_to_string(barcode_path)?
        .lines() // split the lines
        .skip(1) // skip the first line which should be the header
        .map(|line| {
            line.split(',')
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
/// the second needs to be the ID, and the third needs to be the barcode index location
///
/// # Panics
///
/// This panics if the third column of the barcode conversion file does not contain integers.  Also
/// panics if not all integers for barcode numbers is within this columns
pub fn barcode_file_conversion(
    barcode_path: &str,
    barcode_num: usize,
) -> Result<HashMap<usize, HashMap<String, String>>, Box<dyn Error>> {
    // read in the sample barcode file
    let barcode_vecs = fs::read_to_string(barcode_path)?
        .lines() // split the lines
        .skip(1) // skip the first line which should be the header
        .map(|line| {
            line.split(',')
                .take(3) // take only the first three values, or columns
                .map(|value| value.to_string())
                .collect_tuple()
                .unwrap_or(("".to_string(), "".to_string(), "".to_string()))
        }) // comma split the line into a tuple with the first being the key and the last the value
        .collect::<Vec<(String, String, String)>>();
    let mut barcode_data = HashMap::new();
    let mut barcode_num_contained = Vec::new();
    for (barcode, id, barcode_num) in barcode_vecs {
        let barcode_num_usize = barcode_num.parse::<usize>().unwrap_or_else(|err| {
            panic!("Third column of barcode file contains something other than an integer: {}\nError: {}", barcode_num, err)
        });
        if !barcode_num_contained.contains(&barcode_num_usize) {
            barcode_num_contained.push(barcode_num_usize)
        };
        if let std::collections::hash_map::Entry::Vacant(e) = barcode_data.entry(barcode_num_usize)
        {
            let mut intermediate_hash = HashMap::new();
            intermediate_hash.insert(barcode, id);
            e.insert(intermediate_hash);
        } else {
            barcode_data
                .get_mut(&barcode_num_usize)
                .unwrap()
                .insert(barcode, id);
        }
    }
    let mut missing_barcode_num = Vec::new();
    for x in 0..barcode_num {
        let actual_num = x + 1;
        if !barcode_num_contained.contains(&actual_num) {
            missing_barcode_num.push(actual_num)
        }
    }
    if !missing_barcode_num.is_empty() {
        panic!(
            "Barcode conversion file missing barcode numers {:?} in the third column",
            missing_barcode_num
        )
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
    // erors within the counted barcode
    barcode: usize,
    barcode_sizes: Vec<usize>,
}

impl MaxSeqErrors {
    /// Create a new sequence error struct
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// ```
    pub fn new(
        sample_errors_option: Option<usize>,
        sample_barcode_size_option: Option<usize>,
        barcode_errors_option: Option<usize>,
        barcode_sizes: Vec<usize>,
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

        let max_barcode_errors;
        // If max error was set by input arguments, use that value, otherwise calculate 20% of barcode size for max error
        if let Some(barcode_errors) = barcode_errors_option {
            max_barcode_errors = barcode_errors
        } else {
            max_barcode_errors = barcode_sizes.iter().max().unwrap() / 5;
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
            barcode: max_barcode_errors,
            barcode_sizes,
        })
    }

    /// Returns the maximum allowed constant region errors
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 6);
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = Some(3);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 3);
    /// ```
    pub fn max_constant_errors(&self) -> usize {
        self.constant_region
    }

    /// Returns the maximum allowed sample barcode errors
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 2);
    /// let barcode_sizes = vec![8,8,8];
    /// let sample_errors_option = Some(3);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 3);
    /// ```
    pub fn max_sample_errors(&self) -> usize {
        self.sample_barcode
    }

    /// Returns the maximum allowed errors within each counted barcode
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_barcode_errors(), 1);
    /// let barcode_sizes = vec![8,8,8];
    /// let barcode_errors_option = Some(2);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_barcode_errors(), 2);
    /// ```
    pub fn max_barcode_errors(&self) -> usize {
        self.barcode
    }

    /// Print to stdout all maximum sequencing errors
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// max_sequence_errors.display();
    /// ```
    pub fn display(&mut self) {
        let barcode_size_info;
        if self.barcode_sizes.len() > 1 {
            barcode_size_info = format!("Barcode sizes: {:?}", self.barcode_sizes)
        } else {
            barcode_size_info = format!("Barcode size: {}", self.barcode_sizes.first().unwrap())
        }
        println!(
            "
            \n########## Barcode Info ###################################\n\
            Constant region size: {}\n\
            Maximum mismatches allowed per sequence: {}\n\
            -----------------------------------------------------------\n\
            Sample barcode size: {}\n\
            Maximum mismatches allowed per sequence: {}\n\
            -----------------------------------------------------------\n\
            {}\n\
            Maximum mismatches allowed per barcode sequence: {}\n\
            ###########################################################",
            self.constant_region_size,
            self.constant_region,
            self.sample_size,
            self.sample_barcode,
            barcode_size_info,
            self.barcode
        );
        println!();
    }
}

pub enum FormatType {
    RandomBarcode,
    NoRandomBarcode,
}

pub struct Results {
    pub random_hashmap: HashMap<String, HashMap<String, Vec<String>>>,
    pub count_hashmap: HashMap<String, HashMap<String, u32>>,
    pub format_type: FormatType,
    empty_count_hash: HashMap<String, u32>,
    empty_random_hash: HashMap<String, Vec<String>>,
}

impl Results {
    pub fn new(
        samples_hashmap_option: &Option<HashMap<String, String>>,
        random_barcode: bool,
    ) -> Results {
        let format_type;
        if random_barcode {
            format_type = FormatType::RandomBarcode
        } else {
            format_type = FormatType::NoRandomBarcode
        }
        let mut random_hashmap = HashMap::new();
        let mut count_hashmap = HashMap::new();
        let empty_random_hash: HashMap<String, Vec<String>> = HashMap::new();
        let empty_count_hash: HashMap<String, u32> = HashMap::new();
        if let Some(samples_hashmap) = samples_hashmap_option {
            for sample in samples_hashmap.keys() {
                let sample_name = sample.to_string();
                random_hashmap.insert(sample_name.clone(), empty_random_hash.clone());
                count_hashmap.insert(sample_name, empty_count_hash.clone());
            }
        } else {
            random_hashmap.insert("Unknown_sample_name".to_string(), empty_random_hash.clone());
            count_hashmap.insert("Unknown_sample_name".to_string(), empty_count_hash.clone());
        }
        Results {
            random_hashmap,
            count_hashmap,
            format_type,
            empty_count_hash,
            empty_random_hash,
        }
    }
    pub fn add_random(
        &mut self,
        sample_name: &String,
        random_barcode: String,
        barcode_string: &String,
    ) -> bool {
        // If it does not already have the sample name, insert the sample name -> building_block -> random_barcodes
        let barcodes_hashmap = self.random_hashmap.get_mut(sample_name).unwrap();
        // If the random hashmap does not have the building blocks yet, insert building_block -> random_barcodes
        if !barcodes_hashmap.contains_key(barcode_string) {
            let intermediate_vec = vec![random_barcode];
            barcodes_hashmap.insert(barcode_string.clone(), intermediate_vec);
        } else {
            let random_vec = barcodes_hashmap.get_mut(barcode_string).unwrap();
            // else check if the random barcode already used for the sample_name and building_blocks
            // if the random barcode is already in the vector, change add_value to false
            // otherqise add the random barcode to the random_barcodes vector
            if random_vec.contains(&random_barcode) {
                return true;
            } else {
                random_vec.push(random_barcode);
                return false;
            }
        }
        false
    }
    pub fn add_count(&mut self, sample_name: &String, barcode_string: &String) {
        // Insert 0 if the barcodes is not within the sample_name -> barcodes
        // Then add one regardless
        *self
            .count_hashmap
            .get_mut(sample_name)
            .unwrap_or(&mut self.empty_count_hash.clone())
            .entry(barcode_string.to_string())
            .or_insert(0) += 1;
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_sequence_errors_test() {
        let sample_errors_option = None;
        let sample_barcode_size_option = Some(10);
        let barcode_errors_option = None;
        let barcode_sizes = vec![8, 8, 8];
        let constant_errors_option = None;
        let constant_region_size = 30;
        let max_sequence_errors = MaxSeqErrors::new(
            sample_errors_option,
            sample_barcode_size_option,
            barcode_errors_option,
            barcode_sizes,
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
                barcode: 1,
                barcode_sizes: vec![8, 8, 8],
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
                barcode: 0,
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
                barcode: 0,
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
                barcode: 0,
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
                barcode: 0,
                matched: 1,
                duplicates: 0,
            }
        );
        sequence_errors.barcode_error();
        assert_eq!(
            sequence_errors,
            SequenceErrors {
                constant_region: 1,
                sample_barcode: 1,
                barcode: 1,
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
                barcode: 1,
                matched: 1,
                duplicates: 1,
            }
        );
    }

    #[test]
    fn barcode_file_conversion_test() {
        let barcodes = barcode_file_conversion(&"barcode.example.csv".to_string(), 3).unwrap();
        let mut barcode_comparison = HashMap::new();
        for barcode_num in [1usize, 2, 3] {
            if barcode_num == 1 {
                let start_hash: HashMap<String, String> = [
                    ("CAGAGAC".to_string(), "Barcode_name_1".to_string()),
                    ("TGATTGC".to_string(), "Barcode_name_2".to_string()),
                ]
                .iter()
                .cloned()
                .collect();
                barcode_comparison.insert(barcode_num, start_hash);
            }
            if barcode_num == 2 {
                let start_hash: HashMap<String, String> = [
                    ("ATGAAAT".to_string(), "Barcode_name_3".to_string()),
                    ("GCGCCAT".to_string(), "Barcode_name_4".to_string()),
                ]
                .iter()
                .cloned()
                .collect();
                barcode_comparison.insert(barcode_num, start_hash);
            }
            if barcode_num == 3 {
                let start_hash: HashMap<String, String> = [
                    ("GATAGCT".to_string(), "Barcode_name_5".to_string()),
                    ("TTAGCTA".to_string(), "Barcode_name_6".to_string()),
                ]
                .iter()
                .cloned()
                .collect();
                barcode_comparison.insert(barcode_num, start_hash);
            }
        }
        assert_eq!(barcodes, barcode_comparison);
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
