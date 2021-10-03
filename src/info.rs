use itertools::Itertools;
use num_format::{Locale, ToFormattedString};
use regex::Regex;
use std::{
    collections::{HashMap, HashSet},
    error::Error,
    fmt, fs,
    sync::{
        atomic::{AtomicBool, AtomicU32, AtomicUsize, Ordering},
        {Arc, Mutex},
    },
};

// Struct to keep track of sequencing errors and correct matches.  This is displayed at the end of the algorithm for QC measures
#[derive(Debug, Clone)]
pub struct SequenceErrors {
    constant_region: Arc<AtomicU32>, // errors within the constant region
    sample_barcode: Arc<AtomicU32>,  // errors within the sample barcode
    barcode: Arc<AtomicU32>,         // erors within the counted barcode
    matched: Arc<AtomicU32>,         // total matched
    duplicates: Arc<AtomicU32>,      // total random barcode duplicates
    low_quality: Arc<AtomicU32>,     // total random barcode duplicates
}

impl Default for SequenceErrors {
    fn default() -> Self {
        Self::new()
    }
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
    pub fn new() -> Self {
        SequenceErrors {
            constant_region: Arc::new(AtomicU32::new(0)),
            sample_barcode: Arc::new(AtomicU32::new(0)),
            barcode: Arc::new(AtomicU32::new(0)),
            matched: Arc::new(AtomicU32::new(0)),
            duplicates: Arc::new(AtomicU32::new(0)),
            low_quality: Arc::new(AtomicU32::new(0)),
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
        self.constant_region.fetch_add(1, Ordering::Relaxed);
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
        self.sample_barcode.fetch_add(1, Ordering::Relaxed);
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
        self.barcode.fetch_add(1, Ordering::Relaxed);
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
        self.matched.fetch_add(1, Ordering::Relaxed);
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
        self.duplicates.fetch_add(1, Ordering::Relaxed);
    }

    /// Add one to low_quality
    ///
    /// # Example
    /// ```
    /// use barcode::barcode_info::SequenceErrors;
    ///
    /// let mut sequence_errors = SequenceErrors::new();
    /// sequence_errors.low_quality_barcode();
    /// ```
    pub fn low_quality_barcode(&mut self) {
        self.low_quality.fetch_add(1, Ordering::Relaxed);
    }

    pub fn arc_clone(&self) -> SequenceErrors {
        SequenceErrors {
            constant_region: Arc::clone(&self.constant_region),
            sample_barcode: Arc::clone(&self.sample_barcode),
            barcode: Arc::clone(&self.barcode),
            matched: Arc::clone(&self.matched),
            duplicates: Arc::clone(&self.duplicates),
            low_quality: Arc::clone(&self.low_quality),
        }
    }
}

impl fmt::Display for SequenceErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "\
            Correctly matched sequences: {}\n\
            Constant region mismatches:  {}\n\
            Sample barcode mismatches:   {}\n\
            Barcode mismatches:          {}\n\
            Duplicates:                  {}\n\
            Low quality barcodes:        {}",
            self.matched
                .load(Ordering::Relaxed)
                .to_formatted_string(&Locale::en),
            self.constant_region
                .load(Ordering::Relaxed)
                .to_formatted_string(&Locale::en),
            self.sample_barcode
                .load(Ordering::Relaxed)
                .to_formatted_string(&Locale::en),
            self.barcode
                .load(Ordering::Relaxed)
                .to_formatted_string(&Locale::en),
            self.duplicates
                .load(Ordering::Relaxed)
                .to_formatted_string(&Locale::en),
            self.low_quality
                .load(Ordering::Relaxed)
                .to_formatted_string(&Locale::en)
        )
    }
}

// Struct to keep the format information for the sequencing, ie barcodes, regex search etc.
#[derive(Debug)]
pub struct SequenceFormat {
    pub format_string: String,
    pub regions_string: String,
    pub length: usize,
    pub format_regex: Regex,
    pub regex_string: String,
    pub barcode_num: usize,
    format_data: String,
    pub random_barcode: bool,
    starts: Arc<Mutex<Vec<usize>>>,
    start_found: Arc<AtomicBool>,
    start: Arc<AtomicUsize>,
}

impl SequenceFormat {
    /// Creates a new SequenceFormat struct which holds the sequencing format information, such as, where the barcodes are located within the sequence
    pub fn new(format: String) -> Result<Self, Box<dyn Error>> {
        // Read sequenc format file to string
        let format_data = fs::read_to_string(format)?
            .lines() // split into lines
            .filter(|line| !line.starts_with('#')) // remove any line that starts with '#'
            .collect::<String>(); // collect into a String

        let regex_string = build_regex_captures(&format_data)?; // Build the regex string from the input format file information
        let random_barcode = regex_string.contains("random");
        let format_regex = Regex::new(&regex_string)?; // Convert the regex string to a Regex
        let format_string = build_format_string(&format_data)?; // Create the format string replacing 'N's where there is a barcode
        let regions_string = build_regions_string(&format_data)?; // Create the string which indicates where the barcodes are located
        let barcode_num = regex_string.matches("barcode").count(); // Count the number of barcodes.  This is used later for retrieving barcodes etc.
        let length = format_string.chars().count();

        // Create and return the SequenceFormat struct
        Ok(SequenceFormat {
            format_string,
            regions_string,
            length,
            format_regex,
            regex_string,
            barcode_num,
            format_data,
            random_barcode,
            starts: Arc::new(Mutex::new(Vec::new())),
            start_found: Arc::new(AtomicBool::new(false)),
            start: Arc::new(AtomicUsize::new(0)),
        })
    }

    /// Returns a Vec of the size of all counted barcodes within the seqeunce format
    pub fn barcode_lengths(&self) -> Result<Vec<u8>, Box<dyn Error>> {
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
                .parse::<u8>()?;
            // And add to the vector
            barcode_lengths.push(digits)
        }
        Ok(barcode_lengths)
    }

    /// Returns the sample barcode length found in the format file string
    pub fn sample_length_option(&self) -> Result<Option<u8>, Box<dyn Error>> {
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
                .parse::<u8>()?;
            Ok(Some(digits))
        } else {
            Ok(None)
        }
    }

    /// Returns the amount of nucleotides within the constant regions from the format file
    pub fn constant_region_length(&self) -> u8 {
        // Get the full length of the format_string and subtract the amount of 'N's found to get the constant nucleotide count
        let diff = self.format_string.len() - self.format_string.matches('N').count();
        diff as u8
    }

    pub fn add_start(&self, start: usize) {
        self.starts.lock().unwrap().push(start)
    }

    pub fn start_found(&self) -> bool {
        self.start_found.load(Ordering::Relaxed)
    }

    pub fn start(&self) -> usize {
        self.start.load(Ordering::Relaxed)
    }

    pub fn clone_arcs(&self) -> Self {
        SequenceFormat {
            format_string: self.format_string.clone(),
            regions_string: self.regions_string.clone(),
            length: self.length.clone(),
            format_regex: self.format_regex.clone(),
            regex_string: self.regex_string.clone(),
            barcode_num: self.barcode_num.clone(),
            format_data: self.format_data.clone(),
            random_barcode: self.random_barcode.clone(),
            starts: Arc::clone(&self.starts),
            start_found: Arc::clone(&self.start_found),
            start: Arc::clone(&self.start),
        }
    }
}

impl fmt::Display for SequenceFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut key = String::new();
        let mut new_char = HashSet::new();
        for key_char in self.regions_string.chars() {
            if new_char.insert(key_char) {
                let key_info = match key_char {
                    'S' => "\nS: Sample barcode",
                    'B' => "\nB: Counted barcode",
                    'C' => "\nC: Constant region",
                    'R' => "\nR: Random barcode",
                    _ => "",
                };
                key.push_str(key_info);
            }
        }
        write!(
            f,
            "-FORMAT-\n{}\n{}{}",
            self.format_string, self.regions_string, key
        )
    }
}

/// Builds the format string from the file format
///
/// # Example
///
/// ```
/// use barcode::barcode_info::build_format_string;
/// let format_data = "[8]AGCTAGATC{6}TGGA{6}TGGA{6}TGATTGCGC(6)NNNNAT";
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
            let digit_value = digit.as_str().parse::<u8>()?;
            for _ in 0..digit_value {
                final_format.push('N')
            }
        } else {
            final_format.push_str(group_str)
        }
    }
    Ok(final_format)
}

/// Builds the regions string from the file format
///
/// # Example
///
/// ```
/// use barcode::barcode_info::build_regions_string;
/// let format_data = "[8]AGCTAGATC{6}TGGA{6}TGGA{6}TGATTGCGC(6)NNNNAT";
///
/// assert_eq!(build_regions_string(format_data).unwrap(),  "SSSSSSSSCCCCCCCCCBBBBBBCCCCBBBBBBCCCCBBBBBBCCCCCCCCCRRRRRRCCCCCC".to_string())
/// ```
pub fn build_regions_string(format_data: &str) -> Result<String, Box<dyn Error>> {
    let digit_search = Regex::new(r"\d+")?;
    let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|[ATGCN]+")?;
    let mut final_format = String::new();
    for group in barcode_search.find_iter(format_data) {
        let group_str = group.as_str();
        let digits_option = digit_search.find(group_str);
        if let Some(digit) = digits_option {
            let digit_value = digit.as_str().parse::<u8>()?;
            for _ in 0..digit_value {
                if group_str.contains('[') {
                    final_format.push('S')
                } else if group_str.contains('{') {
                    final_format.push('B')
                } else if group_str.contains('(') {
                    final_format.push('R')
                }
            }
        } else {
            for _ in 0..group_str.chars().count() {
                final_format.push('C')
            }
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
                .parse::<u8>()
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

pub struct BarcodeConversions {
    pub samples_barcode_hash: HashMap<String, String>,
    pub sample_seqs: HashSet<String>,
    pub counted_barcodes_hash: Vec<HashMap<String, String>>,
    pub counted_barcode_seqs: Vec<HashSet<String>>,
}

impl Default for BarcodeConversions {
    fn default() -> Self {
        Self::new()
    }
}

impl BarcodeConversions {
    pub fn new() -> Self {
        BarcodeConversions {
            samples_barcode_hash: HashMap::new(),
            sample_seqs: HashSet::new(),
            counted_barcodes_hash: Vec::new(),
            counted_barcode_seqs: Vec::new(),
        }
    }

    /// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
    /// and the second needs to be the ID
    pub fn sample_barcode_file_conversion(
        &mut self,
        barcode_path: &str,
    ) -> Result<(), Box<dyn Error>> {
        // read in the sample barcode file
        for (barcode, sample_id) in fs::read_to_string(barcode_path)?
            .lines() // split the lines
            .skip(1) // skip the first line which should be the header
            .map(|line| {
                line.split(',')
                    .take(2) // take only the first two values, or columns
                    .map(|value| value.to_string())
                    .collect_tuple()
                    .unwrap_or(("".to_string(), "".to_string()))
            })
        {
            self.samples_barcode_hash.insert(barcode, sample_id);
        }
        Ok(())
    }

    /// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
    /// the second needs to be the ID, and the third needs to be the barcode index location
    ///
    /// # Panics
    ///
    /// This panics if the third column of the barcode conversion file does not contain integers.  Also
    /// panics if not all integers for barcode numbers is within this columns
    pub fn barcode_file_conversion(
        &mut self,
        barcode_path: &str,
        barcode_num: usize,
    ) -> Result<(), Box<dyn Error>> {
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
        for _ in 0..barcode_num {
            self.counted_barcodes_hash.push(HashMap::new());
        }
        let mut barcode_num_contained = HashSet::new();
        for (barcode, id, barcode_num) in barcode_vecs {
            let barcode_num_usize = barcode_num.parse::<usize>().unwrap_or_else(|err| {
            panic!("Third column of barcode file contains something other than an integer: {}\nError: {}", barcode_num, err)
        }) - 1;
            barcode_num_contained.insert(barcode_num_usize);
            self.counted_barcodes_hash[barcode_num_usize].insert(barcode, id);
        }
        let mut missing_barcode_num = Vec::new();
        for x in 0..barcode_num {
            if !barcode_num_contained.contains(&x) {
                missing_barcode_num.push(x)
            }
        }
        if !missing_barcode_num.is_empty() {
            panic!(
                "Barcode conversion file missing barcode numers {:?} in the third column",
                missing_barcode_num
            )
        }
        Ok(())
    }
    /// Creates a hashmap of all sample barcode sequences in order to compare for sequencing errors
    pub fn get_sample_seqs(&mut self) {
        if !self.samples_barcode_hash.is_empty() {
            for sample_barcode in self.samples_barcode_hash.keys() {
                self.sample_seqs.insert(sample_barcode.to_string());
            }
        }
    }

    /// Creates a hashmap of all counted barcode sequences in order to compare for sequencing errors
    pub fn get_barcode_seqs(&mut self) {
        if !self.counted_barcodes_hash.is_empty() {
            self.counted_barcode_seqs = self
                .counted_barcodes_hash
                .iter()
                .map(|hash| {
                    hash.keys()
                        .map(|key| key.to_string())
                        .collect::<HashSet<String>>()
                }) // creates a hashset for each sequential barcode, then collects into a vector with the index being each sequential counted barcode
                .collect::<Vec<HashSet<String>>>();
        }
    }
}

// Struct of how many sequencing errrors are allowed
#[derive(Debug, Clone, PartialEq)]
pub struct MaxSeqErrors {
    // errors within the constant region
    constant_region: u8,
    constant_region_size: u8,
    // errors within the sample barcode
    sample_barcode: u8,
    sample_size: u8,
    // erors within the counted barcode
    barcode: Vec<u8>,
    barcode_sizes: Vec<u8>,
    min_quality: f32,
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
        sample_errors_option: Option<u8>,
        sample_barcode_size_option: Option<u8>,
        barcode_errors_option: Option<u8>,
        barcode_sizes: Vec<u8>,
        constant_errors_option: Option<u8>,
        constant_region_size: u8,
        min_quality: f32,
    ) -> Result<Self, Box<dyn Error>> {
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

        let mut max_barcode_errors = Vec::new();
        // If max error was set by input arguments, use that value, otherwise calculate 20% of barcode size for max error
        for barcode_size in &barcode_sizes {
            if let Some(barcode_errors) = barcode_errors_option {
                max_barcode_errors.push(barcode_errors);
            } else {
                max_barcode_errors.push(barcode_size / 5);
            }
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
            min_quality,
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
    pub fn max_constant_errors(&self) -> u8 {
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
    pub fn max_sample_errors(&self) -> u8 {
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
    /// assert_eq!(max_sequence_errors.max_barcode_errors(), vec![1,1,1]);
    /// let barcode_sizes = vec![8,8,8];
    /// let barcode_errors_option = Some(2);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size).unwrap();
    /// assert_eq!(max_sequence_errors.max_barcode_errors(), vec![2,2,2]);
    /// ```
    pub fn max_barcode_errors(&self) -> &[u8] {
        &self.barcode
    }
}

impl fmt::Display for MaxSeqErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let barcode_size_info;
        let barcode_error_info;
        if self.barcode_sizes.len() > 1 {
            barcode_size_info = format!("Barcode sizes: {:?}", self.barcode_sizes);
            barcode_error_info = format!(
                "Maximum mismatches allowed per barcode sequence: {:?}",
                self.barcode
            );
        } else {
            barcode_size_info = format!("Barcode size: {}", self.barcode_sizes.first().unwrap());
            barcode_error_info = format!(
                "Maximum mismatches allowed per barcode sequence: {}",
                self.barcode.first().unwrap()
            );
        }
        write!(
            f,
            "\
            -BARCODE INFO-\n\
            Constant region size: {}\n\
            Maximum mismatches allowed per sequence: {}\n\
            --------------------------------------------------------------\n\
            Sample barcode size: {}\n\
            Maximum mismatches allowed per sequence: {}\n\
            --------------------------------------------------------------\n\
            {}\n\
            {}\n\
            --------------------------------------------------------------\n\
            Minimum allowed average read quality score per barcode: {}\n\
            ",
            self.constant_region_size,
            self.constant_region,
            self.sample_size,
            self.sample_barcode,
            barcode_size_info,
            barcode_error_info,
            self.min_quality
        )
    }
}

// An enum for whether or not the sequencing format contains a random barcode.  This changes how the barcodes are counted and which are kept
#[derive(Debug)]
pub enum FormatType {
    RandomBarcode,
    NoRandomBarcode,
}

// A struct which holds the count results, whether that is for a scheme which contains a random barcode or not
#[derive(Debug)]
pub struct Results {
    pub random_hashmap: HashMap<String, HashMap<String, HashSet<String>>>, // The counts for for schemes that contain random barcodes
    pub count_hashmap: HashMap<String, HashMap<String, u32>>, // The counts for the schemes which don't contain random barcodes.  Right now it is either on or the other.  Too much memory may be needed otherwise
    pub format_type: FormatType, // Whether it is with a random barcode or not
    empty_count_hash: HashMap<String, u32>, // An empty hashmap that is used a few times and therefor stored within the struct
}

impl Results {
    /// Create a new Results struct
    pub fn new(samples_barcode_hash: &HashMap<String, String>, random_barcode: bool) -> Self {
        // Record and keep the format type of whether or not the random barcode is included
        let format_type;
        if random_barcode {
            format_type = FormatType::RandomBarcode
        } else {
            format_type = FormatType::NoRandomBarcode
        }

        let mut random_hashmap = HashMap::new(); // create the hashmap that is used to count with random barcode scheme
        let mut count_hashmap = HashMap::new(); // create the hashmap that is used ot count when a random barcode is not used.  One or the other stays empty depending on the scheme

        // create empty hashmaps to insert and have the sample name included.  This is so sample name doesn't need to be searched each time
        let empty_random_hash: HashMap<String, HashSet<String>> = HashMap::new();
        let empty_count_hash: HashMap<String, u32> = HashMap::new();

        // If sample name conversion was included, add all sample names to the hashmaps used to count
        if !samples_barcode_hash.is_empty() {
            for sample in samples_barcode_hash.keys() {
                let sample_barcode = sample.to_string();
                random_hashmap.insert(sample_barcode.clone(), empty_random_hash.clone());
                count_hashmap.insert(sample_barcode, empty_count_hash.clone());
            }
        } else {
            // If sample names are not included, insert unknown name into the hashmaps
            random_hashmap.insert("Unknown_sample_name".to_string(), empty_random_hash);
            count_hashmap.insert("Unknown_sample_name".to_string(), empty_count_hash.clone());
        }
        // return the Results struct
        Results {
            random_hashmap,
            count_hashmap,
            format_type,
            empty_count_hash,
        }
    }

    /// Adds the random barcode connected to the barcode_ID that is counted and the sample.  When writing these unique barcodes are counted to get a count
    pub fn add_random(
        &mut self,
        sample_barcode: &str,
        random_barcode: &str,
        barcode_string: String,
    ) -> bool {
        // Get the hashmap for the sample
        let barcodes_hashmap_option;
        if sample_barcode.is_empty() {
            barcodes_hashmap_option = self.random_hashmap.get_mut("Unknown_sample_name");
        } else {
            barcodes_hashmap_option = self.random_hashmap.get_mut(sample_barcode);
        }
        if let Some(barcodes_hashmap) = barcodes_hashmap_option {
            // If the barcodes_hashmap is not empty
            // but doesn't contain the barcode
            if !barcodes_hashmap.contains_key(&barcode_string) {
                // insert the hashmap<barcode_id, Set<random_barcodes>>
                let mut intermediate_set = HashSet::new();
                intermediate_set.insert(random_barcode.to_string());
                barcodes_hashmap.insert(barcode_string, intermediate_set);
            } else {
                // if the hashmap<sample_id, hashmap<barcode_id, Set<>> exists, check to see if the random barcode already was inserted
                let random_set = barcodes_hashmap.get_mut(&barcode_string).unwrap();
                return random_set.insert(random_barcode.to_string());
            }
        } else {
            // create the Set<RandomBarcode>
            let mut intermediate_set = HashSet::new();
            intermediate_set.insert(random_barcode.to_string());
            let mut intermediate_hash = HashMap::new();
            // create the HashMap<barcode_id, Set<RandomBarcodes>>
            intermediate_hash.insert(barcode_string.to_string(), intermediate_set);
            // insert this into the random_hashmap connected to the sample_ID
            self.random_hashmap
                .insert(sample_barcode.to_string(), intermediate_hash);
        }
        true
    }

    /// Adds to the count for the barcode_id connected to the sample. This is used when a random barcode is not included in the scheme
    pub fn add_count(&mut self, sample_barcode: &str, barcode_string: String) {
        // Insert 0 if the barcodes are not within the sample_name -> barcodes
        // Then add one regardless
        *self
            .count_hashmap
            .get_mut(sample_barcode)
            .unwrap_or(&mut self.empty_count_hash.clone())
            .entry(barcode_string)
            .or_insert(0) += 1;
    }
}

pub struct ResultsEnrichment {
    pub single_hashmap: HashMap<String, HashMap<String, u32>>, // enrichment of single barcodes hash used at output
    pub double_hashmap: HashMap<String, HashMap<String, u32>>, // enrichment of double barcodes hash used at output
    empty_count_hash: HashMap<String, u32>,
}

impl ResultsEnrichment {
    pub fn new() -> Self {
        let empty_count_hash: HashMap<String, u32> = HashMap::new();
        ResultsEnrichment {
            single_hashmap: HashMap::new(),
            double_hashmap: HashMap::new(),
            empty_count_hash,
        }
    }

    pub fn add_sample_barcodes(&mut self, samples_barcodes: &[String]) {
        for sample_barcode in samples_barcodes {
            self.single_hashmap
                .insert(sample_barcode.to_string(), self.empty_count_hash.clone());
            self.double_hashmap
                .insert(sample_barcode.to_string(), self.empty_count_hash.clone());
        }
    }

    /// Adds the count the the single barcode enrichment hashmap
    pub fn add_single(&mut self, sample_id: &str, barcode_string: &str, count: u32) {
        let barcode_num = barcode_string.split(",").count();
        for (index, single_barcode) in barcode_string.split(",").enumerate() {
            let mut single_barcode_string = String::new();
            for x in 0..barcode_num {
                if x == index {
                    single_barcode_string.push_str(single_barcode);
                }
                if x != (barcode_num - 1) {
                    single_barcode_string.push(',');
                }
            }
            // Insert 0 if the barcodes are not within the sample_name -> barcodes
            // Then add one regardless
            *self
                .single_hashmap
                .get_mut(sample_id)
                .unwrap_or(&mut self.empty_count_hash.clone())
                .entry(single_barcode_string)
                .or_insert(0) += count;
        }
    }

    /// Adds the count to the double barcode enrichment hashmap
    pub fn add_double(&mut self, sample_id: &str, barcode_string: &str, count: u32) {
        let barcode_num = barcode_string.split(",").count();
        let barcode_split = barcode_string.split(",").collect::<Vec<&str>>();
        for index in 0..(barcode_num - 1) {
            let mut double_barcode_string = String::new();
            for second_index in 0..barcode_num {
                if second_index == index {
                    double_barcode_string.push_str(barcode_split[index])
                } else if second_index == (index + 1) {
                    double_barcode_string.push_str(barcode_split[(index + 1)])
                }
                if second_index != (barcode_num - 1) {
                    double_barcode_string.push(',')
                }
            }
            // Insert 0 if the barcodes are not within the sample_name -> barcodes
            // Then add one regardless
            *self
                .double_hashmap
                .get_mut(sample_id)
                .unwrap_or(&mut self.empty_count_hash.clone())
                .entry(double_barcode_string)
                .or_insert(0) += count;
        }
    }
}

impl Default for ResultsEnrichment {
    fn default() -> Self {
        Self::new()
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
                barcode: vec![1, 1, 1],
                barcode_sizes: vec![8, 8, 8],
            }
        )
    }
}
