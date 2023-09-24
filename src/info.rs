use ahash::{AHashSet, HashMap, HashMapExt};
use anyhow::{anyhow, Context, Result};
use itertools::Itertools;
use num_format::{Locale, ToFormattedString};
use regex::Regex;
use std::{
    fmt, fs,
    sync::{
        atomic::{AtomicU32, Ordering},
        Arc,
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
    /// use barcode_count::info::SequenceErrors;
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
    /// use barcode_count::info::SequenceErrors;
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
    /// use barcode_count::info::SequenceErrors;
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
    /// use barcode_count::info::SequenceErrors;
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
    /// use barcode_count::info::SequenceErrors;
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
    /// use barcode_count::info::SequenceErrors;
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
    /// use barcode_count::info::SequenceErrors;
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
            Counted barcode mismatches:  {}\n\
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
#[derive(Debug, Clone)]
pub struct SequenceFormat {
    pub format_string: String,       // sequence with 'N's replacing barcodes
    pub regions_string: String,      // String with each region contain a code
    pub length: usize,               // Total length of format sequence
    pub constant_region_length: u16, // Length of only the consant nucleotides
    pub format_regex: Regex,         // The regex search used to find barcodes
    pub barcode_num: usize,          // Number of counted barcodes.  More for DEL
    pub barcode_lengths: Vec<u16>,   // The length of each counted barcode
    pub sample_length_option: Option<u16>, // Sample barcode length
    pub random_barcode: bool,        // Whether a random barcode is included
    pub sample_barcode: bool,        // Whether a sammple barcode is included
}

impl SequenceFormat {
    /// Creates a new empty SequenceFormat struct
    ///
    /// # Example
    /// ```
    /// use barcode_count::info::SequenceFormat;
    ///
    /// let sequence_format = SequenceFormat::new();
    /// ```
    pub fn new() -> Result<Self> {
        let empty_regex = Regex::new("")?;
        Ok(SequenceFormat {
            format_string: String::new(),
            regions_string: String::new(),
            length: 0,
            constant_region_length: 0,
            format_regex: empty_regex,
            barcode_num: 0,
            barcode_lengths: Vec::new(),
            sample_length_option: None,
            random_barcode: false,
            sample_barcode: false,
        })
    }
    /// Parses the format file into all fields of the SequenceFormat struct, including the regex
    /// search, barcode sizes, and sequence format strings.
    pub fn parse_format_file(format_path: &str) -> Result<Self> {
        let mut sequence_format = SequenceFormat::new()?;
        // Read sequence format file to string
        let format_data = fs::read_to_string(format_path)
            .context(format!("Failed to open {}", format_path))?
            .lines() // split into lines
            .filter(|line| !line.starts_with('#')) // remove any line that starts with '#'
            .collect::<String>(); // collect into a String

        // Starts the string that is used to create the regex search
        let mut regex_string = String::new();
        // Digit search to find the number within any format group
        let digit_search = Regex::new(r"\d+")?;
        // Search groups separated by '|' or statements in order to iterate through each group
        // within the format data from the format file and create the regex search string, along
        // with add the other needed information.  Uses the {#}, [#], (#), [ATGC], and 'N's as
        // groups
        let barcode_search = Regex::new(r"(?i)(\{\d+\})|(\[\d+\])|(\(\d+\))|N+|[ATGC]+")?;
        for group in barcode_search.find_iter(&format_data) {
            let group_str = group.as_str();
            // Holds the capture group name.  Is non-barcode regions
            let mut group_name_option = None;

            // If the group is a barcode group, add the capture group name, and set barcode
            // included fields to true
            if group_str.contains('[') {
                group_name_option = Some("sample".to_string());
                sequence_format.sample_barcode = true;
            } else if group_str.contains('{') {
                sequence_format.barcode_num += 1;
                group_name_option = Some(format!("barcode{}", sequence_format.barcode_num));
            } else if group_str.contains('(') {
                group_name_option = Some("random".to_string());
                sequence_format.random_barcode = true;
            }

            if let Some(group_name) = group_name_option {
                let digits = digit_search
                    .captures(group_str)
                    .unwrap()
                    .get(0)
                    .unwrap()
                    .as_str()
                    .parse::<u16>()
                    .unwrap();

                // Create the capture group with the group name for the barcode and add it to the
                // string created for the regex search
                let mut capture_group = format!("(?P<{}>.", group_name);
                capture_group.push('{');
                capture_group.push_str(&digits.to_string());
                capture_group.push_str("})");
                regex_string.push_str(&capture_group);

                // Add lengths of any of the barcodes to the sequence_format struct fields.  Also
                // set the code for the regions_string
                let mut push_char = '\0';
                if group_name == "sample" {
                    sequence_format.sample_length_option = Some(digits);
                    push_char = 'S'
                } else if group_name.contains("barcode") {
                    sequence_format.barcode_lengths.push(digits);
                    push_char = 'B'
                } else if group_name == "random" {
                    push_char = 'R'
                }
                // For the number of nucleotides of the barcode add 'N's to format string and the
                // push_char just set to regions_string
                for _ in 0..digits {
                    sequence_format.regions_string.push(push_char);
                    sequence_format.format_string.push('N')
                }
            } else if group_str.contains('N') {
                // Used to handle if 'N's are added to the format file.  These will be treated as
                // 'any' nucleotide for error handling and matching
                let num_of_ns = group_str.matches('N').count();
                let mut n_group = "[AGCT]{".to_string();
                n_group.push_str(&num_of_ns.to_string());
                n_group.push('}');
                regex_string.push_str(&n_group);
                sequence_format.format_string.push_str(group_str);
            } else {
                // Any A,G,C, or T is treated as constant region here
                regex_string.push_str(&group_str.to_uppercase());
                sequence_format.format_string.push_str(group_str);
                let constant_group_length = group_str.chars().count();
                for _ in 0..constant_group_length {
                    sequence_format.regions_string.push('C');
                }
                sequence_format.constant_region_length += constant_group_length as u16;
            }
        }
        sequence_format.length = sequence_format.format_string.chars().count();
        sequence_format.format_regex = Regex::new(&regex_string)?;
        Ok(sequence_format)
    }
}

impl fmt::Display for SequenceFormat {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut key = String::new();
        let mut new_char = AHashSet::new();
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

/// Contains all possible barcode sequences for error handling and barcode to ID conversion
pub struct BarcodeConversions {
    pub samples_barcode_hash: HashMap<String, String>,
    pub sample_seqs: AHashSet<String>,
    pub counted_barcodes_hash: Vec<HashMap<String, String>>,
    pub counted_barcode_seqs: Vec<AHashSet<String>>,
}

impl Default for BarcodeConversions {
    fn default() -> Self {
        Self::new()
    }
}

impl BarcodeConversions {
    /// Creates an empty BarcodeConversions struct
    pub fn new() -> Self {
        BarcodeConversions {
            samples_barcode_hash: HashMap::new(),
            sample_seqs: AHashSet::new(),
            counted_barcodes_hash: Vec::new(),
            counted_barcode_seqs: Vec::new(),
        }
    }

    /// Reads in comma separated barcode file (CSV).  The columns need to have headers.  The first column needs to be the nucleotide barcode
    /// and the second needs to be the ID
    pub fn sample_barcode_file_conversion(&mut self, barcode_path: &str) -> Result<()> {
        // read in the sample barcode file
        for (barcode, sample_id) in fs::read_to_string(barcode_path)
            .context(format!("Failed to open {}", barcode_path))?
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
    ) -> Result<()> {
        // read in the sample barcode file
        let barcode_vecs = fs::read_to_string(barcode_path)
            .context(format!("Failed to read {}", barcode_path))?
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
        let mut barcode_num_contained = AHashSet::new();
        for (barcode, id, barcode_num) in barcode_vecs {
            let barcode_num_usize = barcode_num.parse::<usize>().context(format!(
                "Third column of barcode file contains something other than an integer: {}",
                barcode_num
            ))? - 1;
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
            return Err(anyhow!(format!(
                "Barcode conversion file missing barcode numers {:?} in the third column",
                missing_barcode_num
            )));
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
                        .collect::<AHashSet<String>>()
                }) // creates a hashset for each sequential barcode, then collects into a vector with the index being each sequential counted barcode
                .collect::<Vec<AHashSet<String>>>();
        }
    }
}

/// Struct of how many sequencing errrors are allowed
#[derive(Debug, Clone, PartialEq)]
pub struct MaxSeqErrors {
    // errors within the constant region
    constant_region: u16,
    constant_region_size: u16,
    // errors within the sample barcode
    sample_barcode: u16,
    sample_size: u16,
    // erors within the counted barcode
    barcode: Vec<u16>,
    barcode_sizes: Vec<u16>,
    min_quality: f32,
}

impl MaxSeqErrors {
    /// Create a new sequence error struct
    ///
    /// # Example
    /// ```
    /// use barcode_count::info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let min_quality = 0.0;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// ```
    pub fn new(
        sample_errors_option: Option<u16>,
        sample_barcode_size_option: Option<u16>,
        barcode_errors_option: Option<u16>,
        barcode_sizes: Vec<u16>,
        constant_errors_option: Option<u16>,
        constant_region_size: u16,
        min_quality: f32,
    ) -> Self {
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

        MaxSeqErrors {
            constant_region: max_constant_errors,
            constant_region_size,
            sample_barcode: max_sample_errors,
            sample_size,
            barcode: max_barcode_errors,
            barcode_sizes,
            min_quality,
        }
    }

    /// Returns the maximum allowed constant region errors
    ///
    /// # Example
    /// ```
    /// use barcode_count::info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let min_quality = 0.0;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 6);
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = Some(3);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// assert_eq!(max_sequence_errors.max_constant_errors(), 3);
    /// ```
    pub fn max_constant_errors(&self) -> u16 {
        self.constant_region
    }

    /// Returns the maximum allowed sample barcode errors
    ///
    /// # Example
    /// ```
    /// use barcode_count::info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let min_quality = 0.0;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 2);
    /// let barcode_sizes = vec![8,8,8];
    /// let sample_errors_option = Some(3);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// assert_eq!(max_sequence_errors.max_sample_errors(), 3);
    /// ```
    pub fn max_sample_errors(&self) -> u16 {
        self.sample_barcode
    }

    /// Returns the maximum allowed errors within each counted barcode
    ///
    /// # Example
    /// ```
    /// use barcode_count::info::MaxSeqErrors;
    ///
    /// let sample_errors_option = None;
    /// let sample_barcode_size_option = Some(10);
    /// let barcode_errors_option = None;
    /// let barcode_sizes = vec![8,8,8];
    /// let constant_errors_option = None;
    /// let constant_region_size = 30;
    /// let min_quality = 0.0;
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// assert_eq!(max_sequence_errors.max_barcode_errors(), vec![1,1,1]);
    /// let barcode_sizes = vec![8,8,8];
    /// let barcode_errors_option = Some(2);
    /// let mut max_sequence_errors = MaxSeqErrors::new(sample_errors_option, sample_barcode_size_option, barcode_errors_option, barcode_sizes, constant_errors_option, constant_region_size, min_quality);
    /// assert_eq!(max_sequence_errors.max_barcode_errors(), vec![2,2,2]);
    /// ```
    pub fn max_barcode_errors(&self) -> &[u16] {
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

#[derive(Debug)]
pub enum ResultsHashmap {
    RandomBarcode(HashMap<String, HashMap<String, AHashSet<String>>>),
    NoRandomBarcode(HashMap<String, HashMap<String, usize>>),
}

// A struct which holds the count results, whether that is for a scheme which contains a random barcode or not
#[derive(Debug)]
pub struct Results {
    pub results_hashmap: ResultsHashmap, // holds the counted results
    empty_count_hash: HashMap<String, usize>, // An empty hashmap that is used a few times and therefor stored within the struct
    empty_random_hash: HashMap<String, AHashSet<String>>,
    sample_conversion_omited: bool,
}

impl Results {
    /// Create a new Results struct
    pub fn new(
        samples_barcode_hash: &HashMap<String, String>,
        random_barcode: bool,
        sample_barcode: bool,
    ) -> Self {
        let mut results_hashmap;
        // Create an empty hashmap into the enum depending on whether or not a random barcode is
        // included
        if random_barcode {
            results_hashmap = ResultsHashmap::RandomBarcode(HashMap::new());
        } else {
            results_hashmap = ResultsHashmap::NoRandomBarcode(HashMap::new());
        }

        // If sample name conversion was included, add all sample names to the hashmaps used to count
        let mut sample_conversion_omited = false;
        // create empty hashmaps to insert and have the sample name included.  This is so sample name doesn't need to be searched each time
        let empty_random_hash: HashMap<String, AHashSet<String>> = HashMap::new();
        let empty_count_hash: HashMap<String, usize> = HashMap::new();
        // If there is a sample barcode file included, add these as keys in the relevant count hashmap
        if !samples_barcode_hash.is_empty() {
            for sample in samples_barcode_hash.keys() {
                let sample_barcode = sample.to_string();
                match results_hashmap {
                    ResultsHashmap::RandomBarcode(ref mut random_hashmap) => {
                        random_hashmap.insert(sample_barcode.clone(), empty_random_hash.clone());
                    }
                    ResultsHashmap::NoRandomBarcode(ref mut count_hashmap) => {
                        count_hashmap.insert(sample_barcode, empty_count_hash.clone());
                    }
                }
            }
        } else if !sample_barcode {
            // If there is not a sample barcode within the format, add 'barcode' as key
            match results_hashmap {
                ResultsHashmap::RandomBarcode(ref mut random_hashmap) => {
                    random_hashmap.insert("barcode".to_string(), empty_random_hash.clone());
                }
                ResultsHashmap::NoRandomBarcode(ref mut count_hashmap) => {
                    count_hashmap.insert("barcode".to_string(), empty_count_hash.clone());
                }
            }
        } else {
            // If there is a sample barcode in the format but no sample barcode conversion file,
            // set the following to true to make sample DNA barcodes into keys later on
            sample_conversion_omited = true;
        }
        // return the Results struct
        Results {
            results_hashmap,
            empty_count_hash,
            empty_random_hash,
            sample_conversion_omited,
        }
    }

    /// Adds the count to results hashmap
    pub fn add_count(
        &mut self,
        sample_barcode: &str,
        random_barcode: Option<&String>,
        barcode_string: String,
    ) -> bool {
        // If conversion file does not exist, add the barcode as a key value
        if self.sample_conversion_omited {
            match self.results_hashmap {
                ResultsHashmap::NoRandomBarcode(ref mut count_hashmap) => {
                    if !count_hashmap.contains_key(sample_barcode) {
                        count_hashmap
                            .insert(sample_barcode.to_string(), self.empty_count_hash.clone());
                    };
                }
                ResultsHashmap::RandomBarcode(ref mut random_hashmap) => {
                    if !random_hashmap.contains_key(sample_barcode) {
                        random_hashmap
                            .insert(sample_barcode.to_string(), self.empty_random_hash.clone());
                    };
                }
            }
        };

        match self.results_hashmap {
            // If random barcode is not included, add the count to this hashmap
            ResultsHashmap::NoRandomBarcode(ref mut count_hashmap) => {
                *count_hashmap
                    .get_mut(sample_barcode)
                    .unwrap_or(&mut self.empty_count_hash.clone())
                    .entry(barcode_string)
                    .or_insert(0) += 1;
            }
            // If a random barcode is included, add the random barcode and later use the number of
            // random barcodes as the count
            ResultsHashmap::RandomBarcode(ref mut random_hashmap) => {
                // Get the hashmap for the sample
                let barcodes_hashmap_option = if sample_barcode.is_empty() {
                    random_hashmap.get_mut("barcode")
                } else {
                    random_hashmap.get_mut(sample_barcode)
                };
                if let Some(barcodes_hashmap) = barcodes_hashmap_option {
                    // If the barcodes_hashmap is not empty
                    // but doesn't contain the barcode
                    if let std::collections::hash_map::Entry::Vacant(e) = barcodes_hashmap.entry(barcode_string.clone()) {
                        // insert the hashmap<barcode_id, Set<random_barcodes>>
                        let mut intermediate_set = AHashSet::new();
                        intermediate_set
                            .insert(random_barcode.unwrap_or(&"".to_string()).to_string());
                        e.insert(intermediate_set);
                    } else {
                        // if the hashmap<sample_id, hashmap<barcode_id, Set<>> exists, check to see if the random barcode already was inserted
                        let random_set = barcodes_hashmap.get_mut(&barcode_string).unwrap();
                        return random_set
                            .insert(random_barcode.unwrap_or(&"".to_string()).to_string());
                    }
                } else {
                    // create the Set<RandomBarcode>
                    let mut intermediate_set = AHashSet::new();
                    intermediate_set.insert(random_barcode.unwrap_or(&"".to_string()).to_string());
                    let mut intermediate_hash = HashMap::new();
                    // create the HashMap<barcode_id, Set<RandomBarcodes>>
                    intermediate_hash.insert(barcode_string.to_string(), intermediate_set);
                    // insert this into the random_hashmap connected to the sample_ID
                    random_hashmap.insert(sample_barcode.to_string(), intermediate_hash);
                }
            }
        }

        // Return that a count was added.  An earlier return value is used for when a random
        // barcode is already within the results
        true
    }
}

/// A struct which holds hte enriched single and double counted barcodes.  Useful for DEL.  This struct is used during output.
pub struct ResultsEnrichment {
    pub single_hashmap: HashMap<String, HashMap<String, usize>>, // enrichment of single barcodes hash used at output
    pub double_hashmap: HashMap<String, HashMap<String, usize>>, // enrichment of double barcodes hash used at output
    empty_count_hash: HashMap<String, usize>,
}

impl ResultsEnrichment {
    pub fn new() -> Self {
        let empty_count_hash: HashMap<String, usize> = HashMap::new();
        ResultsEnrichment {
            single_hashmap: HashMap::new(),
            double_hashmap: HashMap::new(),
            empty_count_hash,
        }
    }

    /// Adds sample barcodes for keys within the hashmaps.  This is added later in order to first initiate the struct then add sample barcodes later
    pub fn add_sample_barcodes(&mut self, samples_barcodes: &[String]) {
        // For each sample barcode, create a sample barcode key to empty hashmap into single and double enrichment hashmaps
        for sample_barcode in samples_barcodes {
            self.single_hashmap
                .insert(sample_barcode.to_string(), self.empty_count_hash.clone());
            self.double_hashmap
                .insert(sample_barcode.to_string(), self.empty_count_hash.clone());
        }
    }

    /// Adds the count the the single barcode enrichment hashmap
    pub fn add_single(&mut self, sample_id: &str, barcode_string: &str, count: usize) {
        // get the number of barcodes to know homu much to iterate
        let barcode_num = barcode_string.split(',').count();
        // For each single barcode in the comma separate barcodes, create a new string with just one barcode and empty other columns
        for (index, single_barcode) in barcode_string.split(',').enumerate() {
            let mut single_barcode_string = String::new();
            // Recreate the new comma separated barcode with only one barcode
            for x in 0..barcode_num {
                // If the index is sthe same as x, add the single barcode.  This should put it in the right column
                if x == index {
                    single_barcode_string.push_str(single_barcode);
                }
                // Don't add a comma at the end
                if x != (barcode_num - 1) {
                    single_barcode_string.push(',');
                }
            }
            // Insert 0 if the barcodes are not within the single_hashmap -> barcodes
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
    pub fn add_double(&mut self, sample_id: &str, barcode_string: &str, count: usize) {
        // get the number of barcodes to know homu much to iterate
        let barcode_num = barcode_string.split(',').count();
        // split the barcodes into a vec from their comma separated form
        let barcode_split = barcode_string.split(',').collect::<Vec<&str>>();
        // iterate through the number of barcode_num - 1, and take this index for the first barcode
        for first_barcode_index in 0..(barcode_num - 1) {
            // Get the amount needed to add to the first index in order to get the second index.  This is iterated to account for the second being the next barcode or two away etc. Eg from 1,2,3 = 1,2,, and 1,,3
            for next_barcode_add in 1..(barcode_num - first_barcode_index) {
                // Initiate the new barcode string
                let mut double_barcode_string = String::new();
                // Iterate over each comma separated column and insert the barcode if needed
                for column_index in 0..barcode_num {
                    // If it is either the first or second barcode, add comma separated to the new string
                    if column_index == first_barcode_index {
                        double_barcode_string.push_str(barcode_split[first_barcode_index])
                    } else if column_index == (first_barcode_index + next_barcode_add) {
                        double_barcode_string
                            .push_str(barcode_split[first_barcode_index + next_barcode_add])
                    }
                    // If we are not on the last barcode, add a comma
                    if column_index != (barcode_num - 1) {
                        double_barcode_string.push(',')
                    }
                }
                // Insert 0 if the barcodes are not within the double_hashmap -> barcodes
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
}

impl Default for ResultsEnrichment {
    fn default() -> Self {
        Self::new()
    }
}
