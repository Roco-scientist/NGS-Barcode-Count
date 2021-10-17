pub mod info;
pub mod input;
pub mod output;
pub mod parse;

use chrono::Local;
use clap::{App, Arg};
// use rayon::prelude::*;
use std::error::Error;

/// A struct that contains and initiates all input arguments
pub struct Args {
    pub fastq: String,                          // fastq file path
    pub format: String,                         // format scheme file path
    pub sample_barcodes_option: Option<String>, // sample barcode file path.  Optional
    pub barcodes_option: Option<String>,        // building block barcode file path. Optional
    pub output_dir: String,                     // output directory.  Deafaults to './'
    pub threads: u8, // Number of threads to use.  Defaults to number of threads on the machine
    pub prefix: String, // Prefix string for the output files
    pub merge_output: bool, // Whether or not to create an additional output file that merges all samples
    pub barcodes_errors_option: Option<u8>, // Optional input of how many errors are allowed in each building block barcode.  Defaults to 20% of the length
    pub sample_errors_option: Option<u8>, // Optional input of how many errors are allowed in each sample barcode.  Defaults to 20% of the length
    pub constant_errors_option: Option<u8>, // Optional input of how many errors are allowed in each constant region barcode.  Defaults to 20% of the length
    pub min_average_quality_score: f32,
    pub double_barcode_enrichment: bool,
    pub single_barcode_enrichment: bool,
}

impl Args {
    pub fn new() -> Result<Self, Box<dyn Error>> {
        let total_cpus = num_cpus::get().to_string();
        let today = Local::today().format("%Y-%m-%d").to_string();
        // parse arguments
        let args = App::new("NGS-Barcode-Count")
        .version("0.8.2")
        .author("Rory Coffey <coffeyrt@gmail.com>")
        .about("Counts barcodes located in sequencing data")
        .arg(
            Arg::with_name("fastq")
                .short("f")
                .long("fastq")
                .takes_value(true)
                .required(true)
                .help("FASTQ file unzipped"),
        )
        .arg(
            Arg::with_name("sequence_format")
                .short("q")
                .long("sequence_format")
                .takes_value(true)
                .required(true)
                .help("Sequence format file"),
        )
        .arg(
            Arg::with_name("sample_barcodes")
                .short("s")
                .long("sample_barcodes")
                .takes_value(true)
                .help("Sample barcodes file"),
        )
        .arg(
            Arg::with_name("barcodes_file")
                .short("b")
                .long("barcodes_file")
                .takes_value(true)
                .help("Building block barcodes file"),
        )
        .arg(
            Arg::with_name("threads")
                .short("t")
                .long("threads")
                .takes_value(true)
                .default_value(&total_cpus)
                .help("Number of threads"),
        )
        .arg(
            Arg::with_name("output_dir")
                .short("o")
                .long("output_dir")
                .takes_value(true)
                .default_value("./")
                .help("Directory to output the counts to"),
        )
        .arg(
            Arg::with_name("prefix")
                .short("p")
                .long("prefix")
                .takes_value(true)
                .default_value(&today)
                .help("File prefix name.  THe output will end with '_<sample_name>_counts.csv'"),
        )
        .arg(
            Arg::with_name("merge_output")
                .short("m")
                .long("merge_output")
                .takes_value(false)
                .help("Merge sample output counts into a single file.  Not necessary when there is only one sample"),
        )
        .arg(
            Arg::with_name("double")
                .long("double")
                .takes_value(false)
                .help("Output double barcode enrichment"),
        )
        .arg(
            Arg::with_name("single")
                .long("single")
                .takes_value(false)
                .help("Output single barcode enrichment"),
        )
        .arg(
            Arg::with_name("barcodes_errors")
                .long("barcodes_errors")
                .takes_value(true)
                .help("Maximimum number of sequence errors allowed within each counted barcode. Defaults to 20% of the total."),
        )
        .arg(
            Arg::with_name("sample_errors")
                .long("sample_errors")
                .takes_value(true)
                .help("Maximimum number of sequence errors allowed within sample barcode. Defaults to 20% of the total."),
        )
        .arg(
            Arg::with_name("contant_errors")
                .long("constant_errors")
                .takes_value(true)
                .help("Maximimum number of sequence errors allowed within constant region. Defaults to 20% of the total."),
        )
        .arg(
            Arg::with_name("min_quality")
                .long("min_quality")
                .takes_value(true)
                .default_value("0")
                .help("Minimum average read quality score per barcode"),
        )
        .get_matches();

        let sample_barcodes_option;
        if let Some(sample) = args.value_of("sample_barcodes") {
            sample_barcodes_option = Some(sample.to_string())
        } else {
            sample_barcodes_option = None
        }

        let barcodes_option;
        if let Some(barcodes) = args.value_of("barcodes_file") {
            barcodes_option = Some(barcodes.to_string())
        } else {
            barcodes_option = None
        }

        let barcodes_errors_option;
        if let Some(barcodes) = args.value_of("barcodes_errors") {
            barcodes_errors_option = Some(barcodes.parse::<u8>()?)
        } else {
            barcodes_errors_option = None
        }

        let sample_errors_option;
        if let Some(sample) = args.value_of("sample_errors") {
            sample_errors_option = Some(sample.parse::<u8>()?)
        } else {
            sample_errors_option = None
        }

        let constant_errors_option;
        if let Some(constant) = args.value_of("constant_errors") {
            constant_errors_option = Some(constant.parse::<u8>()?)
        } else {
            constant_errors_option = None
        }

        let merge_output;
        if args.is_present("merge_output") {
            merge_output = true
        } else {
            merge_output = false
        }
        let double_barcode_enrichment;
        if args.is_present("double") {
            double_barcode_enrichment = true
        } else {
            double_barcode_enrichment = false
        }
        let single_barcode_enrichment;
        if args.is_present("single") {
            single_barcode_enrichment = true
        } else {
            single_barcode_enrichment = false
        }
        let fastq = args.value_of("fastq").unwrap().to_string();
        let format = args.value_of("sequence_format").unwrap().to_string();
        let output_dir = args.value_of("output_dir").unwrap().to_string();
        let threads = args.value_of("threads").unwrap().parse::<u8>().unwrap();
        let prefix = args.value_of("prefix").unwrap().to_string();
        let min_average_quality_score = args
            .value_of("min_quality")
            .unwrap()
            .parse::<f32>()
            .unwrap();

        Ok(Args {
            fastq,
            format,
            sample_barcodes_option,
            barcodes_option,
            output_dir,
            threads,
            prefix,
            merge_output,
            barcodes_errors_option,
            sample_errors_option,
            constant_errors_option,
            min_average_quality_score,
            double_barcode_enrichment,
            single_barcode_enrichment,
        })
    }
}
