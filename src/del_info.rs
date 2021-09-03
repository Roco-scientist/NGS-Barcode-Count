use itertools::Itertools;
use std::{collections::HashMap, error::Error, fs};

pub fn regex_search(format: String) -> Result<String, Box<dyn Error>> {
    let format_data = fs::read_to_string(format)?
        .lines()
        .filter(|line| !line.starts_with("#"))
        .map(|line| {
            let new_line = reformat_line(line.to_string());
            new_line
        })
        .collect::<Vec<String>>()
        .join("");

    let mut final_format = String::new();
    let mut bb_num = 0;
    for letter in format_data.chars() {
        if letter == '#' {
            bb_num += 1;
            for bb_char in bb_num.to_string().chars() {
                final_format.push(bb_char);
            }
        } else {
            final_format.push(letter);
        }
    }
    println!("Format: {}", &final_format);
    Ok(final_format)
}

fn reformat_line(mut line: String) -> String {
    if line.contains("{") {
        line = line.replace("{", "(?P<bb#>.{").replace("}", "})");
    }
    if line.contains("[") {
        line = line.replace("[", "(?P<sample>.{").replace("]", "})");
    }
    line
}

pub fn sample_barcodes(barcode_path: String) -> Result<HashMap<String, String>, Box<dyn Error>> {
    let barcode_data: HashMap<String, String> = fs::read_to_string(barcode_path)?
        .lines()
        .skip(1)
        .map(|line| {
            line.split(",")
                .take(2)
                .map(|value| value.to_string())
                .collect_tuple()
                .unwrap_or(("".to_string(), "".to_string()))
        })
        .collect::<Vec<(String, String)>>()
        .iter()
        .cloned()
        .collect();
    Ok(barcode_data)
}
