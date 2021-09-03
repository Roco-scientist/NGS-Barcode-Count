use std::{error::Error, fs};

pub fn regex_search(format: String) -> Result<String, Box<dyn Error>> {
    let format_data = fs::read_to_string(format)?.replace("\n", "");
    println!("Format: {}", &format_data);
    Ok(String::new())
}
