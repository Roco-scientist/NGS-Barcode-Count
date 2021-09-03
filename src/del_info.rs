use std::{error::Error, fs};

// [11]GGAGGTCTTGCAGACAGAGGA{8}TCGT{8}CCAT{8}GCAAGCGATT
pub fn regex_search(format: String) -> Result<String, Box<dyn Error>> {
    let mut bb_num = 0;
    let format_data = fs::read_to_string(format)?
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.contains("#"))
        .map(|(line_num, line)| {
            let new_line = reformat_line(line.to_string(), line_num);
            new_line
        })
        .collect::<Vec<String>>()
        .join("");
    println!("Format: {}", &format_data);
    Ok(String::new())
}

fn reformat_line(mut line: String, bb_num: usize) -> String {
    if line.contains("{") {
        let mut replacement = format!("(?P<bb_{}>.", bb_num);
        replacement.push('{');
        line = line.replace("{", &replacement).replace("}", "})");
    }
    if line.contains("[") {
        line = line.replace("[", "(?P<sample>.{").replace("]", "})");
    }
    line
}
