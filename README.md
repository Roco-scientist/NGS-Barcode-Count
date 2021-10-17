# NGS-Barcode-Count
Counts barcodes from next generation sequencing data.  Works for DEL (DNA encoded libraries), high throughput CRISPR sequencing, barcode sequencing.  If the barcode file is included, the program will convert to barcode names and correct for errors. If a random barcode is included to collapse PCR duplicates, these duplicates will not be counted.<br>
<br>
<br>
For DEL analysis, a companion python package was created: <a href=https://github.com/Roco-scientist/DEL-Analysis>DEL-Analysis</a>
<br>
<br>
Multithreaded and low resource use.  Uses one thread to read and the rest to process the data, so at least a 2 threaded machine is essential.
This program does not store all data within RAM but instead sequentially processes the sequencing data in order to remain memory efficient.  
<br>
<br>
Error handling is defaulted at 20% maximum sequence error per constant region and barcode. This can be changed through CLI arguments.  The algorithm fixes any sequenced constant region or barcode with the best match possible.  If there are two or more best matches,
it is not counted.
<br>
<br>
Filtering by read quality score is also an option.  If used, each barcode has its read quality average calculated and if it is below the set threshold, the read is not counted.  The algorithm is defaulted to not filter unless the --min_quality argument is called.  See fastq documentation to understand read quality scores. The scores used are after ascii conversion and 33 subtraction.
<br><br>
Inspired by and some ideas adopted from <a href=https://github.com/sunghunbae/decode target="_blank" rel="noopener noreferrer">decode</a>

## Table of Contents
<ul>
<li><a href=#installation>Installation</a></li>
<li><a href=#files-needed>Files Needed</a></li>
<li><a href=#run>Run</a></li>
<li><a href=#uses>Uses</a></li>
<li><a href=#test-results>Test Results</a></li>
<li><a href=#notes>Notes</a></li>
</ul>

## Installation

### Rust installed locally: <a href=https://www.rust-lang.org/tools/install target="_blank" rel="noopener noreferrer">instructions here</a>

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### NGS-Barcode-Count downoaded and compiled

```
git clone https://github.com/Roco-scientist/NGS-Barcode-Count
cd NGS-Barcode-Count
cargo build --release
mv ./target/release/barcode ./
```

## Files Needed
Currently supports FASTQ, sequence format, sample barcode conversion, and building block barcode conversion.
<ul>
<li><a href=#fastq-file>FASTQ</a></li>
<li><a href=#sequence-format-file>Sequence format file</a></li>
<li><a href=#sample-barcode-file>Sample barcode file (optional)</a></li>
<li><a href=#counted-barcode-conversion-file>Counted barcode conversion file (optional)</a></li>
</ul>


### Fastq File
Accepts unzipped fastq files.<br>
Accepts gzipped fastq files, but if the program stops before the expected number of sequencing reads, unzip and rerun.

### Sequence Format File
The sequence format file should be a text file that is line separated by the type of format.  The following is supported where the '#' should be replaced by the number of nucleotides corresponding to the barcode:<br>
<table>
<tr>
<th>Sequence Type</th>
<th>File Code</th>
<th>Number Needed/Allowed</th>
</tr>
<td>Constant</td>
<td>ATGCN</td>
<td>1 or more</td>
<tr>
<td>Sample Barcode</td>
<td>[#]</td>
<td>0-1</td>
</tr>
<tr>
<td>Barcode for counting</td>
<td>{#}</td>
<td>1 or more</td>
</tr>
<tr>
<td>Random Barcode</td>
<td>(#)</td>
<td>0-1</td>
</tr>
</table>

An example can be found in [scheme.example.txt](scheme.example.txt).  Since the algorthm uses a regex search to find the scheme, the scheme can exist anywhere within the sequence read.

### Sample Barcode File
<b>Optional</b><br>
The sample_barcode_file is a comma separate file with the following format:<br>
<table>
<tr>
<th>Barcode</th>
<th>Sample_ID</th>
</tr>
<tr>
<td>AGCATAC</td>
<td>Sample_name_1</td>
</tr>
<tr>
<td>AACTTAC</td>
<td>Sample_name_2</td>
</tr>
</table>

An example can be found in [sample_barcode.example.csv](sample_barcode.example.csv).

### Counted Barcode Conversion File
<b>Optional</b><br>
The barcode_file is a comma separate file with the following format:<br>
<table>
<tr>
<th>Barcode</th>
<th>Barcode_ID</th>
<th>Barcode_Number</th>
</tr>
<tr>
<td>CAGAGAC</td>
<td>Barcode_name_1</td>
<td>1</td>
</tr>
<tr>
<td>TGATTGC</td>
<td>Barcode_name_2</td>
<td>1</td>
</tr>
<tr>
<td>ATGAAAT</td>
<td>Barcode_name_3</td>
<td>2</td>
</tr>
<tr>
<td>GCGCCAT</td>
<td>Barcode_name_4</td>
<td>2</td>
</tr>
<tr>
<td>GATAGCT</td>
<td>Barcode_name_5</td>
<td>3</td>
</tr>
<tr>
<td>TTAGCTA</td>
<td>Barcode_name_6</td>
<td>3</td>
</tr>
</table>

An example can be found in [barcode.example.csv](barcode.example.csv).<br><br>

Where the first column is the DNA barcode, the second column is the barcode ID which can be a smile string for DEL, CRISPR target ID, etc. but cannot contain commas. 
The last column is the barcode number as an integer.  The barcode numbers are in the same order as the sequence format file and starting
at 1. For example, if there are a total of 3 barcodes, which may be the case with DEL, you would only have 1, 2, or 3 within this column for each row, with each number
representing one of the three barcodes. For CRISPR or barcode seq, where there may only be one barcode to count, this column would be all 1s.

## Run
After compilation, the `barcode` binary can be moved anywhere.
<br>
<br>
Run DEL-Decode<br>

```
barcode --fastq <fastq_file> \
	--sample_barcodes <sample_barcode_file> \
	--sequence_format <sequence_format_file> \
	--barcodes_file <barcodes_file> \
	--output_dir <output_dir> \
	--prefix <file_prefix> \
	--threads <num_of_threads> \
	--merge_output \
	--min_quality <min_barcode_read_quality>\
	--single\
	--double\
```

<br>
<ul>
<li>
--barcodes_file is optional.  If it is not used, the output counts uses the DNA barcode to count with no error handling on these barcodes.
</li>
<li>
--sample_barcodes is optional.  If it is not used, all samples are marked as unknown.
</li>
<li>
--output_dir defaults to the current directory if not used.
</li>
<li>
--prefix defaults to the current date.  All files end with _sample_name_counts.csv
</li>
<li>
--threads defaults to the number of threads on the machine if not used.
</li>
<li>
--merge_output flag that merges the output csv file so that each sample has one column
</li>
<li>
--min_quality will filter out reads where any of the barcodes have an average quality score below the threshold set here.  Default is 0 and no filtering.
</li>
<li>
--single argument flag that will find the counts for each barcode if there are 2 or more counted barcodes included, and output the file. Useful for DEL
</li>
<li>
--double argument flag that will find the counts for each pair of barcodes if there are 3 or more counted barcodes included, and output the file. Useful for DEL
</li>
</ul>

### Output files
Each sample name will get a file in the default format of year-month-day_<sample_name>_counts.csv in the following format (for 3 building blocks):
<table>
<tr>
<th>Barcode_1</th>
<th>Barcode_2</th>
<th>Barcode_3</th>
<th>Count</th>
</tr>
<tr>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>#</td>
</tr>
<tr>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>#</td>
</tr>
</table>

Where Barcode_ID is used if there is a building block conversion file, otherwise the DNA code is used. `#` represents the count number<br><br>
If `--merge_output` is called, an additional file is created with the format (for 3 samples):

<table>
<tr>
<th>Barcode_1</th>
<th>Barcode_2</th>
<th>Barcode_3</th>
<th>Sample_1</th>
<th>Sample_2</th>
<th>Sample_3</th>
</tr>
<tr>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>#</td>
<td>#</td>
<td>#</td>
</tr>
<tr>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>Barcode_ID/DNA code</td>
<td>#</td>
<td>#</td>
<td>#</td>
</tr>
</table>

An additional barcode_stats.txt file is also written/appended to the output folder.  This keeps track of running information.<br><br>
If either `--single` or `--double` arguments are called, single or double barcode count files are ouptut.

## Uses

### DEL
Setup as shown with all example files used throughout this README.  Typically you will use 3 x '[]' for counting barcodes, which represents 3 building blocks, within the format file.

### CRISPR-seq
Same setup as with DEL, but typically with only one '[]' counted barcode in the format file.  As such, within the counted barcode conversion file, the third column will contain all '1's

### Barcode-seq
If the intention is to count the random barcodes and have the counts associated with these random barcodes, which is the case with bar-seq of cell pools for lineage evolution etc., 
then the random barcode, within this situation, is the counted barcode and represented with '[]' in the format file.  A counted barcode conversion file will not be included.  Without the counted barcode conversion file, 
the program will output the counted random barcode sequence and the associated count.  Afterwards, clustering or any other analysis can be applied.

## Tests results
On an 8 threaded i7-4790K CPU @ 4.00GHz with 16gb RAM, this algorithm was able to decode over 400 million sequencing reads in about a half hour.
Results below:
```
Total sequences:             418,770,347
Correctly matched sequences: 257,807,865
Constant region mismatches:  151,955,695
Sample barcode mismatches:   3,324,481
Barcode mismatches:          5,682,306
Duplicates:                  0
Low quality barcodes:        0

Compute time: 0 hours, 27 minutes, 55.717 seconds

Writing counts

Total time: 0 hours, 34 minutes, 7.706 seconds
```

## Notes
In order to remain memory efficient, there are limits on how large each number can get<br>
Any count: u32 with a max of 4,294,967,295 <br>
Barcode lengths, each error max, number of barcodes in a single sequence: u8 with a max of 255 <br><br>
If larger values are needed, edit the script and replace the u32 with u64 or u8 with u16 etc.
