# DEL-Decode
DNA encoded library decoding.  Multithreaded and low resource use.  Uses one thread to read and the rest to process the data, so at least a 2 threaded machine is essential.
This program does not store all data within RAM but instead sequentially processes the sequencing data in order to remain memory efficient.  
With very large DEL libraries, there may be some memory creep.  This is because the more barcodes that exist within the sequencing reads, the more barcodes that need to remain in memory.
<br>
<br>
Error handling is set at 20% maximum sequence error.  The algorithm fixes any sequenced constant region or barcode with the best match possible.  If there are two or more best matches,
it is not counted.
<br>
<br>
Inspired by and some ideas adopted from <a href=https://github.com/sunghunbae/decode target="_blank" rel="noopener noreferrer">decode</a>

## Table of Contents
<ul>
<li><a href=#Requirements>Requirements</a></li>
<li><a href=#files-needed>Files Needed</a></li>
<li><a href=#run>Run</a></li>
<li><a href=#test-results>Test Results</a></li>
</ul>

## Requirements

### Rust installed locally: <a href=https://www.rust-lang.org/tools/install target="_blank" rel="noopener noreferrer">instructions here</a>

```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Del-Decode downoaded and compiled

```
git clone https://github.com/Roco-scientist/DEL-Decode.git
cd Del-Decode
cargo build --release
mv ./target/release/del ./
```

## Files Needed
Currently supports FASTQ, sequence format, sample barcode conversion, and building block barcode conversion.
<ul>
<li><a href=#fastq-file>FASTQ</a></li>
<li><a href=#sequence-format-file>Sequence format file</a></li>
<li><a href=#sample-barcode-file>Sample barcode file (optional)</a></li>
<li><a href=#building-block-barcode-file>Building block barcode file (optional)</a></li>
</ul>


### Fastq File
Only accepts unzipped FASTQ files because the program reads line by line.<br>
Importing gzipped files is currently in the works, but not yet supported.  A gzip inflate stream will be needed so that the whole gzipped file is not placed in RAM.

### Sequence Format File
The sequence format file should be a text file that is line separated by the type of format.  The following is supported where the '#' should be replaced by the number of nucleotides corresponding to the barcode:<br>
<table>
<tr>
<th>Sequence Type</th>
<th>File Code</th>
</tr>
<td>Constant</td>
<td>ATGCN</td>
<tr>
<td>Sample Barcode</td>
<td>[#]</td>
</tr>
<tr>
<td>Building Block Barcode</td>
<td>{#}</td>
</tr>
<tr>
<td>Random Barcode</td>
<td>(#)</td>
</tr>
</table>

An example can be found in [example.scheme.txt](example.scheme.txt)

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

### Building Block Barcode File
<b>Optional</b><br>
The building_block_barcode_file is a comma separate file with the following format:<br>
<table>
<tr>
<th>Barcode</th>
<th>BB_ID</th>
<th>BB_Number</th>
</tr>
<tr>
<td>CAGAGAC</td>
<td>BB_name_1</td>
<td>1</td>
</tr>
<tr>
<td>AACTTAC</td>
<td>BB_name_2</td>
<td>3</td>
</tr>
</table>
Where the first column is the DNA barcode, the second column is the building block ID which can be a smile string (without commas),
and the last column is the building block number as an integer.  The building block numbers are in the same order as the sequence format file and starting
at 1. For example, if there are a total of 3 building block barcodes in each sequence read, you would only have 1, 2, or 3 within this column for each row, with each number
representing one of the three building blocks.

## Run
After compilation, the `del` binary can be moved anywhere.
<br>
<br>
Run DEL-Decode<br>

```
del --fastq <fastq_file> \
	--sample_barcodes <sample_barcode_file> \
	--sequence_format <sequence_format_file> \
	--bb_barcodes <building_block_barcode_file> \
	--output_dir <output_dir> \
	--prefix <file_prefix> \
	--threads <num_of_threads> \
	--merge_output
```

<br>
<ul>
<li>
--bb_barcodes is optional.  If it is not used, the output counts uses the DNA barcode to count with no error handling on these barcodes.
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
</ul>

### Output files
Each sample name will get a file in the default format of year-month-day_<sample_name>_counts.csv in the following format (for 3 building blocks):
<table>
<tr>
<th>BB_1</th>
<th>BB_2</th>
<th>BB_3</th>
<th>Count</th>
</tr>
<tr>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>#</td>
</tr>
<tr>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>#</td>
</tr>
</table>

Where BB_ID is used if there is a building block convernsion file, otherwise the DNA code is used. `#` represents the count number<br><br>
If `--merge_output` is called, an additional file is created with the format (for 3 samples):

<table>
<tr>
<th>BB_1</th>
<th>BB_2</th>
<th>BB_3</th>
<th>Sample_1</th>
<th>Sample_2</th>
<th>Sample_3</th>
</tr>
<tr>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>#</td>
<td>#</td>
<td>#</td>
</tr>
<tr>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>BB_ID/DNA code</td>
<td>#</td>
<td>#</td>
<td>#</td>
</tr>
</table>


## Tests results
On an 8 threaded i7-4790K CPU @ 4.00GHz with 16gb RAM, this algorithm was able to decode over 400 million sequencing reads in just under 1 hour and 20 minutes.
Results below:
```
Total sequences: 418770000
Constant Region Mismatches: 173770206
Sample Barcode Mismatches: 1597170
Building Block Mismatches: 4520082
Writing counts
Total time: 78 minutes
```
