# DEL-Decode
DNA encoded library decoding.  Multithreaded and low resource use.  Uses one thread to read and the rest to process the data, so at least a 2 threaded machine is essential.
This program does not store all data within RAM but instead sequentially processes the sequencing data in order to remain memory efficient.  
With very large DEL libraries, there may be some memory creep.  This is because the more barcodes that exist within the sequencing reads, the more barcodes that need to remain in memory.
<br>
<br>
Error handling is set at 20% maximum sequence error.  The algorithm fixes any sequenced constant region or barcode with the best match possoble.  If there are two or more best matches,
it is not counted.

## Table of Contents
<ul>
<li><a href=#Requirements>Requirements</a></li>
<li><a href=#files-needed>Files Needed</a></li>
<li><a href=#run>Run</a></li>
<li><a href=#test-results>Test Results</a></li>
</ul>

## Requirements

### Rust installed locally: <a href=https://www.rust-lang.org/tools/install>instructions here</a>

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
Currently only accepts unzipped FASTQ files due to this program reading line by line.<br>
Currently tryiing to implement a gzip inflate stream so that the whole file is not placed in RAM.

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
The sample_barcode_file is a comma separate file with the following format:<br>
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
and the last column is the building block number as an integer.  These are in the same order as the sequence format file and starting
at 1.

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
	--threads <num_of_threads>
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
--threads defaults to the number of threads on the machine if not used.
</li>
</ul>


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
