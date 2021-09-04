# DEL-Decode
DNA encoded library decoding

# WORK IN PROGRESS

## Table of Contents
<ul>
<li><a href=#Requirements>Requirements</a></li>
<li><a href=#setup>Setup</a></li>
</ul>

## Requirements
<ul>
<li>Rust install</li>
</ul>

## Setup
`$ cargo build --release`<br>
`$ ./target/release/del --fastq <fastq_file> --sample_barcodes <sample_barcode_file> --sequence_format <sequence_format_file>`<br>
Where the fastq file is unzipped (still working on getting it to work with gzip files) <br>
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

And sequence_format_file is in the format of the example.
