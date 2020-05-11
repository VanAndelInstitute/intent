# intent

Convert fastq files from inDrop v2 libraries to 10X genomics style fastqs to facilitate
analysis by some downstream tools like alevin.

This program takes an R1 and and R2 fastq from an inDrop run, and reorganizes the barcodes. 
R1 and R2 are swapped, and the cellular barcodes 1 and 2 are extracted and put at the start 
of each read and the UMI is put after the cellular barcode.  Note that inDrop v2 uses just 
6 bases for their UMIs. This sofware extends the UMIs to 10-14 basepairs by adding some or all of 
cell barcode 2 on to it. By default, a 12 bp UMI is created in this way. This is shortened 
to 10 bp if the -t flag is given (to simulate v2 10X Genomics chemistry), or lengthened 
to 14 if the -x flag is given. In addition, the combined cell barcode for v2 inDrop chemistry 
is 16-20 basepairs. This is truncated to the first 16 bases of BC2+BC1 when run in default 
mode or with the -t flag. If the -x flag is given, a 20bp cell barcode is provided (which 
in some cases will require padding with up to 4 G's since BC1 can be anywhere between 8 
and 12 bp).

If barcode 1 is shorter than 8 bp or if the UMI is shorter than 6bp, the barcode (and 
corresponding transcript read) is discarded. If barcode 1 os 8-11 basepairs, it is padded
to 12 basepairs with G's. If the W1 portion of R2 aligns poorly to the expected W1 sequence 
(`GAGTGATTGCTTGTGACGCCTT`), then the barcode is discarded. "Poorly" is defined as exceeding 
the Levenshtein distance limit set by the parameter `-d`, which by default is set to '2'.

The quality metrics for the reads are preserved and re-ordered along with the barcodes.

The output files are as follows (C = cell barcode, U = UMI):

    output_r1.fastq: 

        default (Chromium v3 style):
        CCCCCCCCCCCCCCCCUUUUUUUUUUUU

        -t (Chromium v2 style):
        CCCCCCCCCCCCCCCCUUUUUUUUUU

        -x (Entire cell barcode, padded if needed, and augmented UMI)
        CCCCCCCCCCCCCCCCCCCCUUUUUUUUUUUUUU

    output_r2.fastq:

        All transcript reads corresponding to high quality barcodes.

Note that all input files must be gzipped, and the program outputs gzipped files.

**CAVEAT EMPTOR** This software is brand spanking new and minimally tested. Please 
verify results as your mileage may vary.

## Installation

Clone the repository, then:

```bash
cd intent
make
cp bin/intent /somewhere/on/your/PATH
```

## Usage

    intent [options] read1.fastq read2.fastq output    

        read1.fastq.gz      read 1 fastq file (gzipped) from a v2 inDrop sequencing run.
        read2.fastq.gz      read 2 fastq file (gzipped) from a v2 inDrop sequencing run.
        output              root name of the output; transformed fastqs will be  
                            written to output_R1.fastq and output_R2.fastq. R1 
                            and R2 are swapped, and ambigious reads dropped  
                            (i.e. those with low quality w1 read or poorly formed
                            barcodes or UMIs); the resulting R1 fastq will have 16
                            bases of cell barcode follwed by 12 of UMI (unless
                            -x or -t flags are specified); the R2 fastq will
                            contain the corresponding transcript reads from the 
                            original R1 fastq (excluding those corresponding to 
                            dropped low quality barcodes. 

        Options:
            -x              expanded barcodes and umi format (20 / 14) 
            -t              Chromium v2 barcode format (16 / 10) 
            -d 2            maximum W1 alignment distance (default = 2)
            -h              display this help message

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Notably, only V2 of inDrop chemistry is supported. Naively, it should be pretty easy to add support for V1 or V3.
I just don't have access to that fastq files nor deep understanding of their barcoding scheme. I would be happy 
to discuss further or receive pull requests.

## License
[MIT](https://choosealicense.com/licenses/mit/)
