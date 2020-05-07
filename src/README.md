# intent

Convert fastq files from inDrop v2 libraries to 10X genomics style fastqs to facilitate
analysis by some downstream tools like alevin.

This program takes an R1 and and R2 fastq from an inDrop run, and reorganizes the barcodes. 
R1 and R2 are swapped, and the cellular barcodes 1 and 2 are extracted and put at the start 
of each read (the first 20 bases), and the UMI is put after the cellular barcode (the next 
14 bases). Note that inDrop v2 uses just 6 bases for their UMIs which is not enough to 
avoid collisions if the cell identity is ignored. So to help avoid this, I append Barcode 2 
to the UMI (that is why the UMI is 14 bases and not 6). 

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
    
    read1.fastq        read 1 fastq file from inDrop sequencing run.
    read2.fastq        read 2 fastq file from inDrop sequencing run.
    output             root name of the output; transformed fastqs will be  
                        written to output_R1.fastq and output_R2.fastq. R1 
                        and R2 are swapped, and ambigious reads dropped  
                        (i.e. those with low quality w1 read or poorly formed
                        barcodes or UMI)
    Options:
        -h               Print this help

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Notably, only V2 of inDrop chemistry is supported. Naively, it should be pretty easy to add support for V1 or V3.
I just don't have access to that fastq files nor deep understanding of their barcoding scheme. I would be happy 
to discuss further or receive pull requests.

## License
[MIT](https://choosealicense.com/licenses/mit/)
