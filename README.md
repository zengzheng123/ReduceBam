# ReduceBam

##### intersect a BAM file with multiple VCF and/or MAF files, extract informative alignments

### Install:
Run "make" in the program directory to compile

### Usage:
```
[ARGUMENTS]

--input_variant     <string>      Input variant file, support .vcf/.vcf.gz/.maf/.maf.gz format, can be specified multiple times.
--input_bam         <string>      Input BAM file
--output_bam        <string>      Output BAM file
--output_bed        <string>      Optional output BED3 file computed from the input VCF/MAF file and the buffer
--buffer            <int>         Buffer size to be added to each side of the VCF/MAF entry. Default is 0
--help                            Print command line usage
```

### Notes:
```
The optional output bed3 file specified by --output_bed is the actual interval regions used to intersect with the BAM file,
computed by merging all the input varaint files plus the buffer region, all interval regions are limited to [1, CHROM_LEN]

For SNP in VCF file, the region is defined as [POS, POS + len(REF) - 1], plus the buffer region on each side
For INS in VCF file, the region is defined as [POS, POS + len(ALT) - 1], plus the buffer region on each side
For DEL in VCF file, the region is defined as [POS, POS + len(REF) - 1], plus the buffer region on each side

For SNP in MAF file, the region is defined as [START_POS, END_POS], plus the buffer region on each side
For INS in MAF file, the region is defined as [START_POS, END_POS + len(ALT) - 1], plus the buffer region on each side
For DEL in MAF file, the region is defined as [START_POS - 1, END_POS], plus the buffer region on each side
```

### FAQ:

```






``` 

### This software uses the following library

bamtools https://github.com/pezmaster31/bamtools

gzstream http://www.cs.unc.edu/Research/compgeom/gzstream/

