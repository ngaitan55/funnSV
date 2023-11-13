# FunnSV
Efficient and reliable annotation of Structural Variants (SVs) from a VCF file and a GFF3 annotation source.

## Installation

FunnSV can be installed using pip:

`pip install funnSV`

## Requirements

- Python >= 3.7
- pysam
- variant-extractor

#### Tested on:

- Python 3.10
- pysam 0.21.0

## Common use cases
### Gene only annotation
To annotate gene elements from a GFF3 file in the form of their fields Name and Symbol:

`funnSV -i <INPUT_VCF> -g <INPUT_GFF3> -r <REFERENCE_GENOME> -o <OUTPUT_VCF> -f ID,Name`

### Gene, Transcript, and transcript element annotations
To annotate all possible common elements in a GFF3 file, including all fields:

`funnSV -i <INPUT_VCF> -g <INPUT_GFF3> -r <REFERENCE_GENOME> -o <OUTPUT_VCF> -m balanced`

### All types annotations
This mode is not recommended but can be useful if there are pseudogenes or other uncommon types, but it is prone to errors due to inconsistencies in GFF3 files:

`funnSV -i <INPUT_VCF> -g <INPUT_GFF3> -r <REFERENCE_GENOME> -o <OUTPUT_VCF> -m complete -f ID,Name`

### Specific use case: Protein coding genes

This mode allows to only annotate protein coding genes, if gene records contain the `biotype=protein_coding` field in the attributes column (This has been tested using the Human GRCh37 Release 87 GFF3):

`funnSV -i <INPUT_VCF> -g <INPUT_GFF3> -r <REFERENCE_GENOME> -o <OUTPUT_VCF> -m minimal_pc -f ID,Name`

Finally, this software was implemented by using the GFF3 specification from https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md. Possible errors or bugs coming from deviations in the format are to be expected. Please refer any such occurrances to the GitHub issues section. We are open to support other annotation source formats in the future.