#!/usr/bin/env python3
#
# FunnSV
# Fast and effective functional annotation of SVs
#
# Copyright 2022 - Barcelona Supercomputing Center
# Author:  Nicolas Gaitan
# Contact: nicolas.gaitan@bsc.es
# MIT License
import sys
from argparse import ArgumentParser
import logging
from structural_variants_functional_annotator import run_sv_annotation

if not sys.version_info >= (3, 7):
    raise SystemError(
        f'Error: FunnSV works with Python version 3.7 or above \
        (detected version: {sys.version_info.major}.{sys.version_info.minor}). Exiting')


def exec_parser():
    parser = ArgumentParser(
        prog='funnSV',
        description='Fast and effective functional annotation of SVs',
        epilog='')
    parser.add_argument('-i', '--input', type=str, help='input SV vcf file (.vcf, .vcf.gz)', required=True)
    parser.add_argument('-g', '--gff_file', type=str, help='input gff to base annotations on (.gff3, .gff, .gtf)',
                        required=True)
    parser.add_argument('-r', '--ref_genome', type=str,
                        help='input fasta reference genome, must be indexed (.fa, .fasta, .fna)', required=True)
    parser.add_argument('-o', '--output_vcf_prefix', type=str, help='output file path prefix for annotated vcf file',
                        required=True)
    parser.add_argument('-m', '--mode', type=str,
                        choices=['minimal', 'balanced', 'complete', 'minimal_pc'], help='mode of annotations, between '
                                                                                        'minimal: annotate only genes '
                                                                                        ', balanced: annotate genes '
                                                                                        'and their transcripts with '
                                                                                        'specific elements ,'
                                                                                        'complete: annotate every '
                                                                                        'element in the gff - '
                                                                                        'to be used with caution -, '
                                                                                        'minimal_pc: Only protein '
                                                                                        'coding genes annotated with'
                                                                                        'the column 9 field '
                                                                                        'biotype=protein_coding',
                        required=False, default='minimal')
    parser.add_argument('-f', '--fields', type=str, help='specific Gene and/or Transcript fields that are annotated '
                                                         'per variant in the vcf. ALL: write all gff fields or '
                                                         'specify the input gff fields separated by comma, example: '
                                                         'ID,Name,Alias Warning: Fields not found in the gff file '
                                                         'should not be used in this parameter',
                        required=False, default='ALL')
    config = parser.parse_args()
    return config


def funnSV_main():
    logging.basicConfig(level=logging.INFO)
    logging.info('Beginning execution')
    config = exec_parser()
    vcf_path = config.input
    gff_path = config.gff_file
    ref_genome_path = config.ref_genome
    mode = config.mode
    fields = config.fields
    vcf_output = config.output_vcf_prefix
    logging.info('Producing annotation of vcf structural variants')
    run_sv_annotation(vcf_path, gff_path, ref_genome_path, mode, fields, vcf_output)
    logging.info('Finished execution of funnSV successfully')


if __name__ == "__main__":
    funnSV_main()
