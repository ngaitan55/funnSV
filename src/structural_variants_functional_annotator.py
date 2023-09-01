from __future__ import annotations

from typing import List, Optional, Union, Sequence, Callable, Dict
import logging
import pysam
from pysam import FastaFile
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantRecord
from transcriptome.transcriptome import Gene, Transcript, TranscriptElement, FunctionalGenomicRegion, GENE, TRANSCRIPT, \
    EXON, CDS, THREE_PRIME_UTR, FIVE_PRIME_UTR
from transcriptome.annotations_handler import load_transcriptome_from_gff3, LOAD_MODE_MINIMAL, LOAD_MODE_BALANCED, \
    LOAD_MODE_COMPLETE

OUTPUT_FILE_EXTENSION = 'vcf'


def _load_structural_variants(vcf_path: str, ref_sequences: Sequence[str]) -> List[VariantRecord]:
    variants: List[VariantRecord] = []
    with VariantExtractor(vcf_path) as variants_iterator:
        for variant_record in variants_iterator:
            if variant_record.contig not in ref_sequences:
                raise ValueError(
                    f'variant record located in invalid sequence: {variant_record.contig}, check vcf and reference genome')
            variants.append(variant_record)
    return variants


def _compare_by_coordinates(variant: VariantRecord,
                            functional_region: Union[Gene, Transcript, TranscriptElement, FunctionalGenomicRegion],
                            sequence_order: Dict[str, int]) -> int:
    """
    Custom comparator for functional regions and variants according to genomic coordinates
    :param variant
    :param functional_region
    :param sequence_order: dictionary that returns the ordered index according to the str sequence
    :return: 0 the regions are exactly the same
             -3 or 3 they are located in different sequences
             -2 or 2 they are located in the same sequence but do not overlap
             -1 or 1 they overlap
             negative and positive shows the position of the variant according to the functional_region
    """

    def _compare(seq_idx1: int, first1: int, last1: int, seq_idx2: int, first2: int, last2: int) -> int:
        overlap = first2 <= last1 and last2 >= first1
        if seq_idx1 < seq_idx2:
            return -3
        if seq_idx1 > seq_idx2:
            return 3
        # For these cases seq_idx1 == seq_idx2
        if last1 < last2:
            return -1 if overlap else -2
        if last2 < last1:
            return 1 if overlap else 2
        # For these cases last1 == last2
        if first1 < first2:
            return -1
        if first2 < first1:
            return 1
        return 0

    return _compare(sequence_order[variant.contig], variant.pos, variant.end,
                    sequence_order[functional_region.sequence_name], functional_region.first, functional_region.last)


# def _compute_regions_intersection(first1: int, last1: int, first2: int, last2: int) -> float:
#   if not _intersects(first1, last1, first2, last2):
#       return float(0)
#   intersected_bases: int
#   if first1 >= first2 and last1 >= last2:
#       intersected_bases = last2 - first1
#   elif last2 > last1 and first2 > first1:
#       intersected_bases = last1 - first2
#   else:
#       intersected_bases = min(last1 - first1, last2 - first2)
#   denominator = max(last1 - first1, last2 - first2)
#   return float(intersected_bases / denominator)


def _annotate_structural_variants(variants: List[VariantRecord], transcriptome: List[Gene, FunctionalGenomicRegion],
                                  compare: Callable[[VariantRecord, Union[
                                      Gene, Transcript, TranscriptElement, FunctionalGenomicRegion]], int]):
    """
    Annotates a structural variant list with a transcriptome coming from the transcriptome module
    :param variants: List with structural variants to annotate - <pre> Must be ordered by reference sequence and coordinates </pre>
    :param transcriptome: List with a transcriptome represented by genes - <pre> Must be ordered by reference sequence and coordinates </pre>
    :param compare: Function that compares a variant with a functional genomic region allowing the
        definition of the comparison to be flexible and dependent on the implementation of the parameter function
        must return an integer type
    """
    gene_annotation_key_prefix = 'GENE_ANNOTATION'
    transcript_annotation_key_prefix = 'TRANSCRIPT_ANNOTATION'
    element_annotation_key_prefix = 'ELEMENT_ANNOTATION'

    def single_annotation(dict_x: Dict[str, str]) -> str:
        return ';'.join([f'{key}={value}' for key, value in dict_x.items()])

    variants_it = iter(variants)
    transcriptome_it = iter(transcriptome)
    left_over_genes: List[Gene] = []
    current_variant: VariantRecord = next(variants_it, None)
    current_gene: Gene = next(transcriptome_it, None)
    while current_variant is not None and current_gene is not None:
        cmp: int = compare(current_variant, current_gene)
        if cmp < -1:
            current_variant = next(variants_it)
        if cmp > 1:
            current_gene = next(transcriptome_it, None)
        else:
            gene_annotation_key = f'{gene_annotation_key_prefix}_{current_gene.ID}'
            gene_annotation_value = f'[{single_annotation(current_gene.info)}]'
            current_variant.info[gene_annotation_key] = gene_annotation_value
            transcripts = current_gene.child_elements
            if transcripts:
                for transcript in transcripts:
                    cmp_tr = compare(current_variant, transcript)
                    if -1 <= cmp_tr <= 1:
                        transcript_annotation_key = f'{transcript_annotation_key_prefix}_{transcript.ID}'
                        transcript_annotation_value = f'[{single_annotation(transcript.info)}]'
                        current_variant.info[transcript_annotation_key] = transcript_annotation_value
                        elements = transcript.child_elements
                        if elements:
                            for element in elements:
                                cmp_elem = compare(current_variant, element)
                                if -1 <= cmp_elem <= 1:
                                    element_annotation_key = f'{element_annotation_key_prefix}_{element.region_type}'
                                    element_annotation_value = f'[{single_annotation(element.info)}]'
                                    current_variant.info[element_annotation_key] = element_annotation_value
            if cmp == 1:
                current_gene = next(transcriptome_it, None)
            if cmp == -1:
                current_variant = next(variants_it, None)
            else:
                current_gene = next(transcriptome_it, None)
                current_variant = next(variants_it, None)


def _assign_index_ordered_sequences(ref_sequences: Sequence[str]) -> Dict[str, int]:
    return {k: v for (k, v) in zip(ref_sequences, range(len(ref_sequences)))}


def _extract_header(vcf_path: str) -> str:
    with open(vcf_path, 'r') as reader:
        with pysam.VariantFile(reader) as vcf_file:
            header_str = str(vcf_file.header)
    return header_str


def run_sv_annotation(vcf_path: str, gff_path: str, ref_genome_path: str, mode: str, vcf_output: str):
    try:
        header_lines: str = _extract_header(vcf_path)
        ref_genome: FastaFile = FastaFile(ref_genome_path)
        ref_sequences: Sequence[str] = ref_genome.references
        indexed_ref_sequences = _assign_index_ordered_sequences(ref_sequences)
        variants: List[VariantRecord] = _load_structural_variants(vcf_path, ref_sequences)
        transcriptome: List[Gene, FunctionalGenomicRegion] = load_transcriptome_from_gff3(gff_path, ref_sequences, mode)
        variants.sort(key=lambda x: (indexed_ref_sequences[x.contig], x.pos))
        transcriptome.sort(key=lambda y: (indexed_ref_sequences[y.sequence_name], y.first))
        _annotate_structural_variants(variants, transcriptome, _compare_by_coordinates)
        with open(''.join([vcf_output, OUTPUT_FILE_EXTENSION])) as writer:
            writer.write(header_lines)
            for annotated_variant in variants:
                writer.write(str(annotated_variant) + '\n')
    except Exception as error:
        logging.error(error)
        raise


# Testing purposes only
if __name__ == "__main__":
    import sys

    ref = sys.argv[1]
    ref_gen: FastaFile = FastaFile(ref)
    sequences: Sequence[str] = ref_gen.references
    for seq in sequences:
        print(seq)
