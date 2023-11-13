from __future__ import annotations
import sys
#sys.path.append('../')

from typing import List, Union, Sequence, Callable, Dict
import logging
import pysam
from pysam import FastaFile
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantRecord
from variant_extractor.variants import VariantType
from src.funnSV.annotation_sources.transcriptome import Gene, Transcript, TranscriptElement, FunctionalGenomicRegion
from src.funnSV.annotation_sources.annotations_handler import load_transcriptome_from_gff3

FIELDS_ALL_PARAM = 'ALL'

FIELDS_INFO_KEY = 'FUNNSV_ANNOTATIONS'
FIELDS_SEPARATOR = '|'
FIELD_NOT_PRESENT = '-'
OUTPUT_FILE_EXTENSION = 'vcf'


def _load_structural_variants(vcf_path: str, ref_sequences: Sequence[str]) -> List[VariantRecord]:
    variants: List[VariantRecord] = []
    variants_extractor = VariantExtractor(vcf_path)
    for variant_record in variants_extractor:
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
                                  ranked_sequences: Dict[str, int],
                                  fields: str,
                                  sv_length: int,
                                  compare: Callable[[VariantRecord, Union[
                                      Gene, Transcript, TranscriptElement, FunctionalGenomicRegion],
                                                     Dict[str, int]], int]):
    """
    Annotates a structural variant list with an annotation_sources coming from the annotation_sources module
    :param variants: List with structural variants to annotate - <pre> Must be ordered by reference sequence and coordinates </pre>
    :param transcriptome: List with an annotation_sources represented by genes - <pre> Must be ordered by reference sequence and coordinates </pre>
    :param fields: Fields to annotate, if == ID, this function will only annotate genes as a list of IDs
    :param sv_length: length threshold that defines if a variant is an SV
    :param compare: Function that compares a variant with a functional genomic region allowing the
        definition of the comparison to be flexible and dependent on the implementation of the parameter function
        must return an integer type
    """

    def make_single_annotation(dict_x: Dict[str, str]) -> str:
        filtered_items = []
        if fields != FIELDS_ALL_PARAM:
            fields_list = fields.split(',')
            for field in fields_list:
                field_value = dict_x.get(field)
                if field_value is None:
                    logging.warning(f"expected field {field} is not present in"
                                    f" record {dict_x['ID']}")
                    field_value = FIELD_NOT_PRESENT
                filtered_items.append(field_value)
        else:
            filtered_items = dict_x.values()
        return FIELDS_SEPARATOR.join(filtered_items)

    def is_indel(variant: VariantRecord) -> bool:
        return (variant.variant_type == VariantType.INS or
                variant.variant_type == VariantType.DEL) and variant.length < sv_length

    i = 0
    j = 0
    j_lower_bound = 0
    variant_interacted = False
    gene_annotations: List[str] = []
    while i < len(variants) and j < len(transcriptome):
        current_variant = variants[i]
        if current_variant.variant_type == VariantType.SNV or is_indel(current_variant):
            i = i + 1
            continue
        current_gene = transcriptome[j]
        cmp: int = compare(current_variant, current_gene, ranked_sequences)
        # print(
        #     f'VARIANT: seq={current_variant.contig} type={current_variant.variant_type} first={current_variant.pos} last={current_variant.end} - GFFRecord: seq={current_gene.sequence_name}'
        #     f' first={current_gene.first} last={current_gene.last} - COMPARE={cmp}')
        if cmp < -1:
            i = i + 1
            j = j_lower_bound
            variant_interacted = False
            if gene_annotations:
                annotation = f"{','.join(gene_annotations)}"
                current_variant.info[FIELDS_INFO_KEY] = annotation
                gene_annotations = []
        elif cmp > 1:
            j_lower_bound = j_lower_bound + 1
            j = j_lower_bound
        else:
            if not variant_interacted:
                variant_interacted = True
                j_lower_bound = j
            gene_annotations.append(make_single_annotation(current_gene.info))
            transcripts = current_gene.child_elements
            if transcripts:
                for transcript in transcripts:
                    cmp_tr = compare(current_variant, transcript, ranked_sequences)
                    if -1 <= cmp_tr <= 1:
                        gene_annotations.append(make_single_annotation(transcript.info))
                        elements = transcript.child_elements
                        if elements:
                            for element in elements:
                                cmp_elem = compare(current_variant, element, ranked_sequences)
                                if -1 <= cmp_elem <= 1:
                                    gene_annotations.append(f"{element.region_type}{FIELDS_SEPARATOR}"
                                                            f"{element.info['Parent']}")
            if cmp == 1:
                j = j + 1
            else:
                i = i + 1
                j = j_lower_bound
                variant_interacted = False
                if gene_annotations:
                    annotation = f"{','.join(gene_annotations)}"
                    current_variant.info[FIELDS_INFO_KEY] = annotation
                    gene_annotations = []


def _assign_index_ordered_sequences(ref_sequences: Sequence[str]) -> Dict[str, int]:
    return {k: v for (k, v) in zip(ref_sequences, range(len(ref_sequences)))}


def _extract_header(vcf_path: str, fields_header_info_field: str) -> str:
    with open(vcf_path, 'r') as reader:
        with pysam.VariantFile(reader) as vcf_file:
            vcf_file.header.add_line(fields_header_info_field)
            header_str = str(vcf_file.header)
    return header_str


def run_sv_annotation(vcf_path: str, gff_path: str, ref_genome_path: str, mode: str, fields: str, sv_length: int,
                      vcf_output: str):
    try:
        fields_header_info_field = (f"##INFO=<ID={FIELDS_INFO_KEY},Number=.,Type=String,Description=Annotated fields "
                                    f"per gene: {fields.replace(',', FIELDS_SEPARATOR)}>")
        header_lines: str = _extract_header(vcf_path, fields_header_info_field)
        ref_genome: FastaFile = FastaFile(ref_genome_path)
        ref_sequences: Sequence[str] = ref_genome.references
        indexed_ref_sequences = _assign_index_ordered_sequences(ref_sequences)
        variants: List[VariantRecord] = _load_structural_variants(vcf_path, ref_sequences)
        transcriptome, _ = load_transcriptome_from_gff3(gff_path, ref_sequences, mode)
        variants.sort(key=lambda x: (indexed_ref_sequences[x.contig], x.pos))
        transcriptome.sort(key=lambda y: (indexed_ref_sequences[y.sequence_name], y.first))
        _annotate_structural_variants(variants, transcriptome, indexed_ref_sequences, fields, sv_length,
                                      _compare_by_coordinates)
        with open('.'.join([vcf_output, OUTPUT_FILE_EXTENSION]), 'w') as writer:
            writer.write(header_lines)
            for annotated_variant in variants:
                writer.write(str(annotated_variant) + '\n')
    except Exception as error:
        logging.error(error)
        raise
