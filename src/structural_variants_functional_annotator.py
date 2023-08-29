from __future__ import annotations
from typing import List, Optional, Union, Sequence
import logging
import pysam
from pysam import FastaFile
from variant_extractor import VariantExtractor
from variant_extractor.variants import VariantRecord
from transcriptome.transcriptome import Gene, Transcript, TranscriptElement, FunctionalGenomicRegion, GENE, TRANSCRIPT, \
    EXON, CDS, THREE_PRIME_UTR, FIVE_PRIME_UTR
from transcriptome.annotations_handler import load_transcriptome_from_gff3, LOAD_MODE_MINIMAL, LOAD_MODE_BALANCED, \
    LOAD_MODE_COMPLETE


def _load_structural_variants(vcf_path: str, ref_sequences: Sequence[str]) -> List[VariantRecord]:
    variants: List[VariantRecord] = []
    with VariantExtractor(vcf_path) as variants_iterator:
        for variant_record in variants_iterator:
            if variant_record.contig not in ref_sequences:
                raise ValueError(
                    f'variant record located in invalid sequence: {variant_record.contig}, check vcf and reference genome')
            variants.append(variant_record)
    return variants


def _compute_regions_intersection(first1: int, last1: int, first2: int, last2: int) -> float:
    intersects = first2 <= last1 and last2 >= first1
    if not intersects:
        return float(0)
    intersected_bases: int
    if first1 >= first2 and last1 >= last2:
        intersected_bases = last2 - first1
    elif last2 > last1 and first2 > first1:
        intersected_bases = last1 - first2
    else:
        intersected_bases = min(last1 - first1, last2 - first2)
    denominator = max(last1 - first1, last2 - first2)
    return float(intersected_bases / denominator)


def _compute_regions_intersection(variant: VariantRecord, functional_region: Union[Gene, Transcript, TranscriptElement, FunctionalGenomicRegion]) -> float:
    return _compute_regions_intersection(variant.pos, variant.end, functional_region.first, functional_region.last)


def _annotate_structural_variants(vcf_path: str, gff_path: str, ref_genome_path: str, mode: str) -> List[VariantRecord]:
    annotated_variants: List[VariantRecord] = []
    ref_genome: FastaFile = FastaFile(ref_genome_path)
    ref_sequences: Sequence[str] = ref_genome.references
    variants: List[VariantRecord] = _load_structural_variants(vcf_path)
    transcriptome: List[Gene, FunctionalGenomicRegion] = load_transcriptome_from_gff3(gff_path, ref_sequences, mode)

    return annotated_variants
