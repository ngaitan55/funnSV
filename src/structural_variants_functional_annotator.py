from typing import List, Optional, Union
from pysam import VariantRecord
from variant_extractor import VariantExtractor
from transcriptome.transcriptome import Gene, Transcript, TranscriptElement, FunctionalGenomicRegion, GENE, TRANSCRIPT, \
    EXON, CDS, THREE_PRIME_UTR, FIVE_PRIME_UTR
from transcriptome.annotations_handler import GFF3FileReader, GFF3Record

# Depth of analysis static fields
MODE_INCLUDE_GENES_ONLY = 'genes'
MODE_INCLUDE_ALL_PROTEIN_CODING = 'balanced'
MODE_INCLUDE_ALL = 'all'


def _read_variants(vcf_path: str) -> List[VariantRecord]:
    variants: List[VariantRecord] = []
    variants_iterator = VariantExtractor(vcf_path)
    for variant_record in variants_iterator:
        variants.append(variant_record)
    return variants


def _read_gff3_to_transcriptome(gff3_path: str, mode: str) -> List[Union[Gene, FunctionalGenomicRegion]]:
    transcriptome: List[Gene] = []
    gff_records_iterator = GFF3FileReader(gff3_path)
    current_gene: Optional[Gene, FunctionalGenomicRegion] = None
    current_transcript: Optional[Transcript] = None
    for gff_record in gff_records_iterator:
        if gff_record is None:
            continue
        if gff_record.type == GENE:
            current_gene = gff_record.make_functional_annotation()
            transcriptome.append(current_gene)
        elif gff_record.type == TRANSCRIPT:
            if MODE_INCLUDE_GENES_ONLY == mode:
                continue
            current_transcript = gff_record.make_functional_annotation(current_gene)
            current_gene.add_transcript(current_transcript)
        elif gff_record.type in (EXON, CDS, FIVE_PRIME_UTR, THREE_PRIME_UTR):
            if MODE_INCLUDE_GENES_ONLY == mode:
                continue
            current_transcript_element = gff_record.make_functional_annotation(current_transcript)
            current_transcript.add_element(current_transcript_element)
        else:
            if MODE_INCLUDE_GENES_ONLY == mode or MODE_INCLUDE_ALL_PROTEIN_CODING == mode:
                continue
            current_other_annotation = gff_record.make_functional_annotation()
            transcriptome.append(current_other_annotation)
    return transcriptome
