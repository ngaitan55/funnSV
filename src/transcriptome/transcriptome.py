from __future__ import annotations
from typing import Optional, List, Dict

# Static values
# Functional common region types
GENE = 'gene'
TRANSCRIPT = 'mRNA'
EXON = 'exon'
CDS = 'CDS'
FIVE_PRIME_UTR = 'five_prime_UTR'
THREE_PRIME_UTR = 'three_prime_UTR'


class FunctionalGenomicRegion:
    ID: str
    sequence_name: str
    first: int
    last: int
    length: int
    region_type: str
    info: Optional[Dict[str, str]]

    def __init__(self, region_id: str, sequence_name: str, first: int, last: int, length: int, region_type: str,
                 info: Optional[Dict[str, str]]):
        self.ID = region_id
        self.sequence_name = sequence_name
        self.first = first
        self.last = last
        self.last = length
        self.region_type = region_type
        self.info = info


class TranscriptElement(FunctionalGenomicRegion):
    transcript: Transcript

    def __init__(self, transcript: Transcript, element_id: Optional[str], first: int, last: int, length: int,
                 region_type: str, info: Optional[Dict[str, str]] = None):
        self.transcript = transcript
        if region_type in (GENE, TRANSCRIPT):
            raise TypeError('Transcript element must be a different type from Gene or Transcript')
        super(TranscriptElement, self).__init__(element_id, transcript.sequence_name, first, last, length, region_type,
                                                info)


class Transcript(FunctionalGenomicRegion):
    gene: Gene
    elements: List[FunctionalGenomicRegion]

    def __init__(self, gene: Gene, transcript_id: str, first: int, last: int, length: int,
                 info: Optional[Dict[str, str]] = None):
        self.gene = gene
        self.elements = []
        super(Transcript, self).__init__(transcript_id, gene.sequence_name, first, last, length, TRANSCRIPT, info)

    def add_element(self, new_element: TranscriptElement):
        self.elements.append(new_element)


class Gene(FunctionalGenomicRegion):
    transcripts: List[Transcript]

    def __init__(self, gene_id: str, sequence_name: str, first: int, last: int, length: int,
                 info: Optional[Dict[str, str]] = None):
        self.transcripts = []
        super(Gene, self).__init__(gene_id, sequence_name, first, last, length, GENE, info)

    def add_transcript(self, new_transcript: Transcript):
        self.transcripts.append(new_transcript)
