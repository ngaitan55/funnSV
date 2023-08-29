from __future__ import annotations
from typing import Optional, List, Dict, Union

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
    child_elements: List[Union[FunctionalGenomicRegion, Transcript, TranscriptElement]]

    def __init__(self, region_id: str, sequence_name: str, first: int, last: int, length: int, region_type: str,
                 info: Optional[Dict[str, str]]):
        self.ID = region_id
        self.sequence_name = sequence_name
        self.first = first
        self.last = last
        self.length = length
        self.region_type = region_type
        self.info = info
        self.child_elements = []

    def add_element(self, new_element: Union[FunctionalGenomicRegion, Transcript, TranscriptElement]):
        self.child_elements.append(new_element)

    def __str__(self):
        return f"{self.ID}'\t'{self.sequence_name}'\t'{self.region_type}'\t'{self.first}'\t'" \
               f"{self.last}'\t'{self.length}'\t'{self.info}"


class TranscriptElement(FunctionalGenomicRegion):

    def __init__(self, element_id: Optional[str], sequence_name: str, first: int, last: int, length: int,
                 region_type: str, info: Optional[Dict[str, str]] = None):
        if region_type in (GENE, TRANSCRIPT):
            raise TypeError('Transcript element must be a different type from Gene or Transcript')
        super(TranscriptElement, self).__init__(element_id, sequence_name, first, last, length, region_type,
                                                info)


class Transcript(FunctionalGenomicRegion):

    def __init__(self, transcript_id: str, sequence_name: str, first: int, last: int, length: int,
                 info: Optional[Dict[str, str]] = None):
        super(Transcript, self).__init__(transcript_id, sequence_name, first, last, length, TRANSCRIPT, info)


class Gene(FunctionalGenomicRegion):

    def __init__(self, gene_id: str, sequence_name: str, first: int, last: int, length: int,
                 info: Optional[Dict[str, str]] = None):
        super(Gene, self).__init__(gene_id, sequence_name, first, last, length, GENE, info)
