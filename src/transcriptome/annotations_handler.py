from __future__ import annotations
import sys
from typing import Dict, Union, Optional
from transcriptome import FunctionalGenomicRegion, Gene, Transcript, TranscriptElement, GENE, TRANSCRIPT, EXON, \
    THREE_PRIME_UTR, FIVE_PRIME_UTR, OTHER_ELEMENT_TYPE

# Static values
GFF_FIELD_SEPARATOR: str = '\t'
GFF_COLUMN_NUMBER: int = 9

ATTRIBUTE_SEPARATOR: str = ';'
ATTRIBUTE_KEY_VALUE: str = '='


class GFF3FileReader:
    # Static value to check GFF format line
    _format_line_id = '#'
    # Attributes
    format_lines = []

    def __init__(self, file_path: str):
        if not file_path.endswith(('.gff', '.gff3', '.gtf')):
            raise NameError('The file extension is not correct, check if it is a gff formatted file')
        self.__path = file_path
        self.__file_object = None

    def __enter__(self):
        self.__file_object = open(self.__path, 'r')
        first_header_line: str = self.__file_object.readline()
        if not first_header_line.startswith(self._format_line_id):
            raise ValueError('gff3 file header is malformed, minimal header line not found')
        else:
            self.format_lines.append(first_header_line)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.__file_object.close()

    def __iter__(self):
        return self

    def __next__(self) -> Optional[GFF3Record]:
        gff_file_line: str = self.__file_object.readline()
        if self.__file_object is None or gff_file_line == '':
            raise StopIteration
        elif gff_file_line.startswith(self._format_line_id):
            self.format_lines.append(gff_file_line)
            return
        else:
            gff_line_elements = gff_file_line.split(GFF_FIELD_SEPARATOR)
            if len(gff_line_elements) != GFF_COLUMN_NUMBER:
                raise ValueError(f'gff3 line is malformed, it contains {len(gff_line_elements)} not 9 columns')
            seqid: str = gff_line_elements[0].strip()
            source: str = gff_line_elements[1].strip()
            genomic_region_type: str = gff_line_elements[2].strip()
            start: int = int(gff_line_elements[3].strip())
            end: int = int(gff_line_elements[4].strip())
            raw_score = gff_line_elements[5].strip()
            score: Union[str, float] = float(raw_score) if raw_score != '.' else '.'
            strand: str = gff_line_elements[6].strip()
            raw_phase = gff_line_elements[7].strip()
            phase: Union[str, int] = int(raw_phase) if raw_phase != '.' else '.'
            info_string: str = gff_line_elements[8].strip()
            return GFF3Record(seqid, source, genomic_region_type, start, end, score, strand, phase, info_string)


# Informal Interface for record classes from different sources
class AnnotationRecord:
    def make_functional_annotation(self, *args) -> FunctionalGenomicRegion:
        pass


def _compute_attributes_from_gff_info_field(attributes_str: str) -> Dict[str, str]:
    attributes = {}
    separated = attributes_str.split(ATTRIBUTE_SEPARATOR)
    for s in separated:
        k, v = s.split(ATTRIBUTE_KEY_VALUE)
        attributes[k] = v
    if len(attributes) == 0:
        raise ValueError('attributes were not consumed for current record, check column 9 of GFF file')
    return attributes


class GFF3Record(AnnotationRecord):
    # Class that represents a typical GFF3 record based on the format specification by Lincoln Stein at https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    ATTRIBUTE_ID = 'ID'
    ATTRIBUTE_NAME = 'Name'
    ATTRIBUTE_ALIAS = 'Alias'

    # Fields
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: Union[str, float]
    strand: str
    phase: Union[str, int]
    attributes: Dict[str, str]

    def __init__(self, seqid: str, source: str, param_type: str, start: int, end: int, score: Union[str, float],
                 strand: str, phase: Union[str, int], attributes_str: str):
        self.seqid = seqid
        self.source = source
        self.type = param_type
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = _compute_attributes_from_gff_info_field(attributes_str)

    def make_functional_annotation(self, parent: Optional[Gene, Transcript] = None) -> FunctionalGenomicRegion:
        ID = self.attributes.get(self.ATTRIBUTE_ID, None)
        sequence_name = self.seqid
        first = self.start
        last = self.end
        length = last - first + 1
        info = self.attributes
        if GENE == self.type:
            if ID is None:
                raise ValueError('gff3 Gene record is malformed, it does not contain an ID')
            return Gene(ID, sequence_name, first, last, length, info)
        elif TRANSCRIPT == self.type:
            if ID is None:
                raise ValueError('gff3 Transcript mRNA record is malformed, it does not contain an ID')
            if parent is None:
                raise ValueError('Transcript has no parent Gene, it cannot be instantiated')
            return Transcript(parent, ID, first, last, length, info)
        else:
            if parent is None:
                raise ValueError('Transcript element has no parent Transcript, it cannot be instantiated')
            return TranscriptElement(parent, ID, first, last, length, self.type, info)

    def __str__(self):
        return f"{self.seqid}{GFF_FIELD_SEPARATOR}{self.source}{GFF_FIELD_SEPARATOR}{self.type}{GFF_FIELD_SEPARATOR}{self.start}{GFF_FIELD_SEPARATOR}" \
               f"{self.end}{GFF_FIELD_SEPARATOR}{self.score}{GFF_FIELD_SEPARATOR}{self.strand}{GFF_FIELD_SEPARATOR}{self.phase}{GFF_FIELD_SEPARATOR}{self.attributes}"


if __name__ == "__main__":
    test_file = sys.argv[1]
    with GFF3FileReader(test_file) as gff_reader:
        for gff_record in gff_reader:
            print(str(gff_record))
