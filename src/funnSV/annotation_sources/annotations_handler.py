from __future__ import annotations
from typing import Dict, Union, Optional, List, Sequence
from src.funnSV.annotation_sources.transcriptome import FunctionalGenomicRegion, Gene, Transcript, TranscriptElement, GENE, TRANSCRIPT, EXON, \
    THREE_PRIME_UTR, FIVE_PRIME_UTR, CDS

# Static values
GFF_FIELD_SEPARATOR: str = '\t'
GFF_COLUMN_NUMBER: int = 9

ATTRIBUTE_SEPARATOR: str = ';'
ATTRIBUTE_KEY_VALUE: str = '='

LOAD_MODE_MINIMAL = 'minimal'
LOAD_MODE_BALANCED = 'balanced'
LOAD_MODE_COMPLETE = 'complete'
LOAD_MODE_PROTEIN_CODING = 'minimal_pc'


def load_transcriptome_from_gff3(gff3_path: str, ref_sequences: Sequence[str], mode: str) -> (List[Union[Gene, FunctionalGenomicRegion]], List[str]):
    """Load annotation_sources from a gff3 file into a List of Gene objects, depending on the mode, different features will be included, as:
    @:param
        gff3_path:str - the string path for the gff3 file
        ref_sequences: Sequence[str] - Contains valid sequences for this annotation_sources from the reference genome
        mode: str - loading mode, as follows:
            LOAD_MODE_MINIMAL: Only gene records will be loaded
            LOAD_MODE_BALANCED: Genes, mRNA and 5UTR, 3UTR, CDS, Exons will be loaded
            LOAD_MODE_COMPLETE: Elements other than balanced mode objects will be loaded as FunctionalGenomicRegion generic instances
    @:returns List of Gene and FunctionalGenomicRegions representing the annotation_sources and a List with the header lines"""
    transcriptome: List[Union[Gene, FunctionalGenomicRegion]] = []
    agenda: Dict[str, Union[Gene, Transcript, FunctionalGenomicRegion]] = {}
    with GFF3FileReader(gff3_path) as gff_records_iterator:
        for gff_rec in gff_records_iterator:
            if gff_rec is None:
                continue
            if gff_rec.seqid not in ref_sequences:
                raise ValueError(
                    f'gff record is annotated on invalid sequence: {gff_rec.seqid}, check gff file and reference genome')
            if gff_rec.type == GENE:
                current_gene: Gene = gff_rec.make_functional_annotation()
                if LOAD_MODE_PROTEIN_CODING == mode:
                    byotype: str = current_gene.info.get('biotype')
                    if byotype is None or byotype != 'protein_coding':
                        continue
                agenda[current_gene.ID] = current_gene
                transcriptome.append(current_gene)
            elif gff_rec.type == TRANSCRIPT:
                if LOAD_MODE_MINIMAL == mode or LOAD_MODE_PROTEIN_CODING == mode:
                    continue
                current_transcript: Transcript = gff_rec.make_functional_annotation()
                agenda[current_transcript.ID] = current_transcript
                parent_id = current_transcript.info.get(gff_rec.ATTRIBUTE_PARENT)
                parent_gene = agenda.get(parent_id)
                if parent_gene is None:
                    continue
                parent_gene.add_element(current_transcript)
            elif gff_rec.type in (EXON, CDS, FIVE_PRIME_UTR, THREE_PRIME_UTR):
                if LOAD_MODE_MINIMAL == mode or LOAD_MODE_PROTEIN_CODING == mode:
                    continue
                current_transcript_element: TranscriptElement = gff_rec.make_functional_annotation()
                parent_id = current_transcript_element.info.get(gff_rec.ATTRIBUTE_PARENT)
                parent_transcript = agenda.get(parent_id, None)
                if parent_transcript is None:
                    continue
                parent_transcript.add_element(current_transcript_element)
            else:
                if mode != LOAD_MODE_COMPLETE:
                    continue
                current_other_annotation: FunctionalGenomicRegion = gff_rec.make_functional_annotation()
                parent_id = current_other_annotation.info.get(gff_rec.ATTRIBUTE_PARENT, None)
                if current_other_annotation.ID != gff_rec.GENERIC_ID:
                    if parent_id is None:
                        transcriptome.append(current_other_annotation)
                        agenda[current_other_annotation.ID] = current_other_annotation
                    else:
                        parent = agenda.get(parent_id)
                        parent.add_element(current_other_annotation)
        format_lines: List[str] = gff_records_iterator.format_lines
    return transcriptome, format_lines


class GFF3FileReader:
    """File Reader class for standard GFF3 files, requires a file path string to initialize"""
    # Static value to check GFF format line
    _format_line_id = '#'
    # Attributes
    format_lines: List[str] = []

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
        """ Returns GFF3Record if the line is formatted correctly or None if it is a header line"""
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


class AnnotationRecord:
    """Informal Interface for record classes from different sources"""

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
    """Class that represents a typical GFF3 record based on the format specification by Lincoln Stein at
    https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md"""
    # Static attributes
    ATTRIBUTE_ID = 'ID'
    ATTRIBUTE_NAME = 'Name'
    ATTRIBUTE_ALIAS = 'Alias'
    ATTRIBUTE_PARENT = 'Parent'

    GENERIC_ID = 'generic_id'
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

    def make_functional_annotation(self) -> Union[Gene, Transcript, TranscriptElement, FunctionalGenomicRegion]:
        ID = self.attributes.get(self.ATTRIBUTE_ID, None)
        parent = self.attributes.get(self.ATTRIBUTE_PARENT, None)
        sequence_name = self.seqid
        first = self.start
        last = self.end
        length = last - first + 1
        info: Dict[str, str] = self.attributes
        if GENE == self.type:
            if ID is None:
                raise ValueError(f'gff3 Gene record {str(self)} is malformed, it does not contain an ID')
            # Specific line for ensembl gff3 non_collateral fix
            info[self.ATTRIBUTE_ID] = info[self.ATTRIBUTE_ID].replace('gene:', '')
            return Gene(ID, sequence_name, first, last, length, info)
        elif TRANSCRIPT == self.type:
            if ID is None:
                raise ValueError(f'gff3 Transcript mRNA record {str(self)} is malformed, it does not contain an ID')
            if parent is None:
                raise ValueError(f'Transcript from {str(self)} has no parent Gene, it cannot be instantiated')
            # Specific line for ensembl gff3 non_collateral fix
            info[self.ATTRIBUTE_ID] = info[self.ATTRIBUTE_ID].replace('transcript:', '')
            return Transcript(ID, sequence_name, first, last, length, info)
        elif self.type in (EXON, FIVE_PRIME_UTR, THREE_PRIME_UTR, CDS):
            if parent is None:
                raise ValueError(f'element from {str(self)} has no parent Transcript or Gene, it cannot be '
                                 f'instantiated')
            return TranscriptElement(ID, sequence_name, first, last, length, self.type, info)
        else:
            if ID is None:
                ID = self.GENERIC_ID
            return FunctionalGenomicRegion(ID, sequence_name, first, last, length, self.type, info)

    def __str__(self):
        return f"{self.seqid}{GFF_FIELD_SEPARATOR}{self.source}{GFF_FIELD_SEPARATOR}{self.type}{GFF_FIELD_SEPARATOR}{self.start}{GFF_FIELD_SEPARATOR}" \
               f"{self.end}{GFF_FIELD_SEPARATOR}{self.score}{GFF_FIELD_SEPARATOR}{self.strand}{GFF_FIELD_SEPARATOR}{self.phase}{GFF_FIELD_SEPARATOR}{self.attributes}"
