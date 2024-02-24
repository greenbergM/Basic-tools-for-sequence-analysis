from typing import List, Union
import scripts.dna_rna_tools as nucl
import scripts.protein_tools as prot
import os
from Bio import SeqIO
from Bio import SeqUtils
from abc import ABC, abstractmethod


# def run_dna_rna_tools(*args: str, seq_type='DNA') -> Union[list[str], str, list[float], float]:
#     """
#     Launch desired operation with nucleic acid sequences. Pass comma-separated sequences,
#     additional argument (if certain function requires it) and specify function name you want to apply to all sequences.
#     Pass arguments strictly in this order, otherwise it won't be parsed.
#
#     :param args:
#     - nucleic acid sequences for analysis (str)
#     - operation name (str): specify procedure you want to apply (str)
#     :param seq_type: type of desired complement sequence (for complement function) (str)
#
#     :return: the result of procedure in list or str format
#     """
#     seqs = args[:-1]
#     tool = args[-1]
#     seqs_identity = nucl.nucl_acid_identity(seqs)
#     if tool == 'nucl_acid_identity':
#         return seqs_identity
#     elif tool == 'complement':
#         return nucl.complement(seqs, seq_type)
#     elif tool == 'reverse_complement':
#         return nucl.reverse_complement(seqs, seq_type)
#     elif tool == 'transcribe':
#         return nucl.transcribe(seqs, seqs_identity)
#     elif tool == 'reverse':
#         return nucl.reverse(seqs)
#     elif tool == 'gc_content':
#         return nucl.gc_content(seqs)
#     else:
#         raise ValueError(f'There is no {tool} tool available!')
#
#
# def run_protein_analysis(*args: str, site_of_interest=None) -> Union[List[str], str, list[float], float]:
#     """
#     Launch desired operation with proteins sequences. Pass comma-separated sequences,
#     additional argument (if certain function requires it) and specify function name you want to apply to all sequences.
#     Pass arguments strictly in this order, otherwise it won't be parsed.
#
#     :param args:
#     - seq (str): amino acids sequences for analysis in 1-letter or 3-letter code (as many as you wish)
#     - curr_encoding (int): type of encoding of given protein sequences
#     - operation name (str): specify procedure you want to apply
#     :param site_of_interest: one letter encoding of desired site (for find_site function)
#     :return: the result of procedure in list or str format
#     """
#     tool = args[-1]
#     curr_encoding = args[-2]
#     seqs = args[:-2]
#     processed_result = []
#     for seq in seqs:
#         if not prot.is_protein(seq, curr_encoding):
#             raise ValueError('Please, use protein sequences!')
#     seqs = prot.change_encoding(seqs, curr_encoding)
#     if tool == 'get_seq_characteristic':
#         for seq in seqs:
#             processed_result.append(prot.get_seq_characteristic(seq))
#     elif tool == 'find_site':
#         for seq in seqs:
#             processed_result.append(prot.find_site(seq, site_of_interest))
#     elif tool == 'calculate_protein_mass':
#         for seq in seqs:
#             processed_result.append(prot.calculate_protein_mass(seq))
#     elif tool == 'calculate_average_hydrophobicity':
#         for seq in seqs:
#             processed_result.append(prot.calculate_average_hydrophobicity(seq))
#     elif tool == 'calculate_isoelectric_point':
#         for seq in seqs:
#             processed_result.append(prot.calculate_isoelectric_point(seq))
#     elif tool == 'get_mrna':
#         for seq in seqs:
#             processed_result.append(prot.get_mrna(seq))
#     else:
#         raise ValueError(f'{tool} operation is not available!')
#     return processed_result[0] if len(processed_result) == 1 else processed_result


class IncorrectNucleotideError(ValueError):
    pass


class IncorrectAminoacidEncodingError(ValueError):
    pass


class IncorrectAminoacidError(ValueError):
    pass



class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass


class NucleicAcidSequence(BiologicalSequence, ABC):
    def __init__(self, sequence):
        self.sequence = sequence
        self.nucleotide_pairs = {}

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return str(self.sequence)

    def check_alphabet(self):
        nucleotides = self.nucleotide_pairs.keys()

        for nucleotide in set(self.sequence):
            if nucleotide not in nucleotides:
                raise IncorrectNucleotideError(f"Invalid nucleotide found: {nucleotide}")
        return True

    def complement(self):
        self.check_alphabet()

        complement_seq = str()
        for nucleotide in self.sequence:
            complement_seq += self.nucleotide_pairs[nucleotide]
        complement_seq_object = type(self)(complement_seq)
        return complement_seq_object

    def get_gc_content(self):
        self.check_alphabet()

        gc_content = SeqUtils.GC(self.sequence.upper())
        return gc_content


class DNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    def __init__(self, seq):
        super().__init__(seq)
        self.complement_pairs = {'A': 'T',
                                 'a': 't',
                                 'G': 'C',
                                 'g': 'c',
                                 'T': 'A',
                                 't': 'a',
                                 'C': 'G',
                                 'c': 'g'}

    def transcribe(self):
        transcribed_seq = self.sequence.replace('T', 'U').replace('t', 'u')
        return RNASequence(transcribed_seq)


class RNASequence(NucleicAcidSequence, BiologicalSequence, ABC):
    def __init__(self, seq):
        super().__init__(seq)
        self.complement_pairs = {'A': 'U',
                                 'a': 'u',
                                 'G': 'C',
                                 'g': 'c',
                                 'U': 'A',
                                 'u': 'a',
                                 'C': 'G',
                                 'c': 'g'}


class AminoAcidSequence(BiologicalSequence, ABC):
    def __init__(self, sequence, encoding):
        self.sequence = sequence
        self.encoding = encoding
        self.residue_names = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                              'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F',
                              'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}
        self.residue_mass = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147, 'G': 75, 'H': 155,
                             'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 'S': 105, 'T': 119, 'W': 204,
                             'Y': 181, 'V': 117}

    def _reformat_based_on_encoding(self):
        if self.encoding == 1:
            return self.sequence
        elif self.encoding == 3:
            amino_acid_list = [self.sequence[letter:letter + 3] for letter in range(0, len(self.sequence), 3)]
            return amino_acid_list
        else:
            raise IncorrectAminoacidEncodingError(f'{self.encoding}-letter encoding is unavailable. '
                                                  f'Please, use 1 or 3 letter encoding ')

    def _make_one_letter(self):
        if self.encoding == 3:
            one_letter_seq = str()
            for aminoacid in self._reformat_based_on_encoding():
                one_letter_seq += self.residue_names[aminoacid]
            return one_letter_seq
        return self.sequence

    def __len__(self):
        return len(self._reformat_based_on_encoding())

    def __getitem__(self, index):
        return self._reformat_based_on_encoding()[index]

    def __str__(self):
        return str(self.sequence)

    def check_alphabet(self):
        aminoacids = set(self._reformat_based_on_encoding())
        for aminoacid in aminoacids:
            if aminoacid not in self.residue_names.keys() and aminoacid not in self.residue_names.values():
                raise IncorrectAminoacidError(f'{aminoacid} is not a supported aminoacid!')

    def get_molecular_mass(self):
        sequence = self._make_one_letter()
        mass = 0
        for aminoacid in sequence:
            mass += self.residue_mass[aminoacid]
        return mass


def filter_fastq(input_path: str, output_filename=None, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32),
                 quality_threshold=0):
    """
        Filter DNA sequences based on the GC-content, length and sequencing quality (phred33) from FASTQ file to a new
        FASTQ file and stores it in fastaq_filtered_results folder in the same directory.

        :param input_path: path to the sequences in fastq format
        :param output_filename: name for output fastq file (str);
        if not given the output file name will be filtered_*input file name*.fastq
        :param gc_bounds: given threshold for GC-content (tuple/int/float)
        :param length_bounds: given threshold for length (tuple/int)
        :param quality_threshold: given threshold for quality (int)
        """

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    output_dir = os.path.join(os.path.dirname(input_path), 'fastq_filtrator_results')
    os.makedirs(output_dir, exist_ok=True)

    if output_filename:
        output_file_path = os.path.join(output_dir, output_filename)
    else:
        output_file_path = os.path.join(output_dir, f'filtered_{os.path.basename(input_path)}')

    with open(output_file_path, "w") as f:
        for record in SeqIO.parse(input_path, "fastq"):
            gc_test = gc_bounds[0] < SeqUtils.GC(record.seq) < gc_bounds[1]
            len_test = length_bounds[0] < len(record.seq) < length_bounds[1]
            quality_test = sum(record.letter_annotations["phred_quality"]) / len(record.seq) >= quality_threshold
            if gc_test and len_test and quality_test:
                SeqIO.write(record, f, "fastq")
