from typing import List, Union

NUCLEOTIDE_SET = {'A', 'T', 'G', 'C', 'U'}
DNA_COMPLEMENTARY_MAP = {'A': 'T', 'a': 't', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'C': 'G', 'c': 'g', 'U': 'A',
                         'u': 'a'}
RNA_COMPLEMENTARY_MAP = {'A': 'U', 'a': 'u', 'G': 'C', 'g': 'c', 'U': 'A', 'u': 'a', 'C': 'G', 'c': 'g', 'T': 'A',
                         't': 'a'}


def nucl_acid_identity(seqs: tuple[str]) -> dict:
    """
    Check if the sequences are DNA, RNA or uncertain nucleic acid
    by identify invalid seq elements, which are not presented in dicts above.
    :param seqs: sequences of nucleic acids (tuple[str])
    :return: dictionary, where the sequence is key and nucleic acid type is its value (dict)
    """
    seqs_identity = dict()
    for seq in seqs:
        nucleotides = set(seq.upper())
        if nucleotides.difference(NUCLEOTIDE_SET):
            raise ValueError("Incorrect input!")
        if 'T' in NUCLEOTIDE_SET.difference(nucleotides) and 'U' in NUCLEOTIDE_SET.difference(nucleotides):
            seqs_identity[seq] = 'uncertain'
        elif 'U' in NUCLEOTIDE_SET.difference(nucleotides):
            seqs_identity[seq] = 'DNA'
        elif 'T' in NUCLEOTIDE_SET.difference(nucleotides):
            seqs_identity[seq] = 'RNA'
        else:
            raise ValueError("There are DNA-RNA hybrids in sequences!")
    return seqs_identity


def complement(seqs: tuple[str], seq_type: str) -> Union[list[str], str]:
    """
    Get complement DNA or RNA sequence
    :param seqs: sequences of nucleic acids (tuple[str])
    :param seq_type: type of desired complement sequence (str)
    :return: list of complement sequences ([list[str]); if given one sequence it will be returned in string (str)
    """
    complement_seqs = []
    for seq in seqs:
        compl_seq = str()
        for nucleotide in seq:
            if seq_type == 'DNA':
                compl_seq += DNA_COMPLEMENTARY_MAP[nucleotide]
            elif seq_type == 'RNA':
                compl_seq += RNA_COMPLEMENTARY_MAP[nucleotide]
            else:
                raise ValueError(f'{seq_type} sequence type is not available. Please, use DNA or RNA')
        complement_seqs.append(compl_seq)
    return complement_seqs[0] if len(seqs) == 1 else complement_seqs


def reverse_complement(seqs: tuple[str], seq_type: str) -> Union[list[str], str]:
    """
    Get reversed complement DNA or RNA sequence
    :param seqs: sequences of nucleic acids (tuple[str])
    :param seq_type: type of desired complement sequence (str)
    :return: list of reversed complement sequences ([list[str]); if given one sequence it will be returned in string (str)
    """
    complement_seqs = complement(seqs, seq_type)
    if isinstance(complement_seqs, str):
        return complement_seqs[::-1]
    else:
        return reverse(complement_seqs)


def transcribe(seqs: tuple[str], seqs_identity: dict) -> Union[list[str], str]:
    """
    Get transcribed sequences of DNA and reversed transcribed sequences of RNA
    :param seqs: sequences of nucleic acids (tuple[str])
    :param seqs_identity:
    :return: list of complement sequences ([list[str]); if given one sequence it will be returned in string (str)
    """
    transcribed_seq_list = []
    for seq in seqs:
        if seqs_identity[seq] == 'DNA':
            transcribed_seq_list.append(seq.replace('T', 'U').replace('t', 'u'))
        elif seqs_identity[seq] == 'RNA':
            transcribed_seq_list.append(seq.replace('U', 'T').replace('u', 't'))
        else:
            transcribed_seq_list.append(seq)
    return transcribed_seq_list[0] if len(seqs) == 1 else transcribed_seq_list


def reverse(seqs: tuple[str]) -> Union[list[str], str]:
    reversed_seq_list = []
    for seq in seqs:
        reversed_seq = seq[::-1]
        reversed_seq_list += [reversed_seq]
    return reversed_seq_list[0] if len(seq_list) == 1 else reversed_seq_list


def gc_content(seqs: tuple[str]) -> Union[list[float], float]:
    gc_content_list = []
    for seq in seqs:
        gc_amount = seq.upper().count('G') + seq.upper().count('C')
        gc_content_list += [round(gc_amount / len(seq) * 100, 2)]
    if len(seqs) == 1:
        return gc_content_list[0]
    else:
        return gc_content_list


def run_dna_rna_tools(*args: str, seq_type='DNA') -> Union[list[str], str, list[float], float]:
    seqs = args[:-1]
    tool = args[-1]
    seqs_identity = nucl_acid_identity(seqs)
    if tool == 'nucl_acid_identity':
        return seqs_identity
    elif tool == 'complement':
        return complement(seqs, seq_type)
    elif tool == 'reverse_complement':
        return reverse_complement(seqs, seq_type)
    elif tool == 'transcribe':
        return transcribe(seqs, seqs_identity)
    elif tool == 'reverse':
        return reverse(seqs)
    elif tool == 'gc_content':
        return gc_content(seqs)
    else:
        raise ValueError(f'There is no {tool} tool available!')
