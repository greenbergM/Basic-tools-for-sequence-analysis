NUCLEOTIDE_SET = {'A', 'T', 'G', 'C', 'U'}
DNA_COMPLEMENTARY_MAP = {'A': 'T', 'a': 't', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'C': 'G', 'c': 'g', 'U': 'A',
                         'u': 'a'}
RNA_COMPLEMENTARY_MAP = {'A': 'U', 'a': 'u', 'G': 'C', 'g': 'c', 'U': 'A', 'u': 'a', 'C': 'G', 'c': 'g', 'T': 'A',
                         't': 'a'}


def is_dna_rna(seq_list):
    dna_rna_identity = dict()
    for seq in seq_list:
        nucleotides = set(seq.upper())
        if nucleotides.difference(NUCLEOTIDE_SET):
            raise ValueError("Incorrect input!")
        if 'T' in NUCLEOTIDE_SET.difference(nucleotides) and 'U' in NUCLEOTIDE_SET.difference(nucleotides):
            dna_rna_identity[seq] = 'uncertain'
        elif 'U' in NUCLEOTIDE_SET.difference(nucleotides):
            dna_rna_identity[seq] = 'DNA'
        elif 'T' in NUCLEOTIDE_SET.difference(nucleotides):
            dna_rna_identity[seq] = 'RNA'
        else:
            raise ValueError("There are DNA-RNA hybrids in sequences!")
    return dna_rna_identity


def complement(seq_list, seq_type):
    complement_seq_list = []
    for seq in seq_list:
        comp_seq = str()
        for nucleotide in seq:
            if seq_type == 'DNA':
                comp_seq += DNA_COMPLEMENTARY_MAP[nucleotide]
            elif seq_type == 'RNA':
                comp_seq += RNA_COMPLEMENTARY_MAP[nucleotide]
            else:
                raise ValueError(f'{seq_type} sequence type is not available. Please, use DNA or RNA')
        complement_seq_list.append(comp_seq)
    return complement_seq_list[0] if len(seq_list) == 1 else complement_seq_list


def reverse_complement(seq_list, seq_type):
    complement_seqs = complement(seq_list, seq_type)
    if isinstance(complement_seqs, str):
        return complement_seqs[::-1]
    else:
        return reverse(complement_seq_list)


def transcribe(seq_list):
    identity_dict = is_dna_rna(seq_list)
    transcribed_seq_list = []
    for seq in seq_list:
        if identity_dict[seq] == 'DNA':
            transcribed_seq_list.append(seq.replace('T', 'U').replace('t', 'u'))
        elif identity_dict[seq] == 'RNA':
            transcribed_seq_list.append(seq.replace('U', 'T').replace('u', 't'))
        else:
            transcribed_seq_list.append(seq)
    return transcribed_seq_list[0] if len(seq_list) == 1 else transcribed_seq_list


def reverse(seq_list):
    reversed_seq_list = []
    for seq in seq_list:
        reversed_seq = seq[::-1]
        reversed_seq_list += [reversed_seq]
    return reversed_seq_list[0] if len(seq_list) == 1 else reversed_seq_list


def gc_content(seq_list):
    gc_content_list = []
    for seq in seq_list:
        gc_amount = seq.upper().count('G') + seq.upper().count('C')
        gc_content_list += [round(gc_amount / len(seq) * 100, 2)]
    if len(seq_list) == 1:
        return gc_content_list[0]
    else:
        return gc_content_list


def run_dna_rna_tools(*args, seq_type='DNA'):
    seqs_n_tool_list = []
    for arg in args:
        seqs_n_tool_list.append(arg)
    seq_list = seqs_n_tool_list[:-1]
    tool = seqs_n_tool_list[-1]
    if tool == 'is_dna_rna':
        return is_dna_rna(seq_list)
    elif tool == 'complement':
        return complement(seq_list, seq_type)
    elif tool == 'reverse_complement':
        return reverse_complement(seq_list, seq_type)
    elif tool == 'transcribe':
        return transcribe(seq_list)
    elif tool == 'reverse':
        return reverse(seq_list)
    elif tool == 'gc_content':
        return gc_content(seq_list)
    else:
        raise ValueError(f'There is no {tool} tool available!')

