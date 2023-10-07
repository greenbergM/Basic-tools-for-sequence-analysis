NUCLEOTIDE_SET = {'A', 'T', 'G', 'C', 'U'}
DNA_COMPLEMENTARY_MAP = {'A': 'T', 'a': 't', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'C': 'G', 'c': 'g', 'U': 'A',
                         'u': 'a'}
RNA_COMPLEMENTARY_MAP = {'A': 'U', 'a': 'u', 'G': 'C', 'g': 'c', 'U': 'A', 'u': 'a', 'C': 'G', 'c': 'g', 'T': 'A',
                         't': 'a'}


def is_dna_rna(seqs):
    dna_rna_identity = dict()
    for seq in seqs:
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


def complement(seqs, seq_type):
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


def reverse_complement(seqs, seq_type):
    complement_seqs = complement(seqs, seq_type)
    if isinstance(complement_seqs, str):
        return complement_seqs[::-1]
    else:
        return reverse(complement_seqs)


def transcribe(seqs):
    identity_dict = is_dna_rna(seqs)
    transcribed_seq_list = []
    for seq in seqs:
        if identity_dict[seq] == 'DNA':
            transcribed_seq_list.append(seq.replace('T', 'U').replace('t', 'u'))
        elif identity_dict[seq] == 'RNA':
            transcribed_seq_list.append(seq.replace('U', 'T').replace('u', 't'))
        else:
            transcribed_seq_list.append(seq)
    return transcribed_seq_list[0] if len(seq_list) == 1 else transcribed_seq_list


def reverse(seqs):
    reversed_seq_list = []
    for seq in seqs:
        reversed_seq = seq[::-1]
        reversed_seq_list += [reversed_seq]
    return reversed_seq_list[0] if len(seq_list) == 1 else reversed_seq_list


def gc_content(seqs):
    gc_content_list = []
    for seq in seqs:
        gc_amount = seq.upper().count('G') + seq.upper().count('C')
        gc_content_list += [round(gc_amount / len(seq) * 100, 2)]
    if len(seqs) == 1:
        return gc_content_list[0]
    else:
        return gc_content_list


def run_dna_rna_tools(*args, seq_type='DNA'):
    seqs = args[:-1]
    tool = args[-1]
    if tool == 'is_dna_rna':
        return is_dna_rna(seqs)
    elif tool == 'complement':
        return complement(seqs, seq_type)
    elif tool == 'reverse_complement':
        return reverse_complement(seqs, seq_type)
    elif tool == 'transcribe':
        return transcribe(seqs)
    elif tool == 'reverse':
        return reverse(seqs)
    elif tool == 'gc_content':
        return gc_content(seqs)
    else:
        raise ValueError(f'There is no {tool} tool available!')