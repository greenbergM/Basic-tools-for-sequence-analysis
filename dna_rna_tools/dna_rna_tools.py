NUCLEOTIDE_SET = {'A', 'T', 'G', 'C', 'U'}
COMPLEMENTARY_PAIRS_DICT = dict(A=['U', 'T'], a=['u', 't'], G=['C', 'C'], g=['c', 'c'], T=['A', 'A'], t=['a', 'a'],
                                C=['G', 'G'], c=['g', 'g'], U=['A', 'A'], u=['a', 'a'])
DNA_COMPLEMENTARY_MAP = {'A': 'T', 'a': 't', 'G': 'C', 'g': 'c', 'T': 'A', 't': 'a', 'C': 'G', 'c': 'g'}
RNA_COMPLEMENTARY_MAP = {'A': 'U', 'a': 'u', 'G': 'C', 'g': 'c', 'U': 'A', 'u': 'a', 'C': 'G', 'c': 'g'}


def is_dna_rna(dna_rna_list):
    dna_rna_identity = []
    for seq in dna_rna_list:
        nucleotides = set(seq.upper())
        if nucleotides.difference(NUCLEOTIDE_SET):
            raise ValueError("Incorrect input!")
        if 'T' in NUCLEOTIDE_SET.difference(nucleotides) and 'U' in NUCLEOTIDE_SET.difference(nucleotides):
            dna_rna_identity.append('uncertain')
        elif 'U' in NUCLEOTIDE_SET.difference(nucleotides):
            dna_rna_identity.append('DNA')
        elif 'T' in NUCLEOTIDE_SET.difference(nucleotides):
            dna_rna_identity.append('RNA')
        else:
            raise ValueError("There are DNA-RNA hybrids in sequences!")
    return dna_rna_identity


def nacid_identity(dna_rna_identity):
    if len(dna_rna_identity) == 1:
        return dna_rna_identity[0]
    else:
        return dna_rna_identity


def complement_dna(dna_rna_list):
    complement_seq_list = []
    for seq in dna_rna_list:
        comp_dna_seq = str()
        for nucleotide in seq:
            comp_dna_seq += COMPLEMENTARY_PAIRS_DICT[nucleotide][1]
        complement_seq_list += [comp_dna_seq]
    if len(dna_rna_list) == 1:
        return complement_seq_list[0]
    else:
        return complement_seq_list


def complement_rna(dna_rna_list):
    result = []
    for seq in dna_rna_list:
        comp_dna_seq = str()
        for nucleotide in seq:
            comp_dna_seq += COMPLEMENTARY_PAIRS_DICT[nucleotide][0]
        result += [comp_dna_seq]
    if len(dna_rna_list) == 1:
        return result[0]
    else:
        return result


def transcribe(dna_rna_list, dna_rna_identity):
    if 'RNA' in dna_rna_identity:
        raise ValueError("RNA sequences cannot be transcribed! Please, use only DNA sequences.")
    transcribed_seq_list = []
    for seq in dna_rna_list:
        transcribed_seq = seq.replace('T', 'U').replace('t', 'u')
        transcribed_seq_list += [transcribed_seq]

    if len(dna_rna_list) == 1:
        return transcribed_seq_list[0]
    else:
        return transcribed_seq_list


def reverse_transcribe(dna_rna_list, dna_rna_identity):
    if 'DNA' in dna_rna_identity:
        raise ValueError("DNA sequences cannot be reversed transcribed! Please, use only RNA sequences.")
    r_transcribed_seq_list = []
    for seq in dna_rna_list:
        r_transcribed_seq = seq.replace('U', 'T').replace('u', 't')
        r_transcribed_seq_list += [r_transcribed_seq]
    if len(dna_rna_list) == 1:
        return r_transcribed_seq[0]
    else:
        return r_transcribed_seq


def reverse(dna_rna_list):
    result = []
    for seq in dna_rna_list:
        reversed_seq = seq[::-1]
        result += [reversed_seq]
    if len(dna_rna_list) == 1:
        return result[0]
    else:
        return result


def complement(dna_rna_list, dna_rna_identity):
    result = []
    for seq in dna_rna_list:
        seq_identity = dna_rna_identity[dna_rna_list.index(seq)]
        if seq_identity == 'DNA' or seq_identity == 'uncertain':
            complementary_seq = str()
            for nucleotide in seq:
                complementary_seq += COMPLEMENTARY_PAIRS_DICT[nucleotide][1]
            result += [complementary_seq]
        else:
            complementary_seq = str()
            for nucleotide in seq:
                complementary_seq += COMPLEMENTARY_PAIRS_DICT[nucleotide][0]
            result += [complementary_seq]
    if len(dna_rna_list) == 1:
        return result[0]
    else:
        return result


def reverse_complement(dna_rna_list, dna_rna_identity):
    result = []
    for seq in dna_rna_list:
        seq_identity = dna_rna_identity[dna_rna_list.index(seq)]
        if seq_identity == 'DNA' or seq_identity == 'uncertain':
            complementary_seq = str()
            for nucleotide in seq:
                complementary_seq += COMPLEMENTARY_PAIRS_DICT[nucleotide][1]
            reversed_complementary_seq = complementary_seq[::-1]
            result += [reversed_complementary_seq]
        else:
            complementary_seq = str()
            for nucleotide in seq:
                complementary_seq += COMPLEMENTARY_PAIRS_DICT[nucleotide][0]
            reversed_complementary_seq = complementary_seq[::-1]
            result += [reversed_complementary_seq]
    if len(dna_rna_list) == 1:
        return result[0]
    else:
        return result


def gc_content(dna_rna_list):
    gc_content_list = []
    for seq in dna_rna_list:
        gc_amount = 0
        gc_amount = seq.upper().count('G') + seq.upper().count('C')
        gc_content_list += [str(round(gc_amount / len(seq) * 100, 2)) + '%']
    if len(dna_rna_list) == 1:
        return gc_content_list[0]
    else:
        return gc_content_list


TOOL_DICT = {
    'transcribe': transcribe,
    'reverse_transcribe': reverse_transcribe,
    'reverse': reverse,
    'complement': complement,
    'reverse_complement': reverse_complement,
    'nacid_identity': nacid_identity,
    'gc_content': gc_content,
    'complement_rna': complement_rna,
    'complement_dna': complement_dna
}


def run_dna_rna_tools(*args):
    seqs = args[:-1]
    tool = args[-1]
    dna_rna_identity = is_dna_rna(seqs)
    if tool in tool_dict.keys():
        if tool in {'reverse', 'gc_content', 'complement_dna', 'complement_rna'}:
            return tool_dict[tool](seqs)
        elif tool == 'nacid_identity':
            return tool_dict[tool](seqs)
        else:
            return tool_dict[tool](seqs, dna_rna_identity)
    else:
        raise ValueError(f'{tool} is not an available tool!')
