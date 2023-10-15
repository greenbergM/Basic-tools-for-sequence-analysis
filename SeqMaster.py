from typing import List, Union
import os
import scripts.dna_rna_tools as nucl
import scripts.protein_tools as prot
import scripts.fastaq_filter as fasq
import scripts.bio_files_processor as bfp


def run_dna_rna_tools(*args: str, seq_type='DNA') -> Union[list[str], str, list[float], float]:
    """
    Launch desired operation with nucleic acid sequences. Pass comma-separated sequences,
    additional argument (if certain function requires it) and specify function name you want to apply to all sequences.
    Pass arguments strictly in this order, otherwise it won't be parsed.

    :param args:
    - nucleic acid sequences for analysis (str)
    - operation name (str): specify procedure you want to apply (str)
    :param seq_type: type of desired complement sequence (for complement function) (str)

    :return: the result of procedure in list or str format
    """
    seqs = args[:-1]
    tool = args[-1]
    seqs_identity = nucl.nucl_acid_identity(seqs)
    if tool == 'nucl_acid_identity':
        return seqs_identity
    elif tool == 'complement':
        return nucl.complement(seqs, seq_type)
    elif tool == 'reverse_complement':
        return nucl.reverse_complement(seqs, seq_type)
    elif tool == 'transcribe':
        return nucl.transcribe(seqs, seqs_identity)
    elif tool == 'reverse':
        return nucl.reverse(seqs)
    elif tool == 'gc_content':
        return nucl.gc_content(seqs)
    else:
        raise ValueError(f'There is no {tool} tool available!')


def run_protein_analysis(*args: str, site_of_interest=None) -> Union[List[str], str, list[float], float]:
    """
    Launch desired operation with proteins sequences. Pass comma-separated sequences,
    additional argument (if certain function requires it) and specify function name you want to apply to all sequences.
    Pass arguments strictly in this order, otherwise it won't be parsed.

    :param args:
    - seq (str): amino acids sequences for analysis in 1-letter or 3-letter code (as many as you wish)
    - curr_encoding (int): type of encoding of given protein sequences
    - operation name (str): specify procedure you want to apply
    :param site_of_interest: one letter encoding of desired site (for find_site function)
    :return: the result of procedure in list or str format
    """
    tool = args[-1]
    curr_encoding = args[-2]
    seqs = args[:-2]
    processed_result = []
    for seq in seqs:
        if not prot.is_protein(seq, curr_encoding):
            raise ValueError('Please, use protein sequences!')
    seqs = prot.change_encoding(seqs, curr_encoding)
    if tool == 'get_seq_characteristic':
        for seq in seqs:
            processed_result.append(prot.get_seq_characteristic(seq))
    elif tool == 'find_site':
        for seq in seqs:
            processed_result.append(prot.find_site(seq, site_of_interest))
    elif tool == 'calculate_protein_mass':
        for seq in seqs:
            processed_result.append(prot.calculate_protein_mass(seq))
    elif tool == 'calculate_average_hydrophobicity':
        for seq in seqs:
            processed_result.append(prot.calculate_average_hydrophobicity(seq))
    elif tool == 'calculate_isoelectric_point':
        for seq in seqs:
            processed_result.append(prot.calculate_isoelectric_point(seq))
    elif tool == 'get_mrna':
        for seq in seqs:
            processed_result.append(prot.get_mrna(seq))
    else:
        raise ValueError(f'{tool} operation is not available!')
    return processed_result[0] if len(processed_result) == 1 else processed_result


def run_filter_fastaq(input_path: str, output_filename=None, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0) -> dict[str: str]:
    """
    Filter DNA sequences based on the GC-content, length and sequencing quality (phred33).
    :param input_path: path to the sequences in fastq format
    :param output_filename: name for output fastq file (str)
    :param gc_bounds: given threshold for GC-content (tuple/int)
    :param length_bounds: given threshold for length (tuple/int)
    :param quality_threshold: given threshold for quality (int)
    :return: given dict only with sequences that fit all the criteria.
    """
    seqs = fasq.get_dict(input_path)
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_seqs = dict()
    for seq_name in seqs.keys():
        gc_result = fasq.gc_test(seqs[seq_name][0], gc_bounds)
        length_result = fasq.length_test(seqs[seq_name][0], length_bounds)
        quality_result = fasq.quality_test(seqs[seq_name][2], quality_threshold)
        if fasq.judge_seq(gc_result, length_result, quality_result):
            filtered_seqs[seq_name] = seqs[seq_name]

    fasq.get_file(filtered_seqs, output_filename, input_path)
    print(f'Filtering completed; {len(filtered_seqs.keys())} sequences were selected from {len(seqs.keys())} given.')


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None):
    """
    Creates fasta file with oneline sequences based on given fasta file with multiline sequences in the same directory.
    :param input_fasta: path to the multiline fasta file (str)
    :param output_fasta: name of oneline fasta file (str)
    """
    if output_fasta is None:
        output_fasta = f'oneline_{os.path.basename(os.path.realpath(input_fasta))}'
    if not output_fasta.endswith('.fasta'):
        output_fasta += '.fasta'

    location = os.path.dirname(os.path.realpath(input_fasta))
    seq = ''
    with open(input_fasta, mode='r') as fr, open(os.path.join(location, output_fasta), mode='w') as fw:
        for line in fr:
            if line.startswith('>'):
                if seq:
                    fw.write(seq + '\n')
                    seq = ''
                fw.write(line)
            else:
                seq += line.strip('\n')
        fw.write(seq)


def select_genes_from_gbk_to_fasta(*genes: str, input_gbk: str, n_before: int, n_after: int, output_fasta=None):
    """
    Creates fasta file with neighbor CDSs to given genes from GBK file and stores it
    in fasta_selected_from_gbk directory.
    :param genes: genes of interest that are used for neighbor CDSs search (str)
    :param input_gbk: path to GBK file (str)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output fasta file
    """
    cds_list = bfp.get_cds_list(input_gbk)
    gene_cds_dict = bfp.get_gene_to_cds_dict(input_gbk)
    translation_dict = bfp.get_cds_translation_dict(input_gbk)
    cds_of_interest = []

    for gene in genes:
        gene_cds = gene_cds_dict[gene]
        gene_position = cds_list.index(gene_cds)
        cds_before = cds_list[gene_position - n_before:gene_position]
        cds_after = cds_list[gene_position + 1:gene_position + n_after + 1]
        cds_of_interest = cds_of_interest + cds_before + cds_after

    fasq.check_dir('fasta_selected_from_gbk')

    if output_fasta is None:
        output_fasta = 'CDS_selected_from_gbk.fasta'
    if not output_fasta.endswith('.fasta'):
        output_fasta = output_fasta + '.fasta'

    with open(os.path.join('fasta_selected_from_gbk', output_fasta), mode='w') as fasta:
        for cds in cds_of_interest:
            fasta.write('>' + cds + ' gene: ' + translation_dict[cds][0] + '\n')
            fasta.write((translation_dict[cds][1].replace('"', '') + '\n'))
