import os


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


def get_cds_list(input_gbk: str) -> list:
    """
    Creates list of all CDS from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: list of all CDSs (list)
    """
    cds_list = []
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith("     CDS"):
                cds_list.append(line.replace(' ', '').strip('\n'))
    return cds_list


def get_gene_to_cds_dict(input_gbk: str) -> dict:
    """
    Creates dict where genes are keys and CDSs are values from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: dict where genes are keys and CDSs are values (dict)
    """
    gene_dict = {}

    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                current_cds = line.replace(' ', '').strip('\n')
            if '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
                gene_dict[current_gene] = current_cds
    return gene_dict


def get_cds_translation_dict(input_gbk: str) -> dict:
    """
    Creates dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file.
    :param input_gbk: path to GBK file (str)
    :return: dict where CDSs are keys and gene names (if they are presented) and translation sequences
    are values from GBK file (dict).
    """
    cds_dict = {}
    current_gene = ''
    current_seq = ''
    translation = False
    with open(input_gbk, mode='r') as gbk:
        for line in gbk:
            if line.startswith('     CDS'):
                translation = False
                if current_seq:
                    cds_dict[current_cds] = (current_gene, current_seq)
                    current_gene = ''
                    current_seq = ''
                current_cds = line.replace(' ', '').strip('\n')
            elif '/gene=' in line:
                current_gene = line.replace(' ', '').strip('\n').strip('/gene=').strip('"')
            elif '/translation' in line:
                current_seq = line.replace(' ', '').strip('\n').strip('/translation=')
                translation = True
            elif 'ORIGIN' in line:
                translation = False
            elif translation:
                current_seq += line.replace(' ', '').strip('\n')
    cds_dict[current_cds] = (current_gene, current_seq)
    return cds_dict
