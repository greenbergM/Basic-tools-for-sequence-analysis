import scripts.bio_files_processor_scripts as bfp
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


def select_genes_from_gbk_to_fasta(*genes: str, input_gbk: str, n_before: int, n_after: int, output_fasta=None):
    """
    Creates fasta file with neighbor CDSs to given genes from GBK file and stores it
    in fasta_selected_from_gbk directory.
    :param genes: genes of interest that are used for neighbor CDSs search (str)
    :param input_gbk: path to GBK file (str)
    :param n_before: number of neighbor CDSs before gene of interest (int)
    :param n_after: number of neighbor CDSs after gene of interest (int)
    :param output_fasta: name of the output fasta file (str)
    """
    cds_list = bfp.get_cds_list(input_gbk)
    gene_cds_dict = bfp.get_gene_to_cds_dict(input_gbk)
    translation_dict = bfp.get_cds_translation_dict(input_gbk)
    cds_of_interest = []

    for gene in genes:
        gene_cds = gene_cds_dict[gene]
        gene_position = cds_list.index(gene_cds)
        if gene_position - n_before < 0:
            raise ValueError(f'CDCs before {gene} are out of bounds! Use other value for n_before.')
        cds_before = cds_list[gene_position - n_before:gene_position]
        cds_after = cds_list[gene_position + 1:gene_position + n_after + 1]
        cds_of_interest = cds_of_interest + cds_before + cds_after

    os.makedirs('fasta_selected_from_gbk', exist_ok=True)

    if output_fasta is None:
        output_fasta = 'CDS_selected_from_gbk.fasta'
    if not output_fasta.endswith('.fasta'):
        output_fasta = output_fasta + '.fasta'

    with open(os.path.join('fasta_selected_from_gbk', output_fasta), mode='w') as fasta:
        for cds in cds_of_interest:
            fasta.write('>' + cds + ' gene: ' + translation_dict[cds][0] + '\n')
            fasta.write((translation_dict[cds][1].replace('"', '') + '\n'))

    print('FASTA file for neighbour CDSs of given genes is created ')
