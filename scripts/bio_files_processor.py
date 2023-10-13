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
