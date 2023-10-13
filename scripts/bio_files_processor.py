import os


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):
    oneline_fasta_dict = {}
    seq = ''
    with open(input_fasta, mode='r') as multiline_fasta:
        for line in multiline_fasta:
            if line.startswith('>'):
                name = line
            else:
                while not line.startswith('>'):
                    seq += (line.strip('\n'))
                oneline_fasta_dict[name] = seq
                seq = ''
    if output_fasta is None:
        output_fasta = 'oneline_' + os.path.basename(os.path.realpath(input_fasta))
    else:
        output_fasta += '.fasta'
    location = os.path.dirname(os.path.realpath(input_fasta))
    with open(os.path.join(location, output_fasta), mode='w') as oneline_fasta:
        for name, seq in oneline_fasta_dict.items():
            oneline_fasta.write(name)
            oneline_fasta.write(seq + '\n')