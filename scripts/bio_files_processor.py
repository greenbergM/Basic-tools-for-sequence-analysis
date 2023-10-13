def convert_multiline_fasta_to_oneline(input_fasta, output_fasta):
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

    with open(output_fasta, mode='w') as oneline_fasta:
        for name, seq in oneline_fasta_dict.items():
            oneline_fasta.write(name)
            oneline_fasta.write(seq + '\n')



