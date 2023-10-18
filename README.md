# SeqMaster

This repository contains an open-source library which make work with protein, DNA and RNA. It consists out of two main scripts. 

SeqMaster.py is capable to process multiple proteins and peptides sequences, calculate physical features, find specific sites. For nucleic acids it can transcribe, find complement sequences, reverse, calculate GC-content. Also SeqMaster supports FASTQ-format and can be used to filter DNA sequences based on their length, GC-content and sequence quality. 

bio_files_proccessor.py is a script that allow user to work with biological files in popular formats. It can convert FASTA files with multiline sequence records to oneline, find neighbour CDSs for given genes from GBK file, change starting position in FASTA file (usefull for prokaryotic genomes) and find descriptions for best BLAST results. All results are saved in new folders in the same directory with input files. 


## Installation

To use this toolbox one need to clone repository

```shell
git clone https://github.com/greenbergM/SeqMaster
cd SeqMaster
```

### System requirements:

Key packages and programs:
- [Python](https://www.python.org/downloads/) (version >= 3.9)

## SeqMaster.py

### Usage

```python
# import main function
from SeqMaster import run_dna_rna_tools, run_protein analysis, run_filter_fastaq

```

### Sequence tools

This section contains description of functions you can find in our library.

#### DNA and RNA tools

This function allows you to perform various operations on nucleic acid sequences. Pass comma-separated sequences as the first arguments, followed by an additional argument (if required by the specific function), and specify the function name you want to apply to all sequences. Make sure to provide arguments strictly in this order.

##### Usage

`run_dna_rna_tools(*args, seq_type = 'DNA')`

##### Arguments

**args**:
- DNA or RNA sequnces in `str` format
- operation name in str format

    **Note**
	- There can be any number of sequnces
	- Sequences can contain upper case and lower case letters
	- Operation always must be the last argument!

**seq_type**:
- type of desired complementary sequence
- by default has value 'DNA'
- only supported values are 'DNA' and 'RNA'

##### Operations

- `nucl_acid_identity` - shows the type of given sequences; if sequence cannot be identified as `'DNA'` or `'RNA'` it will be given name `'uncertain'`

  *Returns* 
	`result: dict` 

- `reverse` - gives the reversed copy of a sequence

  *Returns*
	`result: list[str] or str (if one sequence given)` 

- `complement` - gives complement `DNA` or `RNA` sequence depending on the `seq_type` argument. 

  *Returns*
	`result: list[str] or str (if one sequence given)` 

- `reverse_complement` - gives reversed complement `DNA` or `RNA` sequence depending on the `seq_type` argument. 

  *Returns*
	`result: list[str] or str (if one sequence given)` 

- `transcribe` - transcribes given sequences; DNA sequences are transcribed to RNA, RNA sequences are reversed transcribed to DNA.

  *Returns*
	`result: list[str] or str (if one sequence given)` 


- `gc_content` - gives the GC-content in percent for given sequences

  *Returns*
	`result: list[float] or float (if one sequence given)`

##### Example

```python
seq1= 'AGTC'
seq2 = 'AUGGG'

run_dna_rna_tools(seq1, 'nucl_acid_identity')
run_dna_rna_tools(seq1, 'complement', seq_type = 'RNA')
run_dna_rna_tools(seq1, seq2, 'transcribe')

```

#### Protein tools

This function allows you to perform various operations on protein sequences. Pass comma-separated sequences as the first arguments, followed by an additional argument (if required by the specific function), and specify the function name you want to apply to all sequences. Make sure to provide arguments strictly in this order.

##### Usage

`run_protein_analysis(*args, site_of_interest = None)`

##### Arguments

**args**:
- protein sequnces in `str` format
- sequence encoding type in `int` format (`1` - one-letter encoding or `3` - three-letter encoding)
- operation on `str` format

    **Note**
	- There can be any number of sequnces
	- Sequences can contain upper case and lower case letters
	- Sequences can be given in one-letter format or three-letter format (it must be clarified by sequence encoding type parameter); in three-letter format they can contain spaces between amino acids. 
	- Arguments always must be in the strict order: sequences, sequence encoding type, operation!

**site_of_interest** 
 - sequence of site to find in given protein sequences via `find_site` operation
 - must be in one-letter encoding
 - by default set to be `None`

 ##### Operations

- `get_seq_characteristic` - counts entry of each residue type in given sequences.

  *Returns*
	`result: list[dict] or dict (if one sequence given)` 

- `find_site` - find if sequence contains certain site (given by `site_of_interest` argument) and get positions of this site

  *Returns*
	`result: list[str] or str (if one sequence given)` 

- `calculate_protein_mass` - calculates mass of residues in given sequences in Da

  *Returns*
	`result: list[float] or float (if one sequence given)` 

- `calculate_average_hydrophobicity` - calculates hydrophobicity index for each given sequence as sum of index for each residue in sequence divided by its length

  *Returns*
	`result: list[float] or float (if one sequence given)` 

- `calculate_isoelectric_point` - calculates soelectrinc point for each given sequence as sum of known pI for residues in sequence

  *Returns*
	`result: list[float] or float (if one sequence given)` 

- `get_mrna` - gives encoding mRNA nucleotides in degenerate form for each given sequence

  *Returns*
	`result: list[str] or str (if one sequence given)` 

##### Example

```python
seq1= 'His Val Leu Ile'
seq2 = 'HVGGLLA'

run_protein_analysis(seq1, 3,  'get_mrna')
run_protein_analysis(seq2, 1, 'find_site', site_of_interest = 'GGL')
run_protein_analysis(seq1, seq1, 3, 'calculate_average_hydrophobicity')

```


#### FASTAQ filter

This func filters DNA sequences based on GC-content, length, and sequencing quality (phred33) from a FASTQ file and saves it in the `fastq_filtrator_results` directory within the same location as the input file.


##### Usage

`run_filter_fastaq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_threshold=0)`


##### Arguments

**input_path** 
- Path to the FASTQ file in `str` format.

**output_filename** 
- Name for the output FASTQ file. If not provided, the output file will be named as `filtered_[input file name].fastq`.

**gc_bounds**
- lower and upper boundaries of GC-content in percent orginised in `tuple`
- by default is set to `(0, 100)`
- if one number given it will be considered the upper boundary with the lover boundary = 0

**length_bounds**
- lower and upper boundaries of sequence length orginised in `tuple`
- by default is set to `(0, 2 ** 32)`
- if one number given it will be considered the upper boundary with the lover boundary = 0

**quality_threshold**
- threshold average Q-score (`int`) - lower boundary 
- by default is set to 0

##### Example

```python
run_filter_fastaq('sequences.fastq', output_filename='filtered_sequences.fastq', gc_bounds=(30.5, 70), length_bounds=(50, 200), quality_threshold=30)

# This will create a file named filtered_sequences.fastq in the fastq_filtrator_results directory with filtered sequences.
```


## bio_files_processor.py


### Usage

```python
# import main function
from bio_files_processor import convert_multiline_fasta_to_oneline, select_genes_from_gbk_to_fasta, change_fasta_start_pos, parse_blast_output

```

#### Convert Multiline FASTA to Oneline

This script converts a FASTA file with multiline sequences to a FASTA file with oneline sequences. The resulting file is stored in the `oneline_fasta` directory within the same location as the input file.

##### Usage

`convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta=None)`

##### Arguments

- **input_fasta**  - Path to the multiline FASTA file in `str` format.

- **output_fasta**  - Name of the oneline FASTA file in `str` format. If not provided, the output file will be named as `oneline_[input file name].fasta`.

##### Example

```python
convert_multiline_fasta_to_oneline('input.fasta', 'output.fasta')

#This will create an oneline FASTA file named output.fasta based on the contents of input.fasta
```

#### Select CDSs from GBK to FASTA

This script extracts neighbor CDSs to specified genes from a GBK file and generates a FASTA file containing the sequences. The resulting file is stored in the `fasta_selected_from_gbk` directory within the same location as the input file.

##### Usage

`select_genes_from_gbk_to_fasta(*genes: str, input_gbk: str, n_before: int, n_after: int, output_fasta=None)`

##### Arguments

- **genes** - Genes of interest to be used for neighbor CDSs search in `str` format.

- **input_gbk** - Path to the GBK file in `str` format.

- **n_before** - Number of neighbor CDSs before the gene of interest in `int` format.

- **n_after** - Number of neighbor CDSs after the gene of interest in `int` format.

- **output_fasta** - Name of the output FASTA file. If not provided, the output file will be named as `CDS_selected_from_[input file name].fasta`.

##### Example

```python
select_genes_from_gbk_to_fasta('gene1', 'gene2', input_gbk='input.gbk', n_before=2, n_after=2, output_fasta='output.fasta')

# This will create a FASTA file named output.fasta containing sequences of neighbor CDSs to 'gene1' and 'gene2' from the input.gbk file.
```

#### Shift start position in FASTA Sequence

This script shifts the starting position of a DNA sequence in a FASTA file and saves it in the `shifted_fasta` directory within the same location as the input file.

**Note**
Only works with FASTA files containing one sequence!

##### Usage

`change_fasta_start_pos(input_fasta: str, shift: int, output_fasta=None)`

##### Arguments

- **input_fasta** - Path to the input FASTA file in `str` format.

- **shift** - The number of positions to shift the sequence in `int` format.

- **output_fasta** - Name of the output FASTA file in `str` format. If not provided, the output file will be named as `shifted_[input file name].fasta`.

##### Example

```python
change_fasta_start_pos('input.fasta', 2, output_fasta='output.fasta')

# This will create a FASTA file named output.fasta where the DNA sequence starting position is shifted by 2 positions.

# if input file contained `AAACGT` sequence, the output.fasta will contain 'ACGTAA'

change_fasta_start_pos('input.fasta', -2, output_fasta='output.fasta')

# if input file contained `AAACGT` sequence, the output.fasta will contain 'GTAAAC'

```

#### Parse blast output

This script extracts the descriptions (gene names) of the best BLAST results from a blast results file and saves it in the `best_blast_results` directory within the same location as the input file.

##### Usage

`parse_blast_output(input_file: str, output_file=None, extension='.txt')`

##### Arguments

- **input_file** - Path to the blast results file in `str` format.

- **output_file** - Name for the output file in `str` format. If not provided, the output file will be named as `best_blast_[input file name].[extension]`.

- **extension** - Extension for the output file in `str` format. If not changed, it will be `.txt`.

##### Example

```python
parse_blast_output('blast_results.txt', output_file='output_results', extension='.txt')

# This will create a file named output_results.txt with the descriptions (=names) of the best BLAST results.
```



