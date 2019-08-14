import argparse
import os
import shutil

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(
    description='Convert codons in genetic sequence to have a lower G/C count, or to contain more frequent codons')
parser.add_argument('--codon_table', type=str, nargs='?', default='data/codons.txt',
                    help='Space delimited file that maps codons to their respected amino acids')
parser.add_argument('--freq', type=str, nargs='?', default='data/humanCodon_freqPerThousand.csv',
                    help='CSV file with no headers that defines the frequency of given codons in the human genome'),
parser.add_argument('--seq', type=str, nargs='?', default='data/hngly1-mrna.fasta',
                    help='Sequence to optimize G/C content and codon frequency according to human genome'),
parser.add_argument('--output', type=str, nargs='?', default=None,
                    help='Output of the resulting optimization in fasta format')
parser.add_argument('--weight', type=float, nargs='?', default=0.8,
                    help='How aggressively to decrease the count of G/C nucleotides in the codons')

args = parser.parse_args()

# We assume that the csv has no header and is separated by spaces
codons_to_acid = pd.read_csv(args.codon_table, names=['codon', 'amino_acid'], sep=' ', header=None)
# Use only first and third columns which are codons and frequency count in human genome
codon_freq = pd.read_csv(args.freq, usecols=[0, 2])
ngly1_raw = SeqIO.read(args.seq, 'fasta')
# Get the amino acids from the raw nucleotides so we can reverse search for better codons
ngly1_acids = ngly1_raw.translate()
# Combine the codons to acid table with the frequency count in the human genome for each codon
# This allows us to cross-reference counts, codons, and the amino acids they code for
codon_to_acid_with_freq_and_at_count = pd.merge(codons_to_acid, codon_freq)
# Add the A and T count to each codon, our target is to maximize them instead of having G and C's
codon_to_acid_with_freq_and_at_count['at_count'] = codon_to_acid_with_freq_and_at_count.codon.apply(
    lambda codon: codon.count('A') + codon.count('T'))

better_codons = ""


def get_scores(df, label):
    """
    Annotate a given column in a dataframe with a score.
    1.0 will be the highest value, while 0.0 will be the lowest value. It is essentially a re-map
    :param df: The dataframe to operate on
    :param label: The name of the column
    :return: None, this function modifies the given dataframe.
    """
    max_number = df[label].max()
    min_number = df[label].min()
    range_number = max_number - min_number
    df[f'{label}_score'] = 1.0 if range_number == 0 else df[label].apply(
        lambda val: (val - min_number) / range_number)


# Iterate over each acid, find the best codon sequence, and add it to our new gene builder
for acid in ngly1_acids:
    if acid == '*':
        # biopython represents stop as asterisk while in our codons table it is represented as "Stop"
        acid = "Stop"
    # Find codons that can represent this amino acid. We must drop NaN because by default pandas leaves them in.
    potential_codons = codon_to_acid_with_freq_and_at_count.where(
        codon_to_acid_with_freq_and_at_count['amino_acid'] == acid).dropna()
    get_scores(potential_codons, 'number')
    get_scores(potential_codons, 'at_count')
    # Compute the weighted average between count frequency and lowest G/C content
    potential_codons['overall_score'] = potential_codons.apply(
        lambda row: (row['number_score'] * (1.0 - args.weight) + row['at_count_score'] * args.weight), axis=1)
    idx = potential_codons['overall_score'].idxmax()
    # Search back and find what codons corresponds to this index with maximum overall score
    better_codons += codon_to_acid_with_freq_and_at_count.loc[idx].codon

better_seq = Seq(better_codons)


def print_separator():
    """
    This simply prints a bar that is the width of the current terminal.
    """
    print('=' * shutil.get_terminal_size((80, 20)).columns, end='')


g_count, c_count, a_count, t_count = (
    better_seq.count('G'), better_seq.count('C'), better_seq.count('A'), better_seq.count('T'))

print_separator()
print(
    f'Optimizing for {args.weight:.2%} A/T count (decrease G/C) and {1.0 - args.weight:.2%} frequency in human genome')
print_separator()
print(better_seq)
print_separator()
print(
    f'Total Nucleotides: {len(better_seq)}. G Count: {g_count}, C Count: {c_count}, A Count: {a_count}, T Count: {t_count}')
print_separator()
gc_count = g_count + c_count
print(f'G/C Count: {gc_count}')
print(f'G/C Percentage: {gc_count / len(better_seq):.2%}')
print_separator()
output_path = args.output or f'{args.seq.replace(".fasta", "")}_optimized.fasta'
# Make directories leading up if they do not exist already, writing the record will fail otherwise
os.makedirs(os.path.dirname(output_path), exist_ok=True)
record = SeqRecord(better_seq, id=ngly1_raw.id, description=ngly1_raw.description)
SeqIO.write(record, output_path, 'fasta')
print(f'Wrote output file to {output_path}')
print_separator()
