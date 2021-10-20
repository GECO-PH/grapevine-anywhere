"""
ORIGINALLY FROM DATAFUNK 0.1.0 mask.py
Changed d[regex] to d[pos] in parse_mask_file function
as it better suited our mask file.
"""

from Bio import SeqIO
import re, sys


def parse_mask_file(file):
    """
    input is in the format:
    start (1-based), mask character, regex-format string to match record.id
    e.g.:
    13402,?,^Belgium/
    d is a dictionary with the regex strings as keys and position,
    mask character and compiled regular expression as values.
    it has the same number of entries as lines in file
    """

    d = {}

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(',')
            pos, mask_char, regex = l

            d[pos] = {'pos': int(pos),
                        'mask_char': mask_char,
                        'regex': re.compile(regex)}

    return(d)


def mask(fasta_in, fasta_out, mask_file):

    if fasta_out:
        out = open(fasta_out, 'w')
    else:
        out = sys.stdout

    mask_info = parse_mask_file(mask_file)

    input = SeqIO.parse(fasta_in, 'fasta')

    for record in input:
        ID = record.id
        seq = str(record.seq)

        for entry in mask_info:
            regex = mask_info[entry]['regex']

            if re.search(regex, ID):
                pos = mask_info[entry]['pos']
                mask_char = mask_info[entry]['mask_char']

                seq = seq[:pos - 1] + mask_char + seq[pos:]

        out.write('>' + ID + '\n')
        out.write(seq + '\n')

    if fasta_out:
        out.close()
    pass


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    mask_file = sys.argv[3]
    mask(input_file, output_file, mask_file)
