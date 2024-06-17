#!/usr/bin/env python

import argparse
import textwrap

from Bio import SeqIO

from modeller import Alignment, Environ, log, Model


def clean_seq(seq: str) -> str:
    """Clen sequence from unnecessary characters '-' and '*'"""
    seq = seq.replace('-', '')
    if seq.endswith('*'):
        return seq[:-1]
    return seq


def prepare_target_seq(input: str) -> str:
    record = next(SeqIO.parse(input, 'fasta'))
    seq = clean_seq(str(record.seq))
    seq = seq[27-1:1146+1-7]  # Hardcoded trim values 27-1146 for 7dwz
    seq = textwrap.fill(seq, width=60)
    title = f'>P1;target\n'
    description = f'sequence:target::::::::\n'
    return title + description + seq + '/\n' + seq + '/\n' + seq + '*\n\n'


def write_initial_alignment_file(target, template, output):
    with open(output, 'w') as f:
        f.write(target)
        f.write(template)


def main():
    parser = argparse.ArgumentParser(description='Prepare alignment')
    parser.add_argument('input', type=argparse.FileType('r'), help='Input FASTA file')
    args = parser.parse_args()
    target = prepare_target_seq(args.input)
    template = "\n"
    write_initial_alignment_file(target, template, 'alignment.pir')

    log.level(output=1, notes=1, warnings=1, errors=1, memory=0)
    env = Environ()
    aln = Alignment(env)
    mdl = Model(env)

    aln.append(file='alignment.pir', align_codes='all')
    mdl.read(file='7dwz.pdb', model_segment=('FIRST:A', 'LAST:C'))
    aln.append_model(mdl, align_codes='7dwz', atom_files='7dwz.pdb')
    aln.align()
    aln.write(file='alignment.pir', alignment_format='PIR')


if __name__ == '__main__':
    main()
