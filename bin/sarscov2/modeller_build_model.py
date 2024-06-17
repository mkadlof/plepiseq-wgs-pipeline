#!/usr/bin/env python

import argparse

from modeller import Environ, log
from modeller.automodel import AutoModel, refine


def run_modeller(alignment_pir: str, mode: str):
    env = Environ()
    log.level(output=1, notes=1, warnings=1, errors=1, memory=0)

    a = AutoModel(env, alnfile=alignment_pir,
                  knowns=('7dwz'), sequence='target')
    if mode == 'fast':
        a.very_fast()  # prepare for extremely fast optimization
    else:
        a.md_level = refine.slow
    a.starting_model = 1
    a.ending_model = 1
    a.make()


def main():
    parser = argparse.ArgumentParser(description='Prepare alignment')
    parser.add_argument('alignment_pir', help='Input FASTA file')
    parser.add_argument('--mode', default='fast', choices=['slow', 'fast'],
                        help='Model refinement mode')
    args = parser.parse_args()
    run_modeller(args.alignment_pir, args.mode)


if __name__ == '__main__':
    main()
