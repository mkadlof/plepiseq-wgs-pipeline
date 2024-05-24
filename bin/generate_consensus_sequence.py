#!/usr/bin/env python3

import sys

from Bio import AlignIO
from Bio.Align import AlignInfo

alignment = AlignIO.read(sys.argv[1], "fasta")
summary_align = AlignInfo.SummaryInfo(alignment)
with open('output_consensus.fa', 'w') as f:
    f.write('>consensus\n')
    f.write(str(summary_align.gap_consensus(threshold=0.6, ambiguous='X')))
