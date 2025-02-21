'''
To prosty skrypt do generowania pokrycia contigow
i zachowywania tylko contigow o wywstarczajaco wysokim pokryciu w stosunku do sredniego pokrycia

'''

import pysam
import sys
import numpy as np
from Bio import SeqIO

def read_genome_boundaries(fasta):
    """
    Wprawdzie w influenza tez korzystamy z primerow ale sa one wspolne dla wszystkich
    segmentow i usuwane wczesniej przez cutadapt-a. Na tym etapie chcemy po porstu dostac info jakie koordynaty
    maja segmenty np "chr1_PA1" ma [0,2000] TO JEST PRZYKLAD NIE RZECZYWISTE WARTOSCI
    :param fasta: sciezka do pliku fasta z genmem
    :return:
    """
    dlugosci_segmentow= {}
    genome = sequences = SeqIO.parse(fasta, "fasta")
    for segment in genome:
        dlugosci_segmentow[segment.id] =  len(segment.seq)
    return dlugosci_segmentow

def filter_fastas(input_fasta, output_fasta, slownik_pokryc, global_coverage, threshold, threshold_coverage):
    """
    Skrypt zapisuje do pliku output_fasta tylko te contigi z pliku input_fasta,
    ktorych pokrycie (trzymane w slowniku slownik_pokryc) wynosi co najmniej threshold * global_coverage
    oraz jest wieksze niz threshold_coverage
    De facto wprowadzenie tego warunku wywala wszystkie contigi z pokryciem mniejszym niz
    threshold_coverage
    Outputowa fasta ma dodawana w naglowku _cov_X, gdzie X to wartosc
    """
    record= SeqIO.parse(input_fasta, 'fasta')
    with open(output_fasta, 'w') as plik, open('Rejected_contigs.fa', 'w') as reject:
        for r in record:
            if slownik_pokryc[r.id] > (global_coverage * threshold) and slownik_pokryc[r.id] > threshold_coverage:
                plik.write(f'>{r.id}_cov_{slownik_pokryc[r.id]}\n')
                plik.write(f'{str(r.seq)}\n')
            else:
                reject.write(f'>{r.id}_cov_{slownik_pokryc[r.id]}\n')
                reject.write(f'{str(r.seq)}\n')
                pass
    return True

### main program ###

reference_genome=sys.argv[1] # plik z contigami
plik_bam = sys.argv[2] # bam zmapowany na contigi
threshold = float(sys.argv[3]) # minimalna frakcja sredniego pokrycia contigu w stosunku do pokrycia globalnego
threshold_coverage = float(sys.argv[3])
qual = 5 # ignorujemy slabe zasady, ustawiamy "na sztywno" bo nie jest to bardzo wazny parametr w tym przypadku

dlugosci_segmentow = read_genome_boundaries(reference_genome)

samfile = pysam.AlignmentFile(plik_bam, "rb")

## Liczenie sredniego "globalnego" pokrycia oraz sredniego pokrycia dla kazdego z contigow ##
global_coverage = [] # w tej liscie trzymamy pokrycia kazdej pozycji 
slownik_pokryc = {}

for segment in dlugosci_segmentow:
    segment_coverage = [] # lista do trzymania pokryc danego segmentu
    pokrycia = samfile.count_coverage(contig=segment, start=0, stop=dlugosci_segmentow[segment], quality_threshold=qual)
    # pokrycia to 4 elementowy tupple z wartosciami na danej pozycji A C T i G
    for i in range(len(pokrycia[0])):
        suma = pokrycia[0][i] + pokrycia[1][i]  + pokrycia[2][i]  + pokrycia[3][i]
        global_coverage.append(suma)
        segment_coverage.append(suma)

    slownik_pokryc[segment] = np.mean(segment_coverage)

## Zapisywanie outputu ##
filter_fastas(reference_genome, 'final_scaffold_filtered.fa', slownik_pokryc, np.mean(global_coverage), threshold, threshold_coverage)

