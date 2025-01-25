# skrypt do przygotowania bazy VFDB
# rozbija wielka faste na podfasty umieszczane w katalogu
# organizm/Virulence_class/Virulence_factor/gene.fa
# katalog Virulence_factor moze zawierac wiele genow
# oczywiscie Virulence_class i Virulence_factor to zmienne

from Bio import SeqIO
import re
import sys
import os
import subprocess

def execute_command(polecenie: str):
    """

    :param polecenie: str, polecenie do wykonania w powloce
    :return: Polecenie ma zwrocic zwykle jakis plik, wiec sama funkcja nie zwraca nic
    """
    try:
        wykonanie = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)
        wykonanie.communicate()
        return True
    except:
        return False

# najpierw pobieramy naglowki
def read_fasta(fasta):
    """
    Funkcja do wczytania segemtnu
    :param fasta:
    :return:
    """
    slownik = {}
    fasta = SeqIO.parse(fasta, "fasta")
    for read in fasta:
        slownik[f'>{read.description}'] = str(read.seq)
    return slownik

def extract_info_from_header(naglowek, bis = 0):
    """

    W VFDB naglowki dzieki bogu maja regularna nazwe
    # wyciagamy z nich informacje o organizmie, VF, VFC, VF opis i VFC opis
    # wszystko zapisujemy do slownika gdzie naglowek jest kluczem a kolejne wartosci sa wlscie
    # Generalnie poslugujemy sie troche wyrazeniami regularnymi
    :param naglowek:
    :return:
    """
    # Pattern napisany przy pomocy strony https://regexr.com/7voga
    # dla stringa
    # >VFG030121(gb|WP_013988985) (kefB) cation:proton antiporter [Potassium/proton antiporter (VF0838) - Immune modulation (VFC0258)] [Mycobacterium africanum GM041182]
    # wyciaga 7 grup w kolejnosci
    # match.group(1)
    # 'gb|WP_013988985'
    # match.group(2)
    # 'kefB'
    # match.group(3)
    # 'Potassium/proton antiporter'
    # match.group(4)
    # 'VF0838'
    # match.group(5)
    # 'Immune modulation'
    # match.group(6)
    # 'VFC0258'
    # match.group(7)
    # 'Mycobacterium '
    # testowano tez na
    # >VFG021469(gb|WP_000062035) (stbA) type 1 fimbrial protein [Stb (VF0954) - Adherence (VFC0001)] [Salmonella enterica subsp. enterica serovar Heidelberg str. SL476]
    pattern = '^>\w+\((\w+\|\w+)\)\s\((\w+)\)\s.+\s\[(.+)\s\((.+)\)\s-\s(.+)\s\((.+)\)\]\s\[(\w+)\s.*\]'

    # update patterny gdy naglowek ma gen ww trybie (sycN/vcr2) lub (vgrG-2)
    # pattern = '^>\w+\((\w+\|\w+|\w+\|\w+\.\w+)\)\s\((\w+|\w+\/\w+|\w+-.+)\)\s.+\s\[(.+)\s\((.+)\)\s-\s(.+)\s\((.+)\)\]\s\[(\w+)\s.*\]'
    pattern = '^>\w+\((\w+\|\w+|\w+\|\w+\.\w+)\)\s\((.*)\)\s.+\s\[(.+)\s\((.+)\)\s-\s(.+)\s\((.+)\)\]\s\[(\w+)\s.*\]'
    #obejscie dla sekwencji gdzie id jest bledne np
    #  >VFG000371 (yadA) trimeric autotransporter adhesin YadA [YadA (VF0133) - Effector delivery system (VFC0086)] [Yersinia pestis CO92]
    if bis:
        pattern = '^>(\w)+\s\((.*)\)\s.+\s\[(.+)\s\((.+)\)\s-\s(.+)\s\((.+)\)\]\s\[(\w+)\s.*\]'
    match = re.match(pattern, naglowek)
    return [match.group(1), match.group(2), match.group(3), match.group(4),match.group(5),match.group(6),match.group(7)]

if __name__ == '__main__':

    if os.path.exists('VFDB_setB_nt.fa'):
        # found some previous data remove everything before we do main script
        _ = [shutil.rmtree(dir_path) for dir_path in os.listdir('.') if dir_path not in ['prep_VFDB.py', 'README'] and os.path.isdir(dir_path)]
        
        if os.path.exists('VFs.xls'):
            os.remove('VFs.xls')
    
        os.remove('VFDB_setB_nt.fa')
        
    # download the new data 
    execute_command('curl -u anonymous:anonymous www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz -O')
    execute_command('curl -u anonymous:anonymous www.mgc.ac.cn/VFs/Down/VFs.xls.gz -O')
    
    #unzip
    execute_command('gunzip VFDB_setB_nt.fas.gz')
    execute_command('gunzip VFs.xls.gz')

    fasta_file = 'VFDB_setB_nt.fas'
    slownik_sekwencji = read_fasta(fasta_file)

    # jak to zrobic
    # of to jest glupie, ale zadzial
    # idac po naglowkach tworzymy slownik gdzie kluczem jest  string w postaci
    # {organizm}/{VFC}/{VF}
    # potem jest kolejny klucz czyli {gene}
    # i wartosciami tego slownika jest ostateczny slownik gdzie kluczem jest id sekwencji a wartosci sekwencja
    # na koniec iterujemy po slowniku i zapisujemy faste do pliku ...
    # na przyklad
    # slownik['salmonella/VFC1/VF13'] ma podslowniki
    # slownik['salmonella/VFC1/VF13']['genA'] slownik['salmonella/VFC1/VF13']['genB']
    # i to ma kolejne podslowniki
    # slownik['salmonella/VFC1/VF13']['genA'] = {'some_id':'sekwencja'; 'different_id': 'sekwencja}
    # pierwszy klucz bedzie katalogiem
    # drugi klucz nazwa pliku w tym katalogu
    # a trzeci poslownik jest zapisywany do tego pliku

    slownik_VFDB = {}

    for id,seq in slownik_sekwencji.items():
        try:
            seq_id, gene, VF_name, VF_id, VFC_name, VFC_id, org = extract_info_from_header(id)
        except:
            try:
                seq_id, gene, VF_name, VF_id, VFC_name, VFC_id, org = extract_info_from_header(id, 1)
            except:
                print(f'Omijam sekwencje {id}')
        gene = gene.replace('/', '-')
        pattern = f'{org}/{VFC_id}/{VF_id}'
        if pattern not in slownik_VFDB:
            slownik_VFDB[pattern] = {}

        if ')' in gene:
            gene=gene.split(')')[0]
        if '(' in gene:
            gene=gene.split('(')[0]
        #usuwamy dziwne znaki w nazwach genow
        gene=gene.replace('*', '').replace("'", "").replace(" ","")
        if gene not in slownik_VFDB[pattern]:
            slownik_VFDB[pattern][gene] = {}

        slownik_VFDB[pattern][gene][id] = seq

    for dir in slownik_VFDB.keys():
        os.makedirs(dir)
        for geny in slownik_VFDB[dir]:
            with open(f'{dir}/{geny}.fa', 'w') as f:
                for id,seq in slownik_VFDB[dir][geny].items():
                    f.write(f'{id}\n{seq}\n')

    # running makeblastdb to index ale the file
    print('Running makeblastdb')
    execute_command('find . -name "*fa" | xargs -I {} --max-procs=96  bash -c "makeblastdb -in {} -dbtype nucl"')
