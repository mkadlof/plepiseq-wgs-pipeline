#!/bin/bash

# Initializing parameters
## Stage 1. External files and databases without defaults

reads="" # -r 
primers_id="" # -p
adapters_id="TruSeq3-PE-2" # -a

pangolin_dir="" # -w
nextclade_dir="" # -x
kraken2_dir="" # -y
freyja_dir="" # -z

results_dir="`pwd`/results" # -o
project_dir="" # -m 

## Stage 2. Algorithm parameters with defaults
quality_initial=5 # -b
quality_snp=15 # -j

length=90 # -c
max_number_for_SV=200000 # -d
max_depth=600 # -e
min_cov=20 # -f 
mask=19 # # -i
pval=0.05 # -k
lower_ambig=0.45 # -l
upper_ambig=0.55 # -n
window_size=50 # -s
mapping_quality=30 # -u 

## Stage 3. Other parameters with defaults, only some can be set via context menu
cpu=1 # -t
container_name="" #-g

memory=2024

# Creating help message
echo ${cpu}

usage() {
    echo -e "Skrypt do analizy danych z sekwencjonowania wirusa SARS-CoV-2 z wykorzystaniem platformy Illumina\nNextflow-based pipeline for the analysis of SARS-CoV-2 illumina sequencing data\n"
    echo
    echo -e "Opcje/Options:"
    echo -e "-h                                Wywołanie pomocy/Help message"
    echo
    echo -e "Pliki wejściowe/Input data"
    echo -e "-r\t[format: fastq.gz, example: '/some/directory/*_R{1,2}.fastq.gz']"
    echo -e "  \tŚcieżka do katalogu z wynikami sekwencjonowania wraz z wzorcem nazewnictwa plików"
    echo -e "  \tPath to sequencing data with naming pattern for sequencing files. Use single quotes for this argument"
    echo -e "-p\t[format: string, example: V4.1, expected: one of V1 V1200 V1201 V2 V3 V4.1 V4 V5.3.2 VarSkip2]"
    echo -e "  \tNazwa schematu amplikonów użyta w eksperymencie"
    echo -e "  \tName of amplicon scheme used to amplify genetic materia"
    echo -e "-a\t[format: string, example: TruSeq3-PE-2, default: TruSeq3-PE-2, expected: one of NexteraPE-PE TruSeq2-PE TruSeq2-SE TruSeq3-PE-2 TruSeq3-PE TruSeq3-SE]"
    echo -e "  \tNazwa schematu amplikonów użyta w eksperymencie"
    echo -e "  \tName of amplicon scheme used to amplify genetic materia"
    echo -e "-w\t[format: none, directory, example: /some/directory/]"
    echo -e "  \tŚcieżka do katalogu z baza Pangolin"
    echo -e "  \tPath to a directory with Pangolin data"
    echo -e "-x\t[format: none, directory, example: /some/directory/]"
    echo -e "  \tŚcieżka do katalogu z baza Nextclade"
    echo -e "  \tPath to a directory with Nextclade data"
    echo -e "-y\t[format: none, directory, example: /some/directory/]"
    echo -e "  \tŚcieżka do katalogu z baza dla programu Kraken2"
    echo -e "  \tPath to a directory with data for Kraken2"
    echo -e "-z\t[format: none, directory, example: /some/directory/]"
    echo -e "  \tŚcieżka do katalogu z baza dla programu freyja"
    echo -e "  \tPath to a directory with data for freyja"
    echo -e "-m\t[format: none, directory, example: /some/directory/]"
    echo -e "  \tŚcieżka do katalogu z pobranym z github repozytorium"
    echo -e "  \tPath to a directory with pipeline downloaded from github"
    echo -e "Wyniki/Output data"
    echo -e "-o\t[format: none, directory, example: /some/directory/, default: ./results]"
    echo -e "  \tŚcieżka do katalogu z w, którym zostaną umieszczone wyniki"
    echo -e "  \tPath to a directory for output"
    echo -e "Parametry wewnętrzne, nie trzeba ustawiać"
    echo -e "Program internal parameters. All these parameters have default values"
    echo -e "-b\t[format: int, default: 5, expected 0-60]"
    echo -e "  \tProg jakości nukleotydu (Phred + 33), używany da podprogramów do raportowania jakości danych"
    echo -e "  \tMinimum value for nucleotide quality  (Phred + 33) used to access data quality"
    echo -e "-c\t[format: int, default: 90, expected: 0-151]"
    echo -e "  \tMinimalna długość odczytu wymagana do analizy"
    echo -e "  \tMinimum value for a read to be used in analysis"
    echo -e "-j\t[format: int, default: 15, expected: 0-60]"
    echo -e "  \tProg jakości nukleotydu (Phred + 33), używany podczas identyfikacji SNP/INDEl/SV"
    echo -e "  \tMinimum value for nucleotide quality  (Phred + 33) used during SNPs/INDELs/SVs calling"
    echo -e "-d\t[format: int, default: 200000, expected: at least 1000] "
    echo -e "  \tIlość zmapowanych odczytów używanych do identyfikacji SV"
    echo -e "  \tNumber of mapped pair-end reads used for SV calling"   
    echo -e "-e\t[format: int, default: 600, expected: at least 20]"
    echo -e "  \tOczekiwana wartość pokrycia w genomie po procedurze wygadzania"
    echo -e "  \tExpected coverage after coverage equalization" 
    echo -e "-f\t[format: int, default: 20, expected: at least 5]"
    echo -e "  \tMinimalne pokrycie na danej pozycji wymagane do identyfikacji SNP/INDEL"
    echo -e "  \tMinimum value for coverage at a position required to identify SNPs/INDELs"
    echo -e "-i\t[format int, default: 19, expected: at least 0]"
    echo -e "  \tMaksymalne pokrycie przy którym dana pozycja będzie maskowana symbolem 'N' w genomie, nie może byc większa niz wartość podana dla flagi -f"
    echo -e "  \tMaximum value for coverage value to mask low coverage regions with "N", cannot be greater than value supported with -f"
    echo -e "-k\t[format float, default: 0.05, expected: at most 1]"
    echo -e "  \tWartość p-value poniżej którego wariant traktujemy jako istotny statystycznie"
    echo -e "  \tp-value to determine if a variant is significant"
    echo -e "-l\t[type: float, default: 0.45, expected: 0-1]"
    echo -e "   \tOdsetek odczytów dla pomniejszego wariantu aby na danej pozycji wprowadzić symbol niejednoznacznego nukleotydu"
    echo -e "   \tMinimum ratio of reads carrying a minor allele to introduce ambiguity at a given position"
    echo -e "-n\t[type: float, default: 0.55, expected: 0-1]"
    echo -e "   \tOdsetek odczytów dla pomniejszego wariantu powyżej ktorego nie wprowadzamy symbolu niejednoznacznego nukleotydu"
    echo -e "   \tMinimum ratio of reads carrying a minor allele above which ambiguity symbol is not introduced at a given position"
    echo -e "-s\t[type: int, default: 50, expected: at least 10]"
    echo -e "   \tWielkość okna stosowana podczas wygładzania pokrycia"
    echo -e "   \tWindow size used during coverage equalization"
    echo -e "-u\t[type: int, default: 30, expected: at least 1]"
    echo -e "   \tMinimalna jakość mapowania odczytu do genomu referencyjnego"
    echo -e "   \tMinimum mapping quality"
    echo -e ""
    echo -e "Parametry inne"
    echo -e "-t  [int, default: 1, expected: at least 1]"
    echo -e "  \tIlość CPU/No of CPU threads"
    echo -e "-g  [str, example: sars_illumina_nf:1.0]"
    echo -e "  \tNazwa obrazu docker ze wszystkim programi, z których korzysta program/Docker image name with all the tools used by the pipeline"
}

# Parsing arguments
while getopts ":r:g:p:a:w:x:y:z:m:o:b:c:j:d:e:f:i:k:l:n:s:u:t:h" flag; do
    case  "${flag}" in
        r) reads="${OPTARG}" ;;
        g) container_name="${OPTARG}" ;;
        p) primers_id="${OPTARG}" ;;
        a) adapters_id="${OPTARG}" ;;
        w) pangolin_dir="${OPTARG}" ;;
        x) nextclade_dir="${OPTARG}" ;;
        y) kraken2_dir="${OPTARG}" ;;
        z) freyja_dir="${OPTARG}" ;;
        o) results_dir="${OPTARG}" ;;
        m) project_dir="${OPTARG}" ;;
        b) quality_initial="${OPTARG}" ;;
        c) length="${OPTARG}" ;;
        j) quality_snp="${OPTARG}" ;;
        d) max_number_for_SV="${OPTARG}" ;;
        e) max_depth="${OPTARG}" ;;
        f) min_cov="${OPTARG}" ;;
        i) mask="${OPTARG}" ;;
        k) pval="${OPTARG}" ;;
        l) lower_ambig="${OPTARG}" ;;
        n) upper_ambig="${OPTARG}" ;;
        s) window_size="${OPTARG}" ;;
        u) mapping_quality="${OPTARG}" ;;
        t) cpu=${OPTARG} ;;
        h | *) usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -${OPTARG}"
            exit 1
            ;;
        :)
            echo "Argument -${OPTARG} wymaga podania wartości"
            exit 1
            ;;
    esac
done

# Tests

## Script was called without arguments
if [ $# -eq 0 ]; then
    usage
    exit 1
fi

# Expand the reads pattern and check if there are at least two files
expanded_reads=$(eval ls ${reads} 2> /dev/null)

if [ $(echo "${expanded_reads}" | wc -w) -lt 2 ]; then
    echo -e 'Nie podano ścieżki przy pomocy argumentu -r lub podane pliki nie istnieją\nThere are no files matching pattern provided with -r options'
    exit 1
fi

CORRECT_ID=0
ALL_PRIMERS=(EQA2024.V4_1 EQA2024.V5_3 EQA2023.SARS1 EQA2023.SARS2 V1 V1200 V1201 V2 V3 V4.1 V4 V5.3.2 VarSkip2)
for var in ${ALL_PRIMERS[@]}; do
    if [ ${primers_id} == ${var} ];then
           CORRECT_ID=1
           break
           fi           
done

if [ ${CORRECT_ID} -eq 0 ]; then
    echo -e "Nie podano właściwej wartości dla parametru -p / Dostępne wartości to ${ALL_PRIMERS[@]}\n"
    echo -e "Please specify correct primer scheme name with -p / Available options are ${ALL_PRIMERS[@]}\n"
    exit 1
fi

## Check if directory with pangolin data exists
if [ ! -d "${pangolin_dir}" ]; then
     echo -e 'Nie podano argumentu -w / Podany katalog nie istnieje\n'
     echo -e 'Please specify path to pangolin data with -w / Provided directory does not exist\n'
     exit 1
fi

## Check if directory with nextclade data exists
if [ ! -d "${nextclade_dir}" ]; then
         echo -e 'Nie podano argumentu -x / Podany katalog nie istnieje\n'
         echo -e 'Please specify path to nextclade data with -x / Provided directory does not exist\n'
     exit 1
fi

## Check if a directory with data for kraken2 exists
if [ ! -d "${kraken2_dir}" ]; then
         echo -e 'Nie podano argumentu -y / Podany katalog nie istnieje\n'
         echo -e 'Please specify path to data for kraken2 with -y / Provided directory does not exist\n'
     exit 1
fi

## Check if a directory with data for freyja exists
if [ ! -d "${freyja_dir}" ]; then
         echo -e 'Nie podano argumentu -z / Podany katalog nie istnieje\n'
         echo -e 'Please specify path to data for freyja with -z / Provided directory does not exist\n'
     exit 1
fi

## Print warning if directory where results will be placed already exist
if [ -d "${results_dir}" ]; then
    echo -e 'Uwaga katalog, który podano jako lokalizacje wyjścia istnieje, część plików może zostać nadpisana\n'
    echo -e "Warning the output directory already exists, some files may be overwritten"
fi

## Check if directory with nextflow modules exists
if [ ! -d "${project_dir}" ]; then
         echo -e 'Nie podano argumentu -m / Podany katalog nie istnieje\n'
         echo -e 'Please specify path to directory with modules using -m / Provided directory does not exist\n'
         exit 1
fi

## Check quality initial value
if [ ${quality_initial} -lt 0 ] || [ ${quality_initial} -gt 60 ] ; then
        echo -e 'Wartość podawana jako parametr do flagi -b musi przyjmować wartości miedzy 0 a 60\n'
    echo -e "Value provided with -b must be between 0 and 60" 
        exit 1
fi

# Check read length value
if [ ${length} -lt 0 ] || [ ${length} -gt 151 ] ; then
        echo -e 'Wartość podawana jako parametr do flagi -c musi przyjmować wartości miedzy 0 a 151\n'
    echo -e "Value provided with -c must be between 0 and 151"
        exit 1
fi

## Check quality provided for SNP calling
if [ ${quality_snp} -lt 0 ] || [ ${quality_snp} -gt 60 ] ; then
        echo -e 'Wartość podawana jako parametr do flagi -j musi przyjmować wartości miedzy 0 a 60\n'
    echo -e "Value provided with -j must be between 0 and 60"
        exit 1
fi

## Check if number of reads for SV calling is sufficiently height
if [ ${max_number_for_SV} -le 1000 ]; then
    echo -e 'Wartość podawana jako parametr do flagi -d musi przyjmować wartość co najmniej 1000\n'
    echo -e "Value provided with -d must be between greater than 1000"
        exit 1
fi

## Check if coverage after equalization is above 50
if [ ${max_depth} -le 49 ]; then
        echo -e 'Wartość podawana jako argument do flagi -e musi przyjmować wartości co najmniej 50\n'
    echo -e "Value provided with -e must be at least 50"
        exit 1
fi

## Check if number of reads supporting a variant is a t least 5
if [ ${min_cov} -le 4 ]; then
        echo -e 'Wartość podawana jako parametr do flagi -f musi przyjmować wartości co najmniej 5\n'
    echo -e "Value provided with -e must be at least 5"
        exit 1
fi

## Check if low coverage flag is non negative
if [ ${mask} -lt 0 ]; then
        echo -e 'Wartość podawana jako parametr do flagi -i musi przyjmować wartości większe od 0\n'
    echo -e "Value provided with -i must be at least 0"
        exit 1
else
    mask=`echo "${mask} + 1" | bc -l`
fi

# Check if SNPs are not called in a regions that will be mask in the end
if [ ${min_cov} -lt  ${mask} ]; then
        echo -e 'Wartość podawana jako parametr do flagi -f nie mogą byc mniejsze niz wartości podane do flagi -i'
    echo -e "Value provided with -f must be greater thn one provided with -i"
        exit 1
fi

## Check if ratio for lower ambiguity is between 0-1
if awk "BEGIN {exit !(${lower_ambig} > 1 || ${lower_ambig} < 0) }"; then
        echo -e "Wartość podawana jako parametr do flagi -l musi mieścić sie w przedziale 0-1"
    echo -e "Value provided with -l must be between 0-1"
        exit 1
fi

## Check if ratio for upper ambiguity is between 0-1
if awk "BEGIN {exit !(${upper_ambig} > 1 || ${upper_ambig} < 0) }"; then
        echo -e "Wartość podawana jako parametr do flagi -n musi mieścić sie w przedziale 0-1"
    echo -e "Value provided with -n must be between 0-1"
        exit 1
fi

if awk "BEGIN {exit !(${upper_ambig} < ${lower_ambig}) }"; then
        echo -e "Wartość podawana jako parametr flagi -u nie może byc mniejsza niz wartość podana dla flagi -w"
        exit 1
fi

if [ ${window_size} -lt 10 ]; then
        echo -e "Wartość podawana jako parametr do flagi -s musi przyjmować wartości co najmniej 10"
        echo -e "Value provided with -s must be at least 10"
        exit 1
fi

if [ ${mapping_quality} -lt 1 ]; then
        echo -e "Wartość podawana jako parametr do flagi -u musi przyjmować wartości co najmniej 1"
        echo -e "Value provided with -u must be at least 1"
        exit 1
fi

if [ ${cpu} -lt 1 ]; then
        echo -e "Wartość podawana jako parametr do flagi -t musi przyjmować wartości co najmniej 1\n"
        echo -e "Value provided with -t must be at least 1"
        exit 1
fi

# Test if image exists in the system

# Calling pipeline
nextflow run ${project_dir}/nf_pipeline.nf \
    --reads ${reads} \
    --primers_id ${primers_id} \
    --adapters_id ${adapters_id} \
    --pangolin_db_absolute_path_on_host ${pangolin_dir} \
    --nextclade_db_absolute_path_on_host ${nextclade_dir} \
    --kraken2_db_absolute_path_on_host ${kraken2_dir} \
    --freyja_db_absolute_path_on_host ${freyja_dir} \
    --results_dir ${results_dir} \
    --modules ${project_dir}/modules \
    --quality_initial ${quality_initial} \
    --length ${length} \
    --quality_snp ${quality_snp} \
    --max_number_for_SV ${max_number_for_SV} \
    --max_depth ${max_depth} \
    --min_cov ${min_cov} \
    --mask ${mask} \
    --lower_ambig ${lower_ambig} \
    --upper_ambig  ${upper_ambig} \
    --memory ${memory} \
    --window_size ${window_size} \
    --mapping_quality ${mapping_quality} \
    --threads ${cpu} \
    -with-docker ${container_name} \
    -with-trace

### Example call ###
# $ ./run_nf_pipeline.sh -r '/absolute/path/to/reads/*_{1,2}.fastq.gz' -p EQA2024.V5_3 -w /absolute/path/to/pangolin/ -x /absolute/path/to/nextclade/ -y /absolute/path/to/kraken2/ -z /absolute/path/to/freyja/ -m /absolute/path/to/pipeline/root/nf_illumina_sars -t 10
