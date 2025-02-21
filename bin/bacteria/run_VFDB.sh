#!/bin/bash
# ver 6 oparta

GENOM=$1
MAX_PROC=$2
ORG=$3 # Nazwa organizmu akceptujemy tylko te gatunki ktore obecne sa w bazie VFDB w konkretnych katalogach
       # np. Salmonella. Tylko blastujemy sekwencje z tego katalogu
       # UWAGA baza musi byc zamonotwana w katalogu /db wiec szukamy plikow fasta
       # w /db/${ORG}
PIDENT=$4 # minimalne seq id aby uznac ze target jest ok # powiedzmy 80 ?
EVALUE=$5 # max e-value aby uznac ze target jest ok # powiedzmy 0.01
MIN_COV=$6 # minimalny procent dlugosc i query i target aby uznac ze mapowanie jest ok # powiedzmy 80%

run_blastn () {
	# zwracamy plik tabularyczny w ktorym podajemy kolejnoquery
	# 1. identyfikator query
	# 2. identyfikator target
	# 3. poczatek mapownia quey
	# 4. koniec mapowani query
	# 5. dlugosc query
	# 6. poczatek, 7. koniec mapowania subject i 8. dlugosc subject
	# 9. identycznosc sekwencyjna alignmentu i 10. e-value
	FILENAME=`basename ${1} | cut -d "." -f1`
	VFNAME=`dirname ${1} | cut -d "/" -f5`
	VFCNAME=`dirname ${1} | cut -d "/" -f4`
	
	/blast/bin/blastn -query ${2} -db ${1} -out blastn_${VFCNAME}_${VFNAME}_${FILENAME}.tab -outfmt "6 qseqid saccver qstart qend qlen sstart send slen pident evalue " -max_target_seqs 5 
	return 0
}

export -f run_blastn
find /db/${ORG} -name '*fa'  >> all_files.txt

cat all_files.txt | xargs -I {} --max-procs=${MAX_PROC} bash -c "run_blastn {} ${GENOM}"

echo 'Analizuje wyniki'
for K in `ls blastn_*`
do
	VFC=`echo ${K} | cut -d "_" -f2 | cut -d "." -f1`
	VF=`echo ${K} | cut -d "_" -f3 | cut -d "." -f1`
	GEN=`echo ${K} | cut -d "_" -f4- | cut -d "." -f1`
	mv ${K} blastn_tmp.tab

	# Szukamy hitow ktore spelniaja kryteria identycznosci sekwencyjnrj, pokrycia alignmentem i e-value	
	cat blastn_tmp.tab | awk -v MAX_ID=${PIDENT} -v MIN_COV=${MIN_COV} -v MAX_EVALUE=${EVALUE} '{if($9 >= MAX_ID && ($4-$3)/$5*100 >= MIN_COV && ($7-$6)/$8*100 >= MIN_COV && $10 <= MAX_EVALUE ) print $0}' > tmp

	ILE=`cat tmp | cut -f1 | sort | uniq | wc -l`
	
	if [ ${ILE} -eq 0 ]; then
		echo -e "Dla ${K} zapisuje ${ORG}\t${VF}\t${VFC}\t${GEN}\tBRAK"
		echo -e "${ORG}\t${VF}\t${VFC}\t${GEN}\tNo_hits\tO\t0\t0\tNo_ref\tBRAK" >> VFDB_summary.txt
	elif [ ${ILE} -eq 1 ]; then
		NAZWA=`cat tmp | cut -f1 | head -1`
		REF_NAME=`cat tmp | sort -rnk9 | head -1 | cut -f2 | tr "(" " " | cut -d " " -f2 | tr -d ")"`
		SEQ_ID=`cat tmp | sort -rnk9 | head -1 | cut -f9`
		SLEN=`cat tmp | sort -rnk9 | head -1 | cut -f8`
		QLEN=`cat tmp | sort -rnk9 | head -1 | cut -f5`
		QLEN_END=`cat tmp | sort -rnk9 | head -1 | cut -f4`
		QLEN_START=`cat tmp | sort -rnk9 | head -1 | cut -f3`
		COV=`echo "(${QLEN_END} - ${QLEN_START} + 1)/${QLEN} * 100" | bc -l`
		#echo -e "Dla ${K} zapisuje ${ORG}\t${VF}\t${VFC}\t${GEN}\t${NAZWA}"
		echo -e "${ORG}\t${VF}\t${VFC}\t${GEN}\t${NAZWA}\t${SEQ_ID}\t${COV}\t${REF_NAME}\tONE_HIT" >> VFDB_summary.txt

	elif [ ${ILE} -gt 1 ]; then
		NAZWA=`cat tmp | sort -rnk9 | head -1 | cut -f1` # Zwracamy tylko najlepszy hitm jesli jest wiele hitow "ze 100 seq id" to random
		REF_NAME=`cat tmp | sort -rnk9 | head -1 | cut -f2 | tr "(" " " | cut -d " " -f2 | tr -d ")"`
		SEQ_ID=`cat tmp | sort -rnk9 | head -1 | cut -f9`
                SLEN=`cat tmp | sort -rnk9 | head -1 | cut -f8`
		QLEN=`cat tmp | sort -rnk9 | head -1 | cut -f5`
                QLEN_END=`cat tmp | sort -rnk9 | head -1 | cut -f4`
                QLEN_START=`cat tmp | sort -rnk9 | head -1 | cut -f3`
		COV=`echo "(${QLEN_END} - ${QLEN_START} + 1)/${QLEN} * 100" | bc -l`
		# echo -e "Dla ${K} zapisuje wielohitowej ${ORG}\t${VF}\t${VFC}\t${GEN}\t${NAZWA}"
		echo -e "${ORG}\t${VF}\t${VFC}\t${GEN}\t${NAZWA}\t${SEQ_ID}\t${COV}\t${REF_NAME}\tMULTIPLE_HITS" >> VFDB_summary.txt
	fi
	rm blastn_tmp.tab
	rm tmp
done
