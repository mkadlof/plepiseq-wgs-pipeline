#!/bin/bash
# ver 10

# W wersji 8 dalej opracowywanej na podstawie przykladu ST8 zmieniono if-a przy multimapowani
# Nie liczy sie dlugosc a tylko numer allelu, im nizszy tym lepszy

# W wersji 9
# mapowanie "kolejnego" allelu nie moze zaczynac sie PRZED mapowaniem poprzedniego
# zachowujemy jednak zpojnosc na contig 


# w wersji 10 chyba jednak upraszam kolejny alle ma preferncje na bycie z contgia poprzedniego allelu, ale bez patrzeania czy przed/za

# chyba entero tez sprawdza czy kolejny allel zaczyba sie ZA poprzednim

# skrypt zwraca plik log.log
# w ktrym allele podzielone sa na 4 grupy
# normal - jest tylko jeden allel ze 100% zgodnoscia w genomie, i ten allel jest obecny w genomie w 'pelnej' dlugosci
# short - jest tylko jeden allel ze 100% zgodnoscia w genomie, ale ten allel nie jest caly w genomie /jest uciety/
# multi - jest wiele alleli ze 100 seq ident, wybieramy allel ktory jest pelny i najdluzszy
# unk - nie ma allelu ze 100 seq ident, wybieramy allel ktory ma najwieksza wartosc score = dlugosc allelu w genommie - ilosc mismatch - ilosc gap

# Dla kazdej z tych klas log.log zawiera troche inne dodatkowe kolumny 5+
# normal - dlugosc allelu w cgmlst i w genomie (de facto musza to byc te same liczby)
# multi - lista wszystkich alleli ze 100% seq id
# short dlugosc allelu w cgmlst i w genomie (de facto musza to byc te same liczby, wtedy widac o ile za krtki jest allel w genomie)
# unk -  dlugosc allelu w cgmlst i score 




GENOM=$1
MAX_PROC=$2
SCIEZKA_DO_FASTA=$3
run_blastn () {
	if [ ${1} == 'STY1774' ]; then
		# przerazlwie wolny ale ten allel jest ekstremalnie krotki
		/blast/bin/blastn -query ${2} -db ${3}/${1}.fasta -out blastn_tmp_${1}.tab -outfmt "6 saccver slen qstart qend sstart send pident evalue mismatch gaps qseqid sstrand" -max_target_seqs 5 -word_size 5 
	else
		/blast/bin/blastn -query ${2} -db ${3}/${1}.fasta -out blastn_tmp_${1}.tab -outfmt "6 saccver slen qstart qend sstart send pident evalue mismatch gaps qseqid sstrand" -max_target_seqs 5
	fi
	return 0
}

export -f run_blastn

ls ${SCIEZKA_DO_FASTA}/*fasta | grep -v all | xargs -I {} basename {} | cut -d "." -f1 >> tmp.txt
cat tmp.txt | xargs -I {} --max-procs=${MAX_PROC} bash -c "run_blastn {} ${GENOM} ${SCIEZKA_DO_FASTA}"

PREVIOUS_POS="-1"
PREVIOUS_CONTIG="NONE"
CONTIG_DIR="plus"

# jesli contig jest plus kolejne allele z cg musza byc rowniez plus oraz musza byc ZA poprzednim allelem
# jesli contig jest minus kolejne allele z cg musza byc minus oraz musza byc PRZED poprzednim allelem

touch cgMLST_all_identical_allels.txt # W przypadku gdy mialem genom w ktorym wszystkie loci mialy allel -1, ten plik nie powstawal i wyrzucal pipeline
for K in `cat tmp.txt`
do
	#if [ ${K} == "STMMW_00931" ]; then
	#	exit 1
	#fi
	# zamienic na mv po testach
	cp blastn_tmp_${K}.tab blastn_tmp.tab
	# 1 dostajemy jeden 100 % wynik pod katem seq id, sprawdzamy tylko czy pokrywa sie on w 100% z dlugoscia allelu (nie jest trimowany na koncach)
	#echo ${K}
	HITY=`cat blastn_tmp.tab | awk '{if($7 == 100  &&  ($4 - $3  + 1) == $2) print $0}' | wc -l`
	# Inicjalizacja zmiennych dla allelu
	FINAL_ID=`echo -1`
        FINAL_TRUE_LEN=`echo -1`
        FINAL_QUERY_LEN=`echo -1`
	
	# Inicjalizacja zmiennych ciag dalszy
	FINAL_HIT_START="-1"
	FINAL_HIT_END="-1"
	FINAL_DIR="NONE"
	FINAL_CONTIG="-1"

	PREFERRED_CONTIG="NONE"

	if [ ${HITY} -gt 0 ]; then
		# hity ze 100% seq id tylko te sa wazen
		ALLELS=(`cat blastn_tmp.tab  | awk '{if($7 == 100  &&  ($4 - $3  + 1) == $2) print $1}'`)
		ILEOSC_HITOW=`echo ${ALLELS[@]} | tr " " "\n" | wc -l`
		echo -e "${K}\t${ILEOSC_HITOW}" >> cgMLST_all_identical_allels.txt
		# Musimy sprawdzic czy mapowania sa tylko na jeden contig
		# Jesli na wiele to preferujemy ten contig na ktory mapowal sie POPRZEDNI allel
		# Jesli na jeden to idziemy "klasycznie"

		HIT_CONTIGS_UNIQUE=`cat blastn_tmp.tab | awk '{if($7 == 100 &&  ($4 - $3  + 1) == $2 )  print $11}' | cut -d "_" -f1 | sort | uniq `
		for CONTIG in ${HIT_CONTIGS_UNIQUE}
			do
				if [ ${CONTIG}  ==  ${PREVIOUS_CONTIG} ]; then
					#echo "kontynuacja"
					# wsrod wynikow dla tego allelu sa mapowania na contig ktory byl
					# wybrany jako zawerajacy "poprzedni" allel
					# bedzie to PREFEREOWANY contig i teraz
					PREFERRED_CONTIG="${CONTIG}"
					break
				fi
			done
		
		# Allel na nowym contigu zeruje zmienne contiga
		if [ ${PREFERRED_CONTIG} == "NONE" ]; then
			CONTIG_DIR="NOVEL"
			#PREVIOUS_POS="-1"
			HIT_DIR=`cat blastn_tmp.tab | awk '{if($7 == 100 && ($4 - $3  + 1) == $2)  print $12}' | sort | uniq`
			if [ ${HIT_DIR} == "minus" ]; then
				PREVIOUS_POS="100000000000"
			elif [  ${HIT_DIR} == "plus" ]; then
				PREVIOUS_POS="-1"
			else
				echo "Nie mozna okreslic kierunku startowego"
			fi
		fi

		#echo ${K} prefereowany contig to ${PREFERRED_CONTIG} o kierunku ${CONTIG_DIR}		
		for ID in ${ALLELS[@]}; do
			
			# wybieramy te allele ktore maja 100% pident, a mapowania obejmuja pelna dlugosc sekwencji tego allelu z cgMLST
			# uwaga ID moze byc w wynikach blasta n-razy wybieramy to ktore mapuje sie na pident 100
			
			#uwaga 2 moze byz tak ze ten sam allel ma 2 razy pident 100, ale rozne dlugosci , wybeiramy ten z najdluzszym query len
			
			# Identyfikator allelu
			ALLEL_ID=`cat blastn_tmp.tab |  awk -v id="${ID}" '{if($7 == 100 && $1 == id)  print $1}' |  sed -r  s'/.*_(.+)/\1/' | head -1 `
			# Dlugosc allelu
			ALLEL_LEN=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id)  print $2}' | head -1`

			# Dlugosc mapowania tego allelu na nasz genom
			QUERY_LEN=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id) print $4 - $3  + 1}' | sort -rnk1 | head -1`
			# Gdzie zaczyna sie mapowanie w contigu
			HIT_START=`cat blastn_tmp.tab | awk -v id="${ID}" '{if( $7 == 100 && $1 == id ) print $3}' | head -1` 

			# Gdzie konczy sie mapowanie w contig
			HIT_END=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id ) print $4}' | head -1`
			
			# Na jaki contig mapuje sie allel w naszym genomie
			HIT_CONTIG=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id)  print $11}' | head -1 | cut -d "_" -f1`

			# Kierunkowosc mapowania allelu na nasz genom
			HIT_DIR=`cat blastn_tmp.tab | awk -v id="${ID}" '{if($7 == 100 && $1 == id)  print $12}' | head -1 `
			
			#echo "Dla allelu ${K} zmienne to ${ALLEL_ID} ${ALLEL_LEN} ${QUERY_LEN} ${HIT_START} ${HIT_END} ${HIT_CONTIG} ${HIT_DIR}"
			
			if [ ${HIT_CONTIG} != ${PREFERRED_CONTIG} ] && [ ${PREFERRED_CONTIG} != "NONE" ]; then
				# omijamy allele jesli sa inny niz "Preffered contig" 
				# i preffered contig nie jest ustawiony na NONE
				#echo "omijam allel ${K} jest w zlym contig"
				continue
			fi
		

			# Na tym etapie nie proceduje odczytow nie mapujacych sie na niepreferowany contig 
			# Co musze zdecydowac
			# # MAM preferowany contig 
			# 1. jesli mam n-wersji allelu mapujacych sie na prefereowany contig
			#	a. musze wybrac ten ktory jest ZA (lub PRZED w zaleznosci od kierunku) poprzednim allelem
		      	#	  jesli jest ich wiele wybrac ten o nizszym numerze
			# # NIE MAM prefereowane contigu 
			# 2. to samo co wyzej ale skoro nie mam prefereowanego allelu po prostu wyvieram 
			#    wersje allelu o nizszym mnumerze i okreslam czy kolejne allele beza ZA nim lub PRZED nim
			# 3. mam n-wersji allelu mapujacych sie na rozne contigi, wtedy w sumie robie to samo co
			#    wyzej wybieram nizszy numer i oceniam kierunek ...

			if  [ ${PREFERRED_CONTIG} != "NONE" ]; then
				#echo "PREVIOUS POS to ${PREVIOUS_POS} a HIT start to ${HIT_START}"	
				if [ ${HIT_CONTIG} == ${PREFERRED_CONTIG} ] && [ ${FINAL_ID} -eq -1 ]; then
					FINAL_ID=`echo ${ALLEL_ID}`
                                	FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                                	FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
                                	FINAL_HIT_END="${HIT_END}"
                                	FINAL_HIT_START="${HIT_START}"
					FINAL_DIR="${HIT_DIR}"
					FINAL_CONTIG="${HIT_CONTIG}"


				elif  [ ${HIT_CONTIG} == ${PREFERRED_CONTIG} ] && [ ${ALLEL_ID} -lt ${FINAL_ID} ] && [ ${ALLEL_LEN} -ge ${FINAL_TRUE_LEN} ];then
				# FINAL_ID nie jest -1 wiec nowy allel musi byc o nizszym numerze niz ten juz zidentyfikowany
				# ale reszta warunkow co do pozycji i kierunku musi byc spelniona contig jest plus
				# ALLEL ma nizszy numer niz to co zidentyfikowlaismy, dlugosc tego allelu jest co najmniej taka sama jak allelu zidentyfikowanego
					FINAL_ID=`echo ${ALLEL_ID}`
                                	FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                               		FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
					FINAL_HIT_END="${HIT_END}"
					FINAL_HIT_START="${HIT_START}"
					FINAL_DIR="${HIT_DIR}"
					FINAL_CONTIG="${HIT_CONTIG}"
				
				elif  [ ${HIT_CONTIG} == ${PREFERRED_CONTIG} ] && [ ${ALLEL_ID} -gt ${FINAL_ID} ] && [ ${ALLEL_LEN} -gt ${FINAL_TRUE_LEN} ];then
                                # Allel ma wyzszy numer ale jest DLUZSZY
					FINAL_ID=`echo ${ALLEL_ID}`
                                        FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                                        FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
                                        FINAL_HIT_END="${HIT_END}"
                                        FINAL_HIT_START="${HIT_START}"
                                        FINAL_DIR="${HIT_DIR}"
                                        FINAL_CONTIG="${HIT_CONTIG}"

				fi
			elif  [ ${PREFERRED_CONTIG} == "NONE" ]; then
				# Nie ma prefereowanego contigu wiec nie ma co sprawdzac czy allel mapuje sie za czy przed
				# allelem wczesniejszym
				if [ ${FINAL_ID} -eq -1 ]; then
                                        FINAL_ID=`echo ${ALLEL_ID}`
                                        FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                                        FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
                                        FINAL_HIT_END="${HIT_END}"
                                        FINAL_HIT_START="${HIT_START}"
					FINAL_DIR="${HIT_DIR}"
					FINAL_CONTIG="${HIT_CONTIG}"

                                elif  [ ${ALLEL_ID} -lt ${FINAL_ID} ] && [ ${ALLEL_LEN} -ge ${FINAL_TRUE_LEN} ]; then
                                        FINAL_ID=`echo ${ALLEL_ID}`
                                        FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                                        FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
                                        FINAL_HIT_END="${HIT_END}"
                                        FINAL_HIT_START="${HIT_START}"
					FINAL_DIR="${HIT_DIR}"
					FINAL_CONTIG="${HIT_CONTIG}"
				elif [ ${ALLEL_ID} -gt ${FINAL_ID} ] && [ ${ALLEL_LEN} -gt ${FINAL_TRUE_LEN} ]; then
					FINAL_ID=`echo ${ALLEL_ID}`
                                        FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
                                        FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
                                        FINAL_HIT_END="${HIT_END}"
                                        FINAL_HIT_START="${HIT_START}"
                                        FINAL_DIR="${HIT_DIR}"
                                        FINAL_CONTIG="${HIT_CONTIG}"

                                fi
			fi


		done
		
		# Na tym etapie albo mam wybrany allel (czyli id bedzie inne niz -1) albo nie ..
		# Musze przekazac do koljenej rundy nazwe contigu i jego kierunkwowc
		# czyli czy kolejny allel ma byc PRZED czy ZA tym allelem jesli bedzie na tym samym contigu

		# jesli to nowy allel to musze okreslic kierunkwowsc na ty mcontigu

		if [ ${CONTIG_DIR} == "NOVEL" ] && [ ${FINAL_ID} -gt 0 ]; then
		# ten allel zaczyna contig musze zdecyfdowac o jego kierunkowosci
			CONTIG_DIR="${FINAL_DIR}"
		fi

		# jesli wybralem allel zapisuje jego pozycje
		if [ ${FINAL_ID} -gt 0 ]; then
			PREVIOUS_CONTIG="${FINAL_CONTIG}"
			if [ ${CONTIG_DIR} == "plus" ]; then
			# znalazlem allel zapisuje jego pozycje 
			PREVIOUS_POS="${FINAL_HIT_END}"

			elif [ ${CONTIG_DIR} == "minus" ]; then 
			PREVIOUS_POS="${FINAL_HIT_START}"
			fi
		fi

		echo -e "${K}\t${FINAL_ID}\t100\tmulti\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}" >> log.log

	elif  [ ${HITY} == 0 ]; then
		# Nie ma allelu o 100% pident z allelami w cgMLST
		# Iterujemy po calym pliku i szukamy takiego allelu o najnizszej sumie dlugosc mapowania - (mismatch + gaps)
		# Defaultowe wartosci gdyby analizowany plik byl pusty
		# Na koniec zwrocimy i tak -1, ale zostawiam sobie do testow dodatkowe informacje
		FINAL_ID=`echo -1`
                FINAL_TRUE_LEN=`echo -1`
                FINAL_QUERY_LEN=`echo -1`
		FINAL_PIDENT="-1"	
		FINAL_SCORE="-1"
		while read LINE; do
		
			ALLEL_ID=`echo ${LINE} | cut -d " " -f1 |  sed -r  s'/.*_(.+)/\1/' `
			ALLEL_LEN=`echo ${LINE} | cut -d " " -f2`
                	QUERY_LEN=`echo ${LINE} | awk '{print $4 - $3 + 1}'`
			NO_MISMATCH=`echo ${LINE} | awk '{print $9}'`
			NO_GAPS=`echo ${LINE} | awk '{print $10}'`
			PIDENT=`echo ${LINE} | awk '{print $7}'`
			SCORE=`echo "${QUERY_LEN} - ${NO_MISMATCH} - ${NO_GAPS}" | bc -l | awk '{print int($1)}'`
			if [ ${SCORE} -gt ${FINAL_SCORE} ]; then
				FINAL_ID=`echo ${ALLEL_ID}`
				FINAL_TRUE_LEN=`echo ${ALLEL_LEN}`
				FINAL_QUERY_LEN=`echo ${QUERY_LEN}`
				FINAL_SCORE=`echo ${SCORE}`
				FINAL_PIDENT=`echo ${PIDENT}`

			fi


		done < blastn_tmp.tab 
		
		# echo -e "${K}\t${FINAL_ID}\t${FINAL_PIDENT}\tunk\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}\t${FINAL_SCORE}" >> log.log
		# docelowo entero i tak daje tu brak allelu wiec aby byc spojnym damy -1
		# mozna dac kod jakis inny ...
		echo -e "${K}\t-1\t${FINAL_PIDENT}\tunk\t${FINAL_TRUE_LEN}\t${FINAL_QUERY_LEN}\t${FINAL_SCORE}" >> log.log
	fi

	rm blastn_tmp.tab
done
