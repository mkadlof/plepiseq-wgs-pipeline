process filtering {
    tag "filtering:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path(bam), path(bai), val(QC_status), path(primers), path(pairs)

    output:
    tuple val(sampleId), path('first_pass_sorted.bam'), path('first_pass_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC
    tuple val(sampleId), path('two_amplicons_sorted.bam'), emit: two_amplicon_only
    // tuple val(sampleId), path('Statystyki.txt')

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch first_pass_sorted.bam
      touch first_pass_sorted.bam.bai
      touch two_amplicons_sorted.bam
      QC_exit="nie"
    else
      bed_offset=0
      length=`echo "${params.length} - 20" | bc -l`
      equal_depth=`echo "${params.max_depth} / 2" | bc -l | awk '{print int(\$0)}'` 
      simple_filter_illumina_one_segment.py ${bam} ${primers} \${bed_offset} \${length} ${params.min_mapq} ${params.window_size} \${equal_depth}
      samtools index first_pass_sorted.bam
    
      if [ -e two_amplicons_sorted.bam ]; then
          cp two_amplicons_sorted.bam tmp.bam
      else
          samtools view -b -o two_amplicons_sorted.bam -L /dev/null first_pass_sorted.bam
      fi
      # Komentarz do json-a - Primer_usage.txt zawiera plik z 3 kolumanmi - nazwa segmentu, nazwa primera, uzyciem
    fi
    """
}

process filtering_nanopore {
    // Nanopore filtering and subsequent primer masking does not require ivar specific pairs file
    tag "filtering:${sampleId}"
    container  = params.main_image
    input:
    // emit to_coinfection z minimap2 
    tuple val(sampleId), path(bam), path(bai),  path(ref_genome_with_index), path(primers), val(QC_status)

    output:
    tuple val(sampleId), path('to_classical_masking.bam'), path('to_overshot_masking.bam'), path(ref_genome_with_index), path(primers), val(QC_status), emit: to_masking

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch reads_inner_strict.bam
      touch first_pass_sorted.bam.bai
      touch two_amplicons_sorted.bam
      QC_exit="nie"
    else
      masking_cap=`echo "${params.mask} + 10" | bc -l` #  do jakiej maksymalnej warotosci podbijac coverage w regionach w ktorych brakowalo oczekiwanych jedno-amplikonowych odczytow
      # jako ze korzystamy z odczytow o niejasnym pochodzeniu (najczesciej odczty z fuzji amplikonow)
      # to dobijamy tylk odo wartosci tak aby podbic zeby region nie byl maskowany
      simple_filter_nanopore_final_with_windowstep.py ${bam} primers.bed ${params.bed_offset} ${params.max_depth} ${params.length} ${params.min_mapq} ${params.extra_bed_offset} \${masking_cap}

      # Skrypt wyzej zwraca bardzo duzo plikow, niestety aktualnie ich powstawanie jest zalezne od danych (zawsze bedzie reads_inner_strict.bam, reszta jest opcjonalna)
      # POPRAWIC TO W KOD REVIEW
      # reads_inner_strict.bam odczyty mapujace sie na jeden amplikon (najlepsze jakosciowo)
      # reads_two_amplicons*bam - lista odczytow ktore powstaly na skutek wypadniecia danego amplikonu (a wiec lacza amplikony sasiednie z POOL-a primerow, biologicznie mozliwe)
      # reads_overshot.bam - odczyty mapujace sie na raczej jeden primer ale ktore "przestrzlily" swoj amplikon ( ready nieprwadopodobnie biologicznie, ale w skrajnych przypadkach uzywane)

      # Odczyty pochodzace z dwoch amplikonow sa w oddzielnych plikach nalezy je polaczyc w jeden plik razem z plikiem reads_inner_strict.bam zawierajacym "dobre" odczyty
      # Te odczyty beda filtrowane z wykorzystaniem informacji o granicach primerow z pliku bed
      TO_MERGE_OVER=`ls -l  reads_two_amplicons_*.bam | tr -s " " |  cut -d " " -f9 | tr "\n" " "`
      if [ `ls -l  reads_two_amplicons*bam | wc -l` -gt 0 ];then
        samtools merge -o to_classical_masking.bam \${TO_MERGE_OVER} reads_inner_strict.bam
        rm reads_two_amplicons_*
        rm reads_inner_strict.bam
      else
        mv reads_inner_strict.bam to_classical_masking.bam
        
      fi

      # odczyty reads_overshot.bam nie sa wlaczane bo beda filtrowane tak by uciac ich "nadmiarowa" sekwencje tak by byly przyciete do primerow
      # na potrzeby nextflow musimy stworzyc ten plik jesli go nie ma
      if [ ! -e reads_overshot.bam ]; then
        touch to_overshot_masking.bam
      else
        mv reads_overshot.bam to_overshot_masking.bam
      fi

    fi
    """
}
