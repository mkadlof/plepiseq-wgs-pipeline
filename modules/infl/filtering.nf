process filtering {
    tag "filtering:${sampleId}"
    container  = params.main_image
    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_with_index), val(QC_status), path(primers), path(pairs)
    output:
    tuple val(sampleId), path('to_clip_sorted.bam'), path('to_clip_sorted.bam.bai'), path(primers), path(pairs), val(QC_status), emit: one_amplicon_primers_and_QC
    
    script:
    def final_index = -1
    ref_genome_with_index.eachWithIndex { filename, index ->
        if (filename.toString().endsWith(".fasta")) {
         final_index = index
        }
    }
    """
    # Custom filtering for influenza
    # The script takes in sequence:
    # - BAM for filtering and downsampling
    # - Primer scheme
    # - Target coverage per segment (integer)
    # - Minimum read mapping quality (integer)
    # - Reference genome sequences in FASTA (all segments)

    if [ ${QC_status} == "nie" ]; then
      touch to_clip_sorted.bam
      touch to_clip_sorted.bam.bai
    else
      simple_filter_illumina_INFL.py ${bam} ${primers} ${params.max_depth} ${params.min_mapq} ${params.length} ${ref_genome_with_index[final_index]}
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
      simple_filter_nanopore_INFL_ekstralayer_EQA2024.py  ${bam} primers.bed ${params.bed_offset} ${params.max_depth} ${params.min_mapq}  ${params.length}

${prefix}_sorted.bam ${primers} ${bed_offset} ${max_depth} ${mapq} ${length_min} ${length_max} ${alignment_length_min}  ${input_genome} ${window_size}

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
