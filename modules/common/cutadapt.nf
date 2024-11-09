process cutadapt {
// Na razie nie uzywany w grypie uzywam cutadapta do wywalania koncowek slbych odczytow, alw w SARS i RSV juz tego nie robie
// zostawiam fragment kodu a samm odul mozna uznac za analog trimmomatica z illuminy
cutadapt -j ${cpu} -q ${quality_coverage},${quality_coverage} --trim-n -o tmp.fastq.gz  ${reads}

}
