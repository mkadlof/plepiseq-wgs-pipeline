# Model building: Module modeller

```Bash
    modpy.sh modeller_create_alignment.py ${target_fasta}
    modpy.sh modeller_build_model.py alignment.pir
```

Modeller is a program used for building protein models using homology modeling method (i.e., based on the known structure of a related protein). In our case, we chose protein [7dwz](https://www.rcsb.org/structure/7DWZ) as the basis for modeling. Modeller is software that requires a license key, which needs to be provided during the container build process. As a result, it returns a PDB file containing the trimer of the S protein, as well as a file with sequence alignment. Masked regions are modeled as amino acids labeled with the symbol `UNK` without side chains (effectively as glycines).

Dokumentację formatu PDB można znaleźć pod adresem: [](https://www.wwpdb.org/documentation/file-format)