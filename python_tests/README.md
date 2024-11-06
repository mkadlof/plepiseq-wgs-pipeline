Before running tests you need to run pipeline to generate results. Following structure is assumed:

```
/nf_illumina_sars
├── git_repo
│     ├── (...)
│     ├── plepiseq_json
│     ├── python_tests
│     └── (...)
├── gold_files
│     ├── eqa2023-sars-two-samples
│     │     ├── ESIB_EQA_2023.SARS2.01
│     │     └── ESIB_EQA_2023.SARS2.02
│     └── eqa2024-infl-two-samples
│         ├── INFL2.05
│         └── INFL2.10
├── run_pipeline
├── testing_data
│     ├── EQA2023
│     ├── EQA2024
│     │     └── infl
│     └── RSV
└── work_dirs
    ├── eqa2023-sars-two-samples
    │     ├── ESIB_EQA_2023.SARS2.01
    │     ├── ESIB_EQA_2023.SARS2.02
    │     ├── reports
    │     └── work
    ├── eqa2024-infl-two-samples
    │     ├── INFL2.05
    │     ├── INFL2.10
    │     ├── reports
    │     └── work
    └── rsv-two-samples
        ├── reports
        ├── SRR27383064
        ├── SRR27383065
        └── work
```

To run tests, run the following command:

```bash
PYTHONPATH=./git_repo pytest --tb=no -n 4 ./git_repo/python_tests/test_all.py
```
