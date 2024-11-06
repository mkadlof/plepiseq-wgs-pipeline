import os
from typing import Tuple, List

import pytest
from Bio import SeqIO


def get_test_cases() -> List[Tuple[str, str]]:
    """
    Pobiera testy z pliku tekstowego i zwraca listę krotek zawierających nazwę testu i próbki.

    Funkcja odczytuje plik 'python_tests/samples_and_tests.txt', gdzie każda linia zawiera
    nazwę testu i nazwę próbki oddzielone tabem. Następnie sprawdza, czy istnieje katalog
    odpowiadający nazwie testu w folderze 'work_dirs'. Jeśli tak, dodaje krotkę (nazwa_testu, próbka)
    do listy wynikowej.

    Returns:
        List[Tuple[str, str]]: Lista krotek z nazwami testów i próbkami, gdzie katalogi istnieją.
    """
    test_cases = []
    with open('git_repo/python_tests/samples_and_tests.txt') as f:
        for line in f.readlines():
            test_name, sample = line.split()
            if os.path.isdir(f'work_dirs/{test_name}/'):
                test_cases.append((test_name, sample))
    return test_cases


def get_sequence(path: str) -> str:
    return str(SeqIO.read(path, 'fasta').seq)


# Test cases
@pytest.mark.parametrize('test_name, sample', get_test_cases())
def test_sample_directory_exists(test_name, sample):
    path = f'work_dirs/{test_name}/{sample}'
    assert os.path.isdir(path), f"Directory {path} does not exists."


@pytest.mark.parametrize('test_name, sample', get_test_cases())
def test_fasta_files_exists(test_name, sample):
    if 'sars' in sample.lower():
        file_path = f'work_dirs/{test_name}/{sample}/output_consensus_masked_SV.fa'
        assert os.path.isfile(file_path), f"Fasta file {file_path} does not exists."
    elif 'infl' in sample.lower():
        segments = ['chr1_PB2', 'chr2_PB1', 'chr3_PA', 'chr4_HA', 'chr5_NP', 'chr6_NA', 'chr7_MP', 'chr8_NS']
        file_list = [f'consensus_{segment}.fasta' for segment in segments]
        file_paths = [f'work_dirs/{test_name}/{sample}/{file}' for file in file_list]
        for file_path in file_paths:
            assert os.path.isfile(file_path), f"Fasta file {file_path} does not exists."
    else:
        pytest.skip(f'Test for {sample} are not implemented.')


@pytest.mark.parametrize('test_name, sample', get_test_cases())
def test_fasta_file_is_not_empty(test_name, sample):
    if 'sars' in sample.lower():
        file_path = f'work_dirs/{test_name}/{sample}/output_consensus_masked_SV.fa'
        assert os.path.getsize(file_path) > 0, f'Fasta file {file_path} is empty.'
    elif 'infl' in sample.lower():
        segments = ['chr1_PB2', 'chr2_PB1', 'chr3_PA', 'chr4_HA', 'chr5_NP', 'chr6_NA', 'chr7_MP', 'chr8_NS']
        file_list = [f'consensus_{segment}.fasta' for segment in segments]
        file_paths = [f'work_dirs/{test_name}/{sample}/{file}' for file in file_list]
        for file_path in file_paths:
            assert os.path.getsize(file_path) > 0, f'Fasta file {file_path} is empty.'
    else:
        pytest.skip(f'Test for {sample} are not implemented.')


@pytest.mark.parametrize('test_name, sample', get_test_cases())
def test_DNA_sequence_match(test_name, sample):
    if 'sars' in sample.lower():
        expected_sequence = get_sequence(f"gold_files/{test_name}/{sample}/output_consensus_masked_SV.fa")
        actual_sequence = get_sequence(f"work_dirs/{test_name}/{sample}/output_consensus_masked_SV.fa")
        assert expected_sequence == actual_sequence, f"Sequences do not match for {test_name}/{sample}."
    elif 'infl' in sample.lower():
        segments = ['chr1_PB2', 'chr2_PB1', 'chr3_PA', 'chr4_HA', 'chr5_NP', 'chr6_NA', 'chr7_MP', 'chr8_NS']
        for segment in segments:
            expected_sequence = get_sequence(f"gold_files/{test_name}/{sample}/output_{segment}.fasta")
            actual_sequence = get_sequence(f"work_dirs/{test_name}/{sample}/consensus_{segment}.fasta")
            assert expected_sequence == actual_sequence, f"Sequences do not match for {test_name}/{sample} Segment: {segment}."
    else:
        pytest.skip(f'Test for {sample} are not implemented.')


@pytest.mark.parametrize('test_name, sample', get_test_cases())
def test_output_json_exists(test_name, sample):
    file_path = f'work_dirs/{test_name}/{sample}/output.json'
    assert os.path.isfile(file_path), f"File {file_path} does not exists."


@pytest.mark.parametrize('test_name, sample', get_test_cases())
def test_output_json_is_valid(test_name, sample):
    from plepiseq_json.scheme_validation import validate_json
    schema_file = 'git_repo/plepiseq_json/main_output_schema.json'
    json_file = f'work_dirs/{test_name}/{sample}/output.json'
    validate_json(schema_file, json_file)  # Raises exception if invalid
