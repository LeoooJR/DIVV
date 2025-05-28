import os
import pytest
from sortedcontainers import SortedSet
import sys

tests_dir = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(tests_dir)
if project_dir not in sys.path:
    sys.path.insert(0, project_dir)

import errors
from variants import VariantRepository

@pytest.fixture
def dummy_filters():

    return {"SNP": False, 
            "TRANSITION": False, 
            "MNP": False, 
            "INDELS": False, 
            "SV": False, 
            "VARS": False, 
            "PASS_ONLY": True}

@pytest.mark.variant
def test_setter_getter(dummy_filters):

    input: list = ['1', '10', 'X', '22', '3', 'Y']

    repository: VariantRepository = VariantRepository()

    repository.chromosomes = input

    assert isinstance(repository.chromosomes, SortedSet)
    assert len(repository.chromosomes) == 6

    repository.filters = dummy_filters

    assert repository.filters == {
                                    "exclude_snps": dummy_filters.get("SNP", False),
                                    "exclude_indels": dummy_filters.get("INDELS", False),
                                    "exclude_vars": dummy_filters.get("VARS", False),
                                    "exclude_mnps": dummy_filters.get("MNP", False),
                                    "exclude_transitions": dummy_filters.get("TRANSITION", False),
                                    "exclude_svs": dummy_filters.get("SV", False),
                                    "pass_only": dummy_filters.get("PASS_ONLY", False),
                                }
    
@pytest.mark.parametrize("unsorted_chromosomes, sorted_chromosomes", [
    (['1', '10', '22', '3'], ['1', '3', '10', '22']),
    (['1', '10', 'X', '22', '3', 'Y'], ['1', '3', '10', '22', 'X', 'Y'])
])

@pytest.mark.variant
def test_chromosome_sorting(unsorted_chromosomes, sorted_chromosomes):

    assert sorted(unsorted_chromosomes, key=VariantRepository.chromosome_sort_key) == sorted_chromosomes

    repository: VariantRepository = VariantRepository()

    repository.chromosomes = unsorted_chromosomes

    assert repository.chromosomes == set(sorted_chromosomes)

@pytest.mark.parametrize("genotype, homozygous", [
    ("", False),
    ("1/1", True),
    ("1|1", True),
    ("1.1", False),
    ("1/0", False),
    ("0|1", False),
])

@pytest.mark.variant
def test_homozygous(genotype, homozygous):

    assert VariantRepository.is_homozygous(genotype) == homozygous

def test_homozygous_error():

    input: int = 2

    with pytest.raises(errors.VariantError):

        VariantRepository.is_homozygous(input)

@pytest.mark.parametrize("genotype, heterozygous", [
    ("", False),
    ("1/1", False),
    ("1|1", False),
    ("1.1", False),
    ("1/0", True),
    ("0|1", True),
])
    
@pytest.mark.variant
def test_heterozygous(genotype, heterozygous):

    assert VariantRepository.is_heterozygous(genotype) == heterozygous

@pytest.mark.variant
def test_heterozygous_error():

    input: int = 2

    with pytest.raises(errors.VariantError):

        VariantRepository.is_heterozygous(input)