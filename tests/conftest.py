import pytest
import os
import pickle

from pathlib import Path

from asodesigner.read_human_genome import get_locus_to_data_dict

TESTS_PATH = Path(os.path.dirname(__file__))
TEST_CACHE_PATH = TESTS_PATH / "cache"


@pytest.fixture
def short_mrna():
    target_gene = 'DDX11L1'

    test_cache = TEST_CACHE_PATH / 'gene_to_data_test.pickle'
    if not test_cache.exists():
        gene_to_data = get_locus_to_data_dict(include_introns=True, gene_subset=[target_gene])
        with open(test_cache, 'wb') as f:
            pickle.dump(gene_to_data, f)
    else:
        with open(test_cache, 'rb') as f:
            gene_to_data = pickle.load(f)

    yield gene_to_data[target_gene]
