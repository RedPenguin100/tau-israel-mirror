import pytest
from asodesigner.read_human_genome import get_locus_to_data_dict


def test_sanity():
    gene_to_data = get_locus_to_data_dict(gene_subset=['DDX11L1'])
    gene_to_data['DDX11L1']
