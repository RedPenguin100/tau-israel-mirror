import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from asodesigner.download_and_create_human_genome import get_locus_to_data_dict

def get_all_canonical_gene_names():
    locus_to_data = get_locus_to_data_dict()
    return list(locus_to_data.keys())


if __name__ == "__main__":
    # Use get_locus_to_data_dict to get the locus_to_data mapping
    locus_to_data = get_locus_to_data_dict()
    # Get the first gene's data object
    first_gene = next(iter(locus_to_data.values()))
    # Print the attribute/column names of the data object
    print("Column names in the gene data object:")
    print(list(vars(first_gene).keys()))