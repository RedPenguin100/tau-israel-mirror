from .utils import create_gene_to_data, create_and_load_model, tbl
from .feature_creation import fill_dfs
from .consts_dataframe import *

import os
from .process_utils import LocusInfoOld


# only_exons : true --> genes_lst[0] = Seq
# only_exons : false --> genes_lst[0] = gene_name


selected_features = [TREATMENT_PERIOD, 'at_skew', 'CAI_score_global_CDS', 'stop_codon_count', 'sense_avg_accessibility',
                     'on_target_fold_openness_normalized40_15', 'sense_utr', 'nucleotide_diversity', 'internal_fold',
                     'normalized_start', 'RNaseH1_Krel_score_R7_krel',  # renamed to best
                     'hairpin_score', 'Modification_min_distance_to_3prime', 'at_rich_region_score']

def seq_to_fasta(seq, path):
    seq = str(seq)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(">sequence\n")
        f.write(seq + "\n")


def get_n_best_res(genes_lst, n, mod_type, only_exons, tp=24, full_mRNA_fasta_file=None):
    assert (len(genes_lst) == 1)
    res = {}
    gene_to_data = dict()
    if only_exons:
        gene_to_data['one_exon'] = LocusInfoOld(genes_lst[0])
        genes_lst = ['one_exon']
    else:
        gene_to_data = create_gene_to_data(genes_lst)
    if full_mRNA_fasta_file:
        seq_to_fasta(gene_to_data[genes_lst[0]].full_mrna, full_mRNA_fasta_file)
    dfs = fill_dfs(genes_lst, gene_to_data, tp=tp, mod_type=mod_type)
    model = create_and_load_model()
    for gene, df in dfs.items():
        df = dfs[gene].copy()
        df["score"] = model.predict(df[selected_features].values)
        df["sense"] = df[SEQUENCE].astype(str).str.translate(tbl).str[::-1]
        top_n = df.nlargest(n, 'score')
        top_n['seq_name'] = [f"{mod_type}_{i}" for i in range(len(top_n))]
        res[gene] = top_n.copy()
    return res


if __name__ == "__main__":
    tp = 24
    n = 12
    genes_lst = ['MALAT1']
    res = get_n_best_res(genes_lst, n, tp)
    print(res['MALAT1'])
