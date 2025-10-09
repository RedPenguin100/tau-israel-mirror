import os
import json
import pandas as pd
import random
import string

from .download_assets import _require_samtools, _require_bowtie, _require_seqkit
from .process_utils import LocusInfoOld
from .utils import create_gene_to_data
from .run_pipe import get_n_best_res
from .consts_dataframe import *
from .hybridization.hybridization_features import get_exp_psrna_hybridization
from .off_target.aso_counts import run_aso_counts

features_for_output = [SEQUENCE,
                       'mod_pattern',
                       'exp_ps_hybr',
                       'gc_content',
                       'at_skew',
                       'sense_start',
                       'on_target_fold_openness_normalized40_15',
                       'sense_avg_accessibility'
                       ]


def add_features_for_output(df):
    df.loc[:, 'exp_ps_hybr'] = [
        get_exp_psrna_hybridization(antisense.replace('T', 'U'), temp=37) for
        antisense in df[SEQUENCE]]
    df.loc[:, 'gc_content'] = df['sense'].apply(lambda seq: (seq.count('G') + seq.count('C')) / len(seq))


def df_to_fasta(df):
    os.makedirs("/tmp", exist_ok=True)

    # random file name
    rand_name = ''.join(random.choices(string.ascii_lowercase + string.digits, k=8))
    path = f"/tmp/{rand_name}.fasta"

    # write fasta
    with open(path, "w") as f:
        for name, seq in zip(df['seq_name'], df[SEQUENCE]):
            f.write(f">{name}\n{seq}\n")
    return path


def add_target_scores(df, all_target_json_path, on_target_json_path):
    with open(all_target_json_path) as f:
        all_target = json.load(f)
    with open(on_target_json_path) as f:
        on_target = json.load(f)

    # Compute off_target and on_target
    def calc_scores(seq):
        if seq in all_target and seq in on_target:
            off = [a - b for a, b in zip(all_target[seq], on_target[seq])]
            on = on_target[seq]
            return off, on
        return [], []  # if seq not found

    scores = df[SEQUENCE].map(calc_scores)
    try:
        print(scores.head())
    except Exception:
        print(list(scores)[:5])
    df.loc[:, 'off_target'] = [s[0] for s in scores]
    df.loc[:, 'on_target'] = [s[1] for s in scores]


def add_off_target(df, full_mRNA_fasta_file):
    aso_fasta_path = df_to_fasta(df)
    all_target_json = run_aso_counts(aso_fasta_path)
    on_target_json = run_aso_counts(aso_fasta_path, target_file=full_mRNA_fasta_file)
    add_target_scores(df, all_target_json, on_target_json)


def design_asos(organismName, geneName, geneData, top_k, includeFeatureBreakdown, return_seq=False):
    """
    Main ASO generation function.

    Args:
        organismName: Only 'human' or 'other' is supported.
        geneName: Selected gene identifier
        geneData: Nucleotide sequence or identifier supplied by the client (optional)
        top_k: Number of top results to return
        includeFeatureBreakdown: Whether to include detailed feature breakdown
        return_seq: Whether to return the mRNA sequence of the target

    Returns:
        dataframe with sequence, chemical modification, and optionally additional features ordered by predicted potency.
    """
    _require_seqkit()
    _require_samtools()
    _require_bowtie()

    session_id = random.randint(1, 1_000_000)
    only_exons = True if geneData else False
    gene_lst = [geneData] if geneData else [geneName]
    off_target_flag = True if organismName == 'human' else False
    if off_target_flag:
        if geneData:
            full_mRNA_fasta_path = f'/tmp/{session_id}/full_mrna_{session_id}.fa'  # empty will be created
        else:
            full_mRNA_fasta_path = f'/tmp/{session_id}/{geneName}.fa'  # empty will be created
    else:
        full_mRNA_fasta_path = None

    if only_exons:
        gene_to_data = dict()
        gene_to_data['one_exon'] = LocusInfoOld(gene_lst[0])
        gene_lst = ['one_exon']
    else:
        print("Creating gene to data: ")
        gene_to_data = create_gene_to_data(gene_lst)

    # MOE_DF
    moe_pre_df, full_mrna = get_n_best_res(gene_lst, top_k, 'moe',
                                           full_mRNA_fasta_file=full_mRNA_fasta_path
                                           , gene_to_data=gene_to_data
                                           )
    moe_df = moe_pre_df[gene_lst[0]].copy()

    # LNA_DF
    lna_pre_df, _ = get_n_best_res(gene_lst, top_k, 'lna', gene_to_data=gene_to_data)
    lna_df = lna_pre_df[gene_lst[0]].copy()

    # merge_df
    df = pd.concat([moe_df, lna_df], ignore_index=True)

    if includeFeatureBreakdown:
        add_features_for_output(df)
        if off_target_flag:
            add_off_target(df, full_mRNA_fasta_path)
            features_for_output.append('off_target')
            features_for_output.append('on_target')
        df = df[features_for_output].copy()
    else:
        df = df[[SEQUENCE, 'mod_pattern']].copy()

    if not return_seq:
        return df
    else:
        return [df, full_mrna]


if __name__ == "__main__":
    # show all columns
    pd.set_option('display.max_columns', None)

    # show full content in each cell (no truncation)
    pd.set_option('display.max_colwidth', None)

    # show all rows (optional)
    pd.set_option('display.max_rows', None)

    print("[__main__] starting demo call")
    df = design_asos(organismName='human', geneName='MALAT1', geneData=None, top_k=3, includeFeatureBreakdown=True)
    print("[__main__] result df shape:", getattr(df, "shape", None))
    print(df)
    print("[__main__] done")
