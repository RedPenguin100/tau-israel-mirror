import json

import pytest

import numpy as np
import pandas as pd

from asodesigner.consts import DATA_PATH
from asodesigner.features.cai import calc_CAI_weight
from asodesigner.populate.populate_cai import populate_cai_for_aso_dataframe, load_mrna_by_gene_from_files
from tests.test_sense_accessibility import get_init_df


SUPPORTED_CELL_LINES = ['A431']

def test_cache_weights(override=False):
    weights_flat_dict = {}
    for cell_line in SUPPORTED_CELL_LINES:
        if cell_line == 'A431':
            MRNA_FILENAME = DATA_PATH / 'human' / 'transcripts' / 'ACH-001328_transcriptome.csv'
        else:
            raise ValueError(f'Supporting only {SUPPORTED_CELL_LINES} at the moment.')

        transcript_df = pd.read_csv(MRNA_FILENAME)
        transcript_df.loc[:, "ref sequence"] = (
            transcript_df["Mutated Transcript sequence"].fillna(transcript_df["Original Transcript sequence"]))

        TOP_N = 300
        SEQ_COL = "ref sequence"
        EXPR_COL = "expression_norm"

        # Basic checks
        assert EXPR_COL in transcript_df.columns, f"Missing '{EXPR_COL}' column"
        assert SEQ_COL in transcript_df.columns, f"Missing '{SEQ_COL}' column"

        # 1) Pick top-N by expression_norm
        ref_df = transcript_df.sort_values(EXPR_COL, ascending=False).head(TOP_N).copy()

        # 2) Take their sequences as-is (mRNA with U's; calc_CAI_weight handles U->T internally)
        reference_seqs = ref_df[SEQ_COL].dropna().astype(str).tolist()

        # 3) Build CAI weights
        # TODO: cache weights and then delete this logic
        _, weights_flat = calc_CAI_weight(reference_seqs)
        weights_flat_dict[cell_line] = weights_flat

    if override:
        with open(f'weights_cache.json', 'w') as f:
            json.dump(weights_flat_dict, f)



def test_regression(short_mrna):
    aso_df = get_init_df(short_mrna.full_mrna, [16])
    populate_cai_for_aso_dataframe(aso_df, short_mrna)

    np.testing.assert_allclose(
        aso_df['CAI_score_global_CDS'],
        np.ones_like(aso_df['CAI_score_global_CDS']) * 0.85626
    )
