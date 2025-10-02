import numpy as np

from asodesigner.populate.populate_cai import populate_cai_for_aso_dataframe
from tests.test_sense_accessibility import get_init_df

def test_regression(short_mrna):
    aso_df = get_init_df(short_mrna.full_mrna, [16])
    populate_cai_for_aso_dataframe(aso_df, short_mrna)

    np.testing.assert_allclose(
        aso_df['CAI_score_global_CDS'],
        np.ones_like(aso_df['CAI_score_global_CDS']) * 0.85626
    )
