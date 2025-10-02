import pytest
import numpy as np

from asodesigner.features.vienna_fold import get_weighted_energy, calculate_energies
from tests.conftest import TEST_CACHE_PATH


def test_regression(short_mrna):
    energies = calculate_energies(str(short_mrna.full_mrna), 15, 40)
    print(energies[-3:])
    energy = get_weighted_energy(2525, 16, 15, energies, 40)

    assert pytest.approx(energy, rel=1e-2) == -3.56


def test_sanity(short_mrna):
    gene_length = 2541

    window_size = 40
    step_size = 15
    sample_energies = np.zeros((gene_length - window_size) // step_size + 1 + 1)
    sample_energies[-1:] = -1

    energy = get_weighted_energy(2525, 16, 15, sample_energies, 40)
    # Last window ends in 2540, unique coverage for 2530-2541
    # Before last window is 0, but joint with last.
    # So in total the fold should be the weighted average of the two windows
    assert pytest.approx(energy, rel=1e-2) == (5 * ((0 - 1) / 2) + 11 * (-1)) / 16

    energy = get_weighted_energy(2526, 16, 15, sample_energies, 40)
    assert pytest.approx(energy, rel=1e-2) == (4 * ((0 - 1) / 2) + 12 * (-1)) / 16

    sample_energies[-2:-1] = -2
    energy = get_weighted_energy(2525, 16, 15, sample_energies, 40)
    assert pytest.approx(energy, rel=1e-2) == (5 * ((-2 - 1) / 2) + 11 * (-1)) / 16

    sample_energies[-2:] = -2
    energy = get_weighted_energy(2525, 16, 15, sample_energies, 40)
    assert pytest.approx(energy, rel=1e-2) == -2
