import pytest
from asodesigner.util import get_nucleotide_watson_crick

from asodesigner.util import get_antisense


def test_get_nucleotide_watson_crick():
    assert get_nucleotide_watson_crick('A') == 'T'
    assert get_nucleotide_watson_crick('C') == 'G'
    assert get_nucleotide_watson_crick('G') == 'C'
    assert get_nucleotide_watson_crick('T') == 'A'
    assert get_nucleotide_watson_crick('U') == 'A'
    with pytest.raises(ValueError):
        get_nucleotide_watson_crick('ZZ')
        get_nucleotide_watson_crick('Z')


def test_get_antisense():
    assert get_antisense('AGCT') == 'AGCT'
    assert get_antisense('TTTT') == 'AAAA'