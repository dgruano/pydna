# -*- coding: utf-8 -*-
import pytest
from pydna.crispr import pam


def test_pam_initialization():
    p = pam("NGG")
    assert p.pam == "NGG"
    assert str(p) == "[GATC]G{2}"


def test_pam_to_re():
    p = pam("NNGRRT")
    assert str(p) == "[GATC]{2}G[AG]{2}T"


def test_pam_rc():
    p = pam("NGG")
    assert p.rcre == "C{2}[GATC]"


def test_pam_call():
    p = pam("NNGRRT")
    assert p() == "[GATC]{2}G[AG]{2}T"


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
