# -*- coding: utf-8 -*-
import pytest
from pydna.crispr import pam


def test_pam_initialization():
    p = pam("NGG")
    assert p.pam == "NGG"
    assert str(p) == "[ATGC]GG"


def test_pam_to_re():
    p = pam("NNGRRT")
    assert str(p) == "[ATGC]{2}G[AG]{2}T"


def test_pam_rc():
    p = pam("NGG")
    assert p.rcre == "[ATGC]CC"


def test_pam_call():
    p = pam("NNGRRT")
    assert p() == "[ATGC]{2}G[AG]{2}T"


if __name__ == "__main__":
    pytest.main()
