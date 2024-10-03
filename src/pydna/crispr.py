#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2013-2023 by Björn Johansson.  All rights reserved.
# This code is part of the Python-dna distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.
"""Provides the Dseq class for handling double stranded DNA sequences.

Dseq is a subclass of :class:`Bio.Seq.Seq`. The Dseq class
is mostly useful as a part of the :class:`pydna.dseqrecord.Dseqrecord` class
which can hold more meta data.

The Dseq class support the notion of circular and linear DNA topology.
"""

from abc import ABC, abstractmethod
import re
from pydna.utils import rc
from Bio.Data.IUPACData import ambiguous_dna_values
from pydna.dseqrecord import Dseqrecord


class pam:
    """
    >>> pam("NGG")
    '[ATGC]GG'
    >>> pam("NNGRRT")
    '[ATGC]{2}G[AG]{2}T'
    """

    def __init__(self, pam: str) -> None:
        self.pam = pam.upper()
        self.re = pam_to_re(self.pam)
        self.rcre = pam_to_re(rc(self.pam))

    def __str__(self) -> str:
        return self.re

    def __repr__(self) -> str:
        return f"pam({self.pam})"

    def __call__(self) -> str:
        return self.re


class _cas(ABC):
    scaffold = "ND"
    pam = "ND"
    size = 0
    fst5 = 0
    fst3 = 0

    def __init__(self, protospacer: str):
        self.protospacer = protospacer.upper()
        self.compsite = re.compile(
            f"(?=(?P<watson>{self.protospacer}{self.pam}))|(?=(?P<crick>{rc(self.pam)}{rc(self.protospacer)}))",
            re.UNICODE,
        )

    @abstractmethod
    def search(self, dna: Dseqrecord, linear: bool = True):
        """To override in subclass."""
        pass

    def __repr__(self):
        return f"{type(self).__name__}({self.protospacer[:3]}..{self.protospacer[-3:]})"

    @abstractmethod
    def __str__(self):
        """To override in subclass."""
        pass


class cas9(_cas):
    """docstring.

         |----size----------|

         ---protospacer------
                          -fst3
         fst5             |-|
         |--------------|
                             PAM
    5-NNGGAAGAGTAATACACTA-AAANGGNN-3
      ||||||||||||||||||| ||||||||
    3-NNCCTTCTCATTATGTGAT-TTTNCCNN-5
        ||||||||||||||||| |||
      5-GGAAGAGTAATACACTA-AAAg-u-a-a-g-g  Scaffold
        ---gRNA spacer---    u-a
                             u-a
                             u-a
                             u-a
                             a-u
                             g-u-g
                            a    a
                             g-c-a
                             c-g
                             u-a
                             a-u
                            g   a  tetraloop
                             a-a
    """

    scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
    pam = "NGG"
    size = 20
    fst5 = 17
    fst3 = -3
    ovhg = fst5 - (size + fst3)

    def search(self, dna: Dseqrecord, linear: bool = True):
        """docstring."""
        dna = str(dna).upper()
        if linear:
            dna = dna
        else:
            dna = dna + dna[1 : self.size]
        results = []
        for mobj in self.compsite.finditer(dna):
            w, c = mobj.groups()
            if w:
                results.append(mobj.start("watson") + 1 + self.fst5)
            if c:
                results.append(mobj.start("crick") + len(self.pam) + 1 - self.fst3)
        return results

    def __str__(self):
        """docstring."""
        return f">{type(self).__name__} protospacer scaffold\n{self.protospacer} {self.scaffold}"


class SaCas9(_cas):
    """
    Represents a SaCas9 CRISPR-Cas system.

    Attributes:
        scaffold (str): The scaffold sequence for the SaCas9 system.
        pam (str): The PAM sequence taergeted by the SaCas9 system.
        size (int): The length of the protospacer.
        fst5 (int): The position of the first base of the protospacer from the 5' end.
        fst3 (int): The position of the first base of the protospacer from the 3' end.
        ovhg (int): TODO: ?
    """

    scaffold = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGG"
    pam = "NNGRRT"

    size = 21
    fst5 = 17
    fst3 = -3
    ovhg = fst5 - (size + fst3)

    def __init__(self, protospacer: str):
        self.pamre = pam_to_re(self.pam)
        self.pamrerc = pam_to_re(rc(self.pam))
        self.protospacer = protospacer.upper()
        # TODO: Here we do the rc of a regex instead of the PAM
        self.compsite = re.compile(
            f"(?=(?P<watson>{self.protospacer}{self.pamre}))|(?=(?P<crick>{self.pamrerc}{rc(self.protospacer)}))",
            re.UNICODE,
        )

    def search(self, dna: Dseqrecord, linear: bool = True):
        """docstring."""
        dna = str(dna).upper()
        if linear:
            dna = dna
        else:
            dna = dna + dna[1 : self.size]
        results = []
        for mobj in self.compsite.finditer(dna):
            w, c = mobj.groups()
            if w:
                results.append(mobj.start("watson") + 1 + self.fst5)
            if c:
                results.append(mobj.start("crick") + len(self.pam) + 1 - self.fst3)
        return results

    def __str__(self):
        """docstring."""
        return f">{type(self).__name__} protospacer scaffold\n{self.protospacer} {self.scaffold}"


def protospacer(guide_construct: Dseqrecord, cas=cas9):
    """docstring."""
    in_watson = [
        mobj.group("ps")
        for mobj in re.finditer(f"(?P<ps>.{{{cas.size}}})(?:{cas.scaffold})", str(guide_construct.seq).upper())
    ]
    in_crick = [
        rc(mobj.group("ps"))
        for mobj in re.finditer(f"(?:{rc(cas.scaffold)})(?P<ps>.{{{cas.size}}})", str(guide_construct.seq).upper())
    ]
    return in_watson + in_crick


def pam_to_re(pam: str) -> str:
    """
    Convert a PAM sequence to a regular expression.
    """
    pam_re = ""
    count = 1
    last_nuc = None
    for c, nuc in enumerate(pam):
        new_nuc = ambiguous_dna_values[nuc]
        if new_nuc == last_nuc:
            count += 1
            continue
        else:
            # Check if we need to add a count to the last nuc
            if count > 1:
                pam_re += f"{{{count}}}"

            # Update last nuc
            last_nuc = new_nuc

            # Add to regexp
            if len(new_nuc) == 1:
                pam_re += new_nuc
            else:
                pam_re += f"[{new_nuc}]"

            # Reset count
            count = 1
    else:
        # Add final count if needed
        if count > 1:
            pam_re += f"{{{count}}}"
    return pam_re


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
