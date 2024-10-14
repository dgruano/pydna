# -*- coding: utf-8 -*-
import re


class USER:
    pattern = "([ACGT]{5}U)"
    size = 6
    fst5 = 7  # First 5' cut
    fst3 = None  # (there are no cuts in complementary strand)
    site = "NNNNNU"
    ovhg = fst5 - 1

    def __init__(self):
        self.compsite = re.compile("(?=(?P<USER>[ACGT]{5}U))", re.UNICODE)

    def search(self, dna, linear=True):
        """
        Search function for USER enzyme.

        Parameters
        ----------
        dna : Dseq
            Dseq object representing the DNA sequence to search for USER site.
        linear : bool
            If True, the search is performed on the input sequence.
            If False, the search is performed on the sequence + sequence[1:].

        Returns
        -------
        list
            A list of the positions of the USER target sites.
        """
        # TODO: Deal with circular DNA
        results = []
        for mobj in re.finditer("[ACGT]{5}U", dna.watson):
            print(mobj.start(), mobj.group())
            results.append(mobj.start() + self.fst5)

        for mobj in re.finditer("U[ACGT]{5}", dna.crick):
            print(mobj.end(), mobj.group())
            results.append(mobj.end() + self.fst5 - 1)

        return results

    def __repr__(self):
        return f"USER({self.site})"

    def __str__(self):
        return f"USER({self.site})"


if __name__ == "__main__":
    import os as _os

    cached = _os.getenv("pydna_cached_funcs", "")
    _os.environ["pydna_cached_funcs"] = ""
    import doctest

    doctest.testmod(verbose=True, optionflags=doctest.ELLIPSIS)
    _os.environ["pydna_cached_funcs"] = cached
