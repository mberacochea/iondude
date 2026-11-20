_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def complement(seq):
    """
    Complement a DNA sequence, ie. ACGT gives TGCA.
    :param seq: DNA sequence
    :return: DNA sequence
    :raises ValueError: if the sequence has a base other than A, C, G or T
    """
    if set(seq) - set("ACGTacgt"):
        raise ValueError(f"not a DNA sequence: {seq}")
    return seq.translate(_COMPLEMENT)


def reverse_complement(seq):
    """
    Reverse complement a DNA sequence, ie. ACGTT gives AACGT.
    :param seq: DNA sequence
    :return: DNA sequence
    :raises ValueError: if the sequence has a base other than A, C, G or T
    """
    return complement(seq)[::-1]
