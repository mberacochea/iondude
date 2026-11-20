def print_t(file, *args):
    """
    Print the args on a file
    :param file: output file
    :param args: str list to print
    """
    print('\t'.join(str(a) for a in args), file=file)


def kmer_for_base(kmers, context, base, k):
    """
    Get the count of a kmer for each possible base
    :param kmers: Dict with the kmeners and the count
    :param context: Current context
    :param base: base to consider
    :param k: position
    :return: the count
    """
    if base != context[k]:
        mod_context = context[:k] + base + context[k + 1:]
    else:
        mod_context = context
    return kmers.get(mod_context, 0)


def get_context(read_seq, pos, k):
    """
    Get the context of size lmer k and rightmer k
    :param read_seq: sequence
    :param pos: position
    :param k: lmer size { right and left } + 1 on the middle
    :return: string with the context
    """
    if pos - k < 0 or pos + k >= len(read_seq):
        return None
    return read_seq[pos - k: pos + k + 1]


def get_kmers(read_seq, k):
    """
    Get the kmers
    :param read_seq: sequence
    :param k: size
    :return: the kmers, none if k is longer than the sequence
    """
    for i in range(len(read_seq) + 1 - k):
        yield read_seq[i: i + k]


def get_reference_context(fasta_file, name, position, k):
    """
    Get a context of length end-start from a fasta file
    :param fasta_file: Fasta file (pysam) - pysam.Fastafile
    :param name: sequence name
    :param position: reference position
    :param k: rmer and lmer size
    :return: context
    """
    return fasta_file.fetch(reference=name, start=position - k, end=position + k + 1)
