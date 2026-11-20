# readfq taken from https://github.com/lh3/readfq/blob/master/readfq.py


class Seq:

    @staticmethod
    def readfq(fp):  # this is a generator function
        last = None  # this is a buffer keeping the last unprocessed line
        while True:  # mimic closure; is it a bad idea?
            if not last:  # the first record or a record following a fastq
                for line in fp:  # search for the start of the next record
                    if line[0] in '>@':  # fasta/q header line
                        last = line[:-1]  # save this line
                        break
            if not last:
                break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for line in fp:  # read the sequence
                if line[0] in '@+>':
                    last = line[:-1]
                    break
                seqs.append(line[:-1])
            if not last or last[0] != '+':  # this is a fasta record
                yield Read(name, ''.join(seqs), None) # yield a fasta record
                if not last:
                    break
            else:  # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for line in fp:  # read the quality
                    seqs.append(line[:-1])
                    leng += len(line) - 1
                    if leng >= len(seq):  # have read enough quality
                        last = None
                        yield Read(name, seq, ''.join(seqs))  # yield a fastq record
                        break
                if last:  # reach EOF before reading enough quality
                    yield Read(name, seq, None) # yield a fasta record instead
                    break


class Read:

    def __init__(self, name, seq, qual):
        self.name = name
        self.seq = seq
        self.qual = qual
        self.qual_values = []
        if qual:
            self.qual_values = [ord(q) - 33 for q in qual]

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return f"{self.name}\n{self.seq}"
