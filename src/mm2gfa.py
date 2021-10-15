from Bio import SeqIO
import sys


def main(argc, argv):
    if argc != 4:
        sys.stdout.write("mm2gfa.py <stringmatrix.mm> <readNameMap_0> <reads.fa>\n")
        return 1

    _, mm_fname, readmap_fname, reads_fname = argv

    readmap = {}

    for line in open(readmap_fname, "r"):
        lvec = line.rstrip().split('\t')
        readid = int(lvec[0])
        readname = lvec[1].split('>')[1].split()[0]
        readmap[readname] = readid

    seqmap = {}
    for record in SeqIO.parse(reads_fname, "fasta"):
        readid = readmap[record.id]
        seqmap[readid] = (str(record.seq), record.id)

    string_graph = {}

    f = open(mm_fname, "r")

    f.readline()
    f.readline()

    strand = {}

    for line in f.readlines():
        src,dest,dirv,overlap = (int(v) for v in line.rstrip().split())

        if dest in string_graph and src in string_graph[dest]: continue

        if not src in strand:
            strand[src] = '+'
            string_graph[src] = {}

        if not dest in strand:
            if dirv in [1,2]:
                strand[dest] = strand[src]
            else:
                strand[dest] = '+' if strand[src] == '-' else '-'
            string_graph[dest] = {}

        string_graph[src][dest] = (strand[src], strand[dest])

    f.close()

    for src in string_graph:
        seq, name = seqmap[src]
        sys.stdout.write("S\t{}\t{}\n".format(name, seq))
        dests = string_graph[src]
        for dest in dests:
            sys.stdout.write("L\t{}\t{}\t{}\t{}\t*\n".format(name, strand[src], seqmap[dest][1], strand[dest]))


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
