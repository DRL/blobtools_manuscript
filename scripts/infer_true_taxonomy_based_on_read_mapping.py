#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: generate_table_based_on_read_counts_by_sequence.py         -i <FILE> [-h|--help]

    Options:
        -h --help                                                   show this
        -i, --infile <FILE>                                         Input file '*.read_count_by_reference.txt'

"""
from __future__ import division
import sys
from docopt import docopt
import os

def read_file(infile):
    if not infile or not os.path.exists(infile):
        sys.exit("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")

class Data():
    def __init__(self, infile):
        self.infile = infile
        self.reads_by_organims_by_contig = {}
        self.set_of_organism_ids = set()
        self.parse_infile()
        self.write_table()

    def parse_infile(self):
        for line in read_file(self.infile):
            col = line.lstrip(" ").split()
            read_count = int(col[0])
            contig_id = col[1]
            organism_id = col[2]
            self.set_of_organism_ids.add(organism_id)
            if contig_id not in self.reads_by_organims_by_contig:
                self.reads_by_organims_by_contig[contig_id] = {}
            self.reads_by_organims_by_contig[contig_id][organism_id] = read_count

    def write_table(self):
        organism_ids = sorted(self.set_of_organism_ids)
        output = []
        header = ['#contig_id']
        for organism_id in organism_ids:
            header.append(organism_id)
        header.append('taxonomy')
        header.append('fraction_reads')
        output.append("\t".join(header))
        for contig_id in sorted(self.reads_by_organims_by_contig):
            line = []
            line.append(contig_id)
            read_counts = []
            for organism_id in organism_ids:
                read_count_for_organism = self.reads_by_organims_by_contig[contig_id].get(organism_id, 0)
                read_counts.append(read_count_for_organism)
                line.append(str(read_count_for_organism))
            taxonomy = max(self.reads_by_organims_by_contig[contig_id], key=self.reads_by_organims_by_contig[contig_id].get)
            line.append(taxonomy)
            fraction_reads = self.reads_by_organims_by_contig[contig_id][taxonomy]/sum(read_counts)
            line.append("{0:.2f}".format(fraction_reads))
            output.append("\t".join(line))
        print "\n".join(output)





if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    infile = args['--infile']
    data = Data(infile)
