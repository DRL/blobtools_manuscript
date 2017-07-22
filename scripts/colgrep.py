#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
usage: colgrep.py                   -f <FILE> -i <FILE> -c <INT> [-v]
                                    [-h|--help]

    Options:
        -h --help                   show this
        -f, --file FILE             File in which to search for IDs
        -c, --col INT               Column (1-based) of file in which to search for IDs
        -i, --ids FILE              File containing IDs
        -v, --invert                Invert selection

"""
from __future__ import division
from __future__ import print_function
import sys
from docopt import docopt
import os
from collections import defaultdict

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)
    sys.exit()

def read_file(infile):
    if not infile or not os.path.exists(infile):
        eprint("[X] - File '%s' does not exist." % (infile))
    with open(infile) as fh:
        for line in fh:
            yield line.rstrip("\n")

class Data():
    def __init__(self, infile, column, id_f, invert_flag):
        self.stats = {'total': 0, 'excluded': 0, 'included': 0}
        self.infile = infile
        self.column = column - 1
        self.id_f = id_f
        self.invert_flag = invert_flag
        self.ids = self.parse_ids_f()
        self.filter_infile()

    def parse_ids_f(self):
        ids = []
        for line in read_file(self.id_f):
            ids.append(line.rstrip("\n"))
        return ids

    def filter_infile(self):
        self.status_of_id = {}
        if self.invert_flag:
            self.status_of_id = defaultdict(lambda: True)
            for _id in self.ids:
                self.status_of_id[_id] = False
        else:
            self.status_of_id = defaultdict(lambda: False)
            for _id in self.ids:
                self.status_of_id[_id] = True
        for line in read_file(self.infile):
            self.stats['total'] += 1
            col = line.lstrip(" ").split()
            try:
                val = col[self.column]
            except IndexError:
                eprint("[X] File contains less than %s columns." % (column))
            if self.status_of_id[val] is True:
                self.stats['included'] += 1
                print("\t".join(col))
            else:
                self.stats['excluded'] += 1
        eprint("[+] Total=%s Included=%s Excluded=%s" % (self.stats['total'], self.stats['included'], self.stats['excluded']))

if __name__ == "__main__":
    __version__ = 0.1
    args = docopt(__doc__)
    infile = args['--file']
    try:
        column = int(args['--col'])
    except ValueError:
        eprint("[X] --col has to be an integer.")
    if column < 1:
        eprint("[X] --col has to be greater than 1.")
    id_f = args['--ids']
    invert_flag = args['--invert']
    data = Data(infile, column, id_f, invert_flag)
