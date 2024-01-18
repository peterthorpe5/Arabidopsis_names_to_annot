#!/usr/bin/env python
# Title: plot upset of clusters
# author: Peter Thorpe September 2021. UofSA , UK.

from upsetplot import from_memberships
import upsetplot
import pandas
import matplotlib
import seaborn
from upsetplot import plot
from matplotlib import pyplot
from optparse import OptionParser
import os
from sys import stdin,argv
import sys


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
$ python namne_of_script -h
"""

parser = OptionParser(usage=usage)
parser.add_option("-i", dest="in_file", default=None,
                  help="the orthofinder cluster file. txt",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="outfile", default="default.pdf",
                  help="output filenames")

(options, args) = parser.parse_args()
in_file = options.in_file
outfile = options.outfile

##################################################################


def parse_clusters(in_file):
    """ func to parse orthofinder cluster and return
    from_membership class for plotting"""
    # from_mbership take a list of lists. 
    overall_list = []
    with open(in_file, "r") as fh:
        for line in fh:
            tmp_rename = []
            genes = line.split()
            for gene in genes:
                if gene.startswith("OG"):
                    continue
                prefix =  gene.split("_")[0]
                tmp_rename.append(prefix)
            overall_list.append(tmp_rename)
    example = from_memberships(overall_list)
    return example

def plot_the_data(example, outfile):
    """func to plot the data. Takes in the
    from_membership"""
    upsetplot.plot(example, subset_size="sum", show_counts=False)
    # pyplot.show()

    pyplot.savefig(outfile)

    

# Run as script
if __name__ == '__main__':
    if not os.path.isfile(in_file):
        print("sorry cannot find you %s file" % in_file)
        os._exit(0)
    print("looking at: %s" % in_file)
    example = parse_clusters(in_file)
    if outfile == "default.pdf":
        outfile = in_file + "_upset.pdf"
    plot_the_data(example, outfile)
        
