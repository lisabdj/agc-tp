#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
from operator import itemgetter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "BOUARROUDJ Lisa"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["BOUARROUDJ Lisa"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "BOUARROUDJ Lisa"
__email__ = "lisa.bdj.95@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """Lit un fichier fasta.gz et renvoie les séquences.
      :Paramètres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des séquences
      Retourne:
          Générateur de séquences de longueur minimale
    """
    with gzip.open(amplicon_file, "rt") as filin:
        seq = ""
        for line in filin:
            if line[0] == ">":
                if len(seq) > minseqlen:
                    yield seq
                seq= ""
            elif len(line) > 0:
                seq += line.strip()
        if len(seq) > minseqlen:
            yield seq



def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Identifie les séquences uniques avec un nombre d'occurrence minimal.
      :Paramètres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des séquences
          mincount: nombre minimal d'occurrence
      Retourne:
          Générateur de séquences uniques ayant un nombre d'occurrence minimal
    """
    liste_seq = [] 
    sequences = read_fasta(amplicon_file, minseqlen)
    for seq in sequences:
        liste_seq.append(seq)
    occurrences = [[x, liste_seq.count(x)] for x in set(liste_seq)]
    occurrences = sorted(occurrences, key=itemgetter(1), reverse=True)
    for occ in occurrences:
        if occ[1] >= mincount:
            print(occ)
            #yield occ


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    
    return 

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    #sequence=read_fasta(args.amplicon_file,200)
    #for seq in sequence:
    #    print(seq)
    #read_fasta(args.amplicon_file, 500)
    dereplication_fulllength(args.amplicon_file, 500,1)
    #print(occ)
    #for o in occ:
    #    print(o)
if __name__ == '__main__':
    main()
