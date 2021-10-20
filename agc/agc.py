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
        print("ouvre fichier")
        seq = ""
        #count = 0
        for line in filin:
            if line[0] == ">":
                if len(seq) > minseqlen:
                    yield seq
                seq= ""
            elif len(line) > 0:
                seq += line.strip()
            #count+=1
            #if count > 1000 :
            #    break
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


    seqs = read_fasta(amplicon_file, minseqlen)
    occurrences= {}

    for seq in seqs:
        if seq in occurrences:
            occurrences[seq]+=1
        else:
            occurrences[seq]=1

    occurrences = sorted(occurrences.items(), key=lambda t: t[1] , reverse=True)
    #print("Length :", len(occurrences))
    print("liste crée")
    for i in range(len(occurrences)):
        if occurrences[i][1] >= mincount:
            yield [occurrences[i][0],occurrences[i][1]]



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

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
      :Paramètres:
          kmer_dict: dictionnaire de kmers présents dans les séquences non-chimériques
          sequence: séquence
          id_seq: identifiant de la séquence
          kmer_size: longueur des kmers
      Retourne:
          kmer_dict: liste d'OTU avec le nombre d'occurences
    """
    gene_kmer = cut_kmer(sequence, kmer_size)
    for kmer in gene_kmer:
        if kmer in kmer_dict:
            kmer_dict[kmer].append(id_seq)
        else:
            kmer_dict[kmer]=[id_seq]
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
      :Paramètres:
          kmer_dict: dictionnaire de kmers présents dans les séquences non-chimériques
          sequence: séquence
          kmer_size: longueur des kmers
      Retourne:
          liste_seq: identifiants des séquences parentes
    """
    gene_kmer = cut_kmer(sequence,kmer_size)
    liste_kmer=[]
    liste_seq =[]
    for kmer in gene_kmer:
        if kmer in kmer_dict:
            for i in range(len(kmer_dict[kmer])):
                liste_kmer.append(kmer_dict[kmer][i])
    c = Counter(liste_kmer).most_common(2)
    liste_seq.append(c[0][0])
    liste_seq.append(c[1][0])
    return liste_seq

def detect_chimera(perc_identity_matrix):
    """
      :Paramètre:
          perc_identity_matrix: matrice avec le taux d'identité entre la
          séquence candidate et deux séquences parentes
      Retourne:
          booléen indiquant si la séquence candidate est une chimère ou non
    """

    mat = []
    seq=1
    seq_bool=0
    for i in range (len(perc_identity_matrix)):
        mat.append(statistics.stdev(perc_identity_matrix[i]))
        if statistics.mean(mat)>5:
            if perc_identity_matrix[0][0]<perc_identity_matrix[0][1]:
                seq=2
            for i in range (1,len(perc_identity_matrix)):
                if perc_identity_matrix[i][0]<perc_identity_matrix[i][1]:
                    seq_bool=2
                else:
                    seq_bool=1
                if seq_bool!=seq:
                    return True
    return False


def std(data):
    """Calcule l'écart-type.
      :Paramètre:
          data: donnée
      Retourne:
          std: écart-type des données
    """
    std = statistics.stdev(data)
    return std


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
      :Paramètres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des séquences
          mincount: comptage minimum
          chunk_size: longueur des chunks
          kmer_size: longueur des kmers
      Retourne:
          :
    """
    gene_seq = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    print("dereplication")
    print("Length derep :", len(gene_seq))
    non_chimer = []
    non_chimer.append(gene_seq[0])
    non_chimer.append(gene_seq[1])
    kmer_dict = {}
    mat_ident= []
    for i in range(2,len(gene_seq)):
        kmer_dict=  get_unique_kmer(kmer_dict, gene_seq[i][0], i, kmer_size)
    for i in range(2,len(gene_seq)):
        ide = search_mates(kmer_dict, gene_seq[i][0], kmer_size)
        seq0= get_chunks(gene_seq[i][0], chunk_size)
        seq1= get_chunks(gene_seq[ide[0]][0], chunk_size)
        seq2= get_chunks(gene_seq[ide[1]][0], chunk_size)
        for j in range(4):
            align1 = nw.global_align(seq0[j],seq1[j], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            align2 = nw.global_align(seq0[j],seq2[j], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            mat_ident.append([get_identity([align1[0],align1[1]]),get_identity([align2[0],align2[1]])])
            print()
        if not  detect_chimera(mat_ident):
            non_chimer.append(gene_seq[i])
            print("chimera length :", len(non_chimer))
        mat_ident = []
    for i in range(len(non_chimer)):
        yield non_chimer[i]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
      :Paramètres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des séquences
          mincount: comptage minimum
          chunk_size: longueur des chunks
          kmer_size: longueur des kmers
      Retourne:
          list_otu: liste d'OTU avec le nombre d'occurences
    """

    sequence_non_chimer = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    print("abudance!!!")
    otu_liste = []
    otu_liste.append(sequence_non_chimer[0])
    for i in range(1,len(sequence_non_chimer)):
        flag = True
        for j in range(len(otu_liste)):
            align = nw.global_align(sequence_non_chimer[i][0],otu_liste[j][0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            if get_identity(align) > 97:
                flag = False
                break
            
        if flag is True:
            otu_liste.append(sequence_non_chimer[i])
            print(len(otu_liste))

    return otu_liste


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """
      :Paramètres:
          OTU_list: liste d'OTU
          output_file: chemin vers le fichier de sortie
    """
    with open(output_file,"w") as filout:
        for i in range((len(OTU_list))):
            filout.write(f">OTU_{i+1} occurrence:{OTU_list[i][1]}\n")
            filout.write(fill(OTU_list[i][0], width=80))
            filout.write("\n")


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
    liste_otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(liste_otu, args.output_file)

if __name__ == '__main__':
    main()
