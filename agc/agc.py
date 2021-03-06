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

__author__ = "BOUARROUDJ Lisa & SOULA Sarah"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["BOUARROUDJ Lisa & SOULA Sarah"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "BOUARROUDJ Lisa & SOULA Sarah"
__email__ = "lisa.bdj.95@gmail.com & spequez@gmail.com"
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
    """Lit un fichier fasta.gz et renvoie les s??quences.
      :Param??tres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des s??quences
      Retourne:
          G??n??rateur de s??quences de longueur minimale
    """
    with gzip.open(amplicon_file, "rt") as filin:
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
    """Identifie les s??quences uniques avec un nombre d'occurrence minimal.
      :Param??tres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des s??quences
          mincount: nombre minimal d'occurrence
      Retourne:
          G??n??rateur de s??quences uniques ayant un nombre d'occurrence minimal
    """
    seqs = read_fasta(amplicon_file, minseqlen)
    occurrences= {}
    for seq in seqs:
        if seq in occurrences:
            occurrences[seq]+=1
        else:
            occurrences[seq]=1
    occurrences = sorted(occurrences.items(), key=lambda t: t[1] , reverse=True)
    for i in range(len(occurrences)):
        if occurrences[i][1] >= mincount:
            yield [occurrences[i][0],occurrences[i][1]]



def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2):
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Donne une liste de chunks non chevauchants"""
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
    """Prend en une liste de s??quences align??es au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Recense tous les k-mers pr??sents dans les s??quences non-chim??riques.
      :Param??tres:
          kmer_dict: dictionnaire de kmers pr??sents dans les s??quences non-chim??riques
          sequence: s??quence
          id_seq: identifiant de la s??quence
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
    """Cherche les 2 s??quences parentes.
      :Param??tres:
          kmer_dict: dictionnaire de kmers pr??sents dans les s??quences non-chim??riques
          sequence: s??quence
          kmer_size: longueur des kmers
      Retourne:
          liste_seq: identifiants des s??quences parentes
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
    """Indique si la s??quence est une chim??re, ou non.
      :Param??tre:
          perc_identity_matrix: matrice avec le taux d'identit?? entre la
          s??quence candidate et deux s??quences parentes
      Retourne:
          bool??en indiquant si la s??quence candidate est une chim??re ou non
    """
    mat = []
    seq = 1
    seq_bool = 0
    for i in range (len(perc_identity_matrix)):
        mat.append(statistics.stdev(perc_identity_matrix[i]))
        if statistics.mean(mat) > 5:
            if perc_identity_matrix[0][0] < perc_identity_matrix[0][1]:
                seq = 2
            for i in range (1,len(perc_identity_matrix)):
                if perc_identity_matrix[i][0] < perc_identity_matrix[i][1]:
                    seq_bool = 2
                else:
                    seq_bool = 1
                if seq_bool != seq:
                    return True
    return False


def std(data):
    """Calcule l'??cart-type.
      :Param??tre:
          data: donn??e
      Retourne:
          std: ??cart-type des donn??es
    """
    std = statistics.stdev(data)
    return std


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Retourne les s??quences non-chim??riques
      :Param??tres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des s??quences
          mincount: comptage minimum
          chunk_size: longueur des chunks
          kmer_size: longueur des kmers
      Retourne:
          G??n??rateur de s??quences non-chim??riques
    """
    gene_seq = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
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
        if not  detect_chimera(mat_ident):
            non_chimer.append(gene_seq[i])
        mat_ident = []
    for i in range(len(non_chimer)):
        yield non_chimer[i]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Retourne la liste d'OTU.
      :Param??tres:
          amplicon_file: fichier au format fasta.gz
          minseqlen: longueur minimale des s??quences
          mincount: comptage minimum
          chunk_size: longueur des chunks
          kmer_size: longueur des kmers
      Retourne:
          list_otu: liste d'OTU avec le nombre d'occurences
    """
    sequence_non_chimer = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
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
    return otu_liste


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """Cr??ation du fichier de sortie.
      :Param??tres:
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

    # Programme
    liste_otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(liste_otu, args.output_file)

if __name__ == '__main__':
    main()
