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
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
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
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, 'rt') as amplicon:
        sequence = ""
        for ligne in amplicon:
            if ligne.startswith(">"):
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
            else:
                sequence += ligne.strip()

        # Dernière séquence
        if len(sequence) >= minseqlen:
            yield sequence


def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count]
    of sequence with a count >= mincount and a length >= minseqlen.
    """
    seq_occ = {}
    with gzip.open(amplicon_file, 'rt') as amplicon:
        for sequence in read_fasta(amplicon_file, minseqlen):
            if sequence not in seq_occ.keys():
                seq_occ[sequence] = 1
            else:
                seq_occ[sequence] += 1

    # Ordre décroissant d'occurrence
    seq_occ = dict(sorted(seq_occ.items(), key=lambda item: item[1], reverse=True))

    # Vérification du nombre d'occurence
    for sequence, count in seq_occ.items():
        if count >= mincount:
            yield [sequence, count]

def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the
    format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    seq_a = alignment_list[0]
    seq_b = alignment_list[1]

    len_align = len(seq_a)
    nb_id = 0
    for i in range(len_align):
        if seq_a[i] == seq_b[i]:
            if seq_a[i] != '-':
                nb_id += 1

    percent = nb_id / len_align * 100
    return percent

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int,
                                mincount: int, chunk_size: int,
                                kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    id_percent = 97
    seq_list = []
    for sequence in dereplication_fulllength(amplicon_file, minseqlen, mincount):
        seq_list.append(sequence)

    OTU_list = [seq_list[0]]
    for sequence in seq_list:
        is_OTU = True
        for reference in OTU_list:
            align = nw.global_align(sequence[0], reference[0],
                                    gap_open=-1, gap_extend=-1,
                                    matrix=str(Path(__file__).parent / "MATCH"))
            identity = get_identity(align)

            # Identité avec la référence
            if identity >= id_percent:
                is_OTU = False
                break

        if is_OTU:
            OTU_list.append(sequence)

    return OTU_list


def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, 'w') as output:
        for i, (sequence, count) in enumerate(OTU_list):
            output.write(f">OTU_{i+1} occurrence:{count}\n")

            formatted = textwrap.fill(sequence, width=80)
            output.write(f"{formatted}\n")


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici

    # Lecture des arguments
    amplicon_file = args.amplicon_file
    minseqlen = args.minseqlen
    mincount = args.mincount
    chunk_size = 0 # Valeur nulle car inutile
    kmer_size = 0 # Valeur nulle car inutile
    output_file = args.output_file

    # Récupération de la liste d'OTU
    OTU_list = abundance_greedy_clustering(amplicon_file, minseqlen,
                                           mincount, chunk_size,
                                           kmer_size)

    # Sauvegarde des OTU obtenus
    write_OTU(OTU_list, output_file)


if __name__ == '__main__':
    main()

    # On trouve :
    # Matching unique query sequences: 116 of 117 (99.15%)