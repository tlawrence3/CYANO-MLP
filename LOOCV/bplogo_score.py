import sys
import re
import random
import argparse
from Bio import SeqIO
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--bootstrap", "-b", type=int)
parser.add_argument("--basepair", "-p", action="store_true")
parser.add_argument("--anti", "-a", action="store_true")
parser.add_argument("groups", nargs="+")
args = parser.parse_args()


sequence_offset = 1
groups = args.groups[:-1]
states = ["A", "G", "C", "U"]
bp_states = ["AA", "AU", "AG", "AC", "GA", "GU", "GG", "GC", "CA", "CU", "CG", "CC", "UA", "UU", "UG", "UC"]
logo_scores = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(float))))
bp_logo_scores = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(float))))
anti_logo_scores = defaultdict(lambda : defaultdict(lambda : defaultdict(lambda : defaultdict(float))))

species_scores = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
bp_species_scores = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))
anti_species_scores = defaultdict(lambda : defaultdict(lambda : defaultdict(float)))

total_bits = defaultdict(lambda: defaultdict(float))

for group in groups:
    for state in states:
        with open("{}_{}.eps".format(state, group), "rU") as logo:
            datastart = False
            for line in logo:
                line = line.strip()
                data_start_match = re.match("\%\%EndProlog", line)
                num_match = re.match("numbering\s+\{\((\d+)\)", line)
                aa_height_match = re.match("^(\d\.\d+)\s\(([A-Z])\)\snumchar", line)
                if(data_start_match):
                    datastart = True
                elif (num_match):
                    coord = int(num_match.group(1)) - sequence_offset
                elif (aa_height_match):
                    ht = float(aa_height_match.group(1))
                    aa = aa_height_match.group(2)
                    logo_scores[group][coord][state][aa] = ht
                    total_bits[group][aa] += ht

if (args.basepair):
    for group in groups:
        for state in bp_states:
            with open("anti_{}_{}.eps".format(state, group), "rU") as bplogo:
                datastart = False
                for line in bplogo:
                    line = line.strip()
                    data_start_match = re.match("\%\%EndProlog", line)
                    num_match = re.match("numbering\s+\{\((\d+:\d+)\)", line)
                    aa_height_match = re.match("^(\d\.\d+)\s\(([A-Z])\)\snumchar", line)
                    if(data_start_match):
                        datastart = True
                    elif (num_match):
                        coord = num_match.group(1)
                    elif (aa_height_match):
                        ht = float(aa_height_match.group(1))
                        aa = aa_height_match.group(2)
                        bp_logo_scores[group][coord][state][aa] = ht
                        total_bits[group][aa] += ht


if (args.anti):
    for group in groups:
        for state in bp_states:
            with open("{}_{}.eps".format(state, group), "rU") as bplogo:
                datastart = False
                for line in bplogo:
                    line = line.strip()
                    data_start_match = re.match("\%\%EndProlog", line)
                    num_match = re.match("numbering\s+\{\((\d+)\)", line)
                    aa_height_match = re.match("^(\d\.\d+)\s\(([A-Z])\)\snumchar", line)
                    if(data_start_match):
                        datastart = True
                    elif (num_match):
                        coord = int(num_match.group(1)) - 1
                    elif (aa_height_match):
                        ht = float(aa_height_match.group(1))
                        aa = aa_height_match.group(2)
                        anti_logo_scores[group][coord][state][aa] = ht
                        total_bits[group][aa] += ht

num_tRNA_species = defaultdict(int)
clade_dict = {}
with open(args.groups[-1], "rU") as sequence_file:
    for seq_record in SeqIO.parse(sequence_file, "fasta"):
        id_match = re.match("^(\w).*\|(\w\w*)\|$", seq_record.id)
        aa = id_match.group(1)
        clade = id_match.group(2)
        species = seq_record.id.split("|")[1]
        clade_dict[species] = clade
        num_tRNA_species[species] += 1
        nucs = set("ACGT")
        seq_length = len(seq_record.seq)
        for coord, nucleotide in enumerate(seq_record.seq.upper()):
            if (nucleotide in nucs):
                for group in groups:
                    if (logo_scores[group][coord][nucleotide][aa]):
                        species_scores[species][group][coord] += logo_scores[group][coord][nucleotide][aa]


if (args.anti):
    with open(args.groups[-1], "rU") as sequence_file:
        for seq_record in SeqIO.parse(sequence_file, "fasta"):
            id_match = re.match("^(\w).*\|(\w\w*)\|$", seq_record.id)
            aa = id_match.group(1)
            clade = id_match.group(2)
            species = seq_record.id.split("|")[1]
            nucs = set("ACGT")
            seq_length = len(seq_record.seq)
            for coord, nucleotide in enumerate(seq_record.seq.upper()):
                if (nucleotide in nucs):
                    for group in groups:
                        if (anti_logo_scores[group][coord][nucleotide][aa]):
                            anti_species_scores[species][group][coord] += anti_logo_scores[group][coord][nucleotide][aa]


if (args.basepair):
    with open(args.groups[-1], "rU") as sequence_file:
        for seq_record in SeqIO.parse(sequence_file, "fasta"):
            id_match = re.match("^(\w).*\|(\w\w*)\|$", seq_record.id)
            aa = id_match.group(1)
            clade = id_match.group(2)
            species = seq_record.id.split("|")[1]
            for key in bp_logo_scores:
                for coord in bp_logo_scores[key]:
                    coord1, coord2 = coord.split(":")
                    coord1 = int(coord1) - sequence_offset
                    coord2 = int(coord2) - sequence_offset
                    nuc1 = seq_record.seq.upper()[coord1]
                    nuc2 = seq_record.seq.upper()[coord2]
                    state = str(nuc1 + nuc2)
                    for group in groups:
                        if (bp_logo_scores[group][coord][state][aa]):
                            bp_species_scores[species][group][coord] += bp_logo_scores[group][coord][nucleotide][aa]

for species in sorted(species_scores):
    print("{}\t".format(clade_dict[species]), end = "")
    for group in groups:
        spec_score = sum(species_scores[species][group].values())
        if (args.anti):
            spec_score += sum(anti_species_scores[species][group].values())

        if (args.basepair):
            spec_score += sum(bp_species_scores[species][group].values())

        print("{:0.2f}\t".format(spec_score/num_tRNA_species[species]), end = "")
    print("{}".format(species))

if (args.bootstrap):
    for rep in range(args.bootstrap):
        positions = random.choices(range(seq_length), k = seq_length)
        for species in sorted(species_scores):
            print("boot{}:{}\t".format(rep, clade_dict[species]), end = "")
            for group in groups:
                score = 0
                for pos in positions:
                    score += species_scores[species][group][pos]
                score = score/num_tRNA_species[species]
                print("{:0.2f}\t".format(score), end="")
            print("{}".format(species))
