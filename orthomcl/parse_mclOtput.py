#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


from Bio import SeqIO
from optparse import OptionParser
import re
import os
from Bio import SeqUtils
import sys

def parse_orthomcl_output(orthomcl_file, orthofinder=False):
    orthomcl_groups2proteins= {}
    i = 0
    genome_orthomcl_code2proteins = {}
    protein2genome_ortho_mcl_code = {}
    protein_id2orthogroup_id = {}
    for line in open(orthomcl_file, 'r'):
        group = "group_%s" % str(i)
        orthomcl_groups2proteins[group] = []
        if not orthofinder:
            line = line.rstrip().split("\t")
        else:
            line = line.rstrip().split(' ')
            line = line[1:len(line)]
        for one_protein in line:
            # orthofinder format
            if not orthofinder:
                protein_id = one_protein.split("|")[1]
                genome_orthomcl_code = one_protein.split("|")[0]
                protein2genome_ortho_mcl_code[protein_id] = genome_orthomcl_code
                if genome_orthomcl_code not in genome_orthomcl_code2proteins.keys():
                    genome_orthomcl_code2proteins[genome_orthomcl_code] = [protein_id]
                else:
                    genome_orthomcl_code2proteins[genome_orthomcl_code].append(protein_id)

            else:
                protein_id = one_protein
            protein_id2orthogroup_id[protein_id] = group
            orthomcl_groups2proteins[group].append(protein_id)
        i+=1
    return  protein_id2orthogroup_id, orthomcl_groups2proteins, genome_orthomcl_code2proteins, protein2genome_ortho_mcl_code
