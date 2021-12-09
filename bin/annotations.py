from Bio import Entrez, SeqIO
from Bio.SeqUtils import CheckSum
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

import setup_chlamdb

import pandas as pd
import itertools
import sys
import re
import os

from collections import defaultdict



def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def count_missing_locus_tags(gbk_record):
    count_CDS = 0
    count_no_locus = 0
    for feature in gbk_record.features:
        if feature.type == 'CDS':
            count_CDS += 1
            if "locus_tag" not in feature.qualifiers:
                count_no_locus += 1
    return count_no_locus, count_CDS


def is_annotated(gbk_record):
    return not (len(gbk_record.features) == 1
            and gbk_record.features[0].type == 'source')


def rename_source(record):
    if 'strain' in record.features[0].qualifiers:

        print('--', record.features[0].qualifiers['strain'][0])
        if ';' in record.features[0].qualifiers['strain'][0]:
            print('ACHRTUNG: record has 2 strain names! \t%s\t --> check and edit source manually' % record.name)
            # put everythink lower size
            strain = record.features[0].qualifiers['strain'][0].split(';')[1]
        else:
            strain = record.features[0].qualifiers['strain'][0]
        if strain == 'strain':
            return (False, False)
        if strain.lower() not in record.annotations['source'].lower():
            msg = '%s' % record.annotations['source']
            print("ACHTUNG changing source\t%s\t--> %s " % (msg, record.annotations['source'] + strain))


        return strain, "%s %s" % (record.annotations['source'], strain)
    else:
        return (False, False)

def clean_description(description):
    import re
    description = re.sub(", complete genome\.", "", description)
    description = re.sub(", complete genome", "", description)
    description = re.sub(", complete sequence\.", "", description)
    description = re.sub("strain ", "", description)
    description = re.sub("str\. ", "", description)
    description = re.sub(" complete genome sequence\.", "", description)
    description = re.sub(" complete genome\.", "", description)
    description = re.sub(" chromosome", "", description)
    description = re.sub(" DNA", "", description)
    description = re.sub("Merged record from ", "", description)
    description = re.sub(", wgs", "", description)
    description = re.sub("Candidatus ", "", description)
    description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
    description = re.sub("complete genome, isolate", "", description)
    description = re.sub(" complete", "", description)
    description = re.sub(" genome assembly.*", "", description)
    description = re.sub("Chlamydophila", "Chlamydia", description)
    description = re.sub(", whole genome shotgun sequence", "", description)

    return description


def orthogroups_to_fasta(genomes_list):
    fasta_list = genomes_list.split(' ')

    sequence_data = {}
    for fasta_file in fasta_list:
        sequence_data.update(SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")))

      # write fasta
    with open("Orthogroups.txt") as f:
        all_grp = [i for i in f]
        for n, line in enumerate(all_grp):
            groups = line.rstrip().split(' ')
            group_name = groups[0][:-1]
            groups = groups[1:]
            if len(groups)>1:
                new_fasta = [sequence_data[i] for i in groups]
                out_path = "%s.faa" % group_name
                out_handle = open(out_path, "w")
                SeqIO.write(new_fasta, out_handle, "fasta")


def check_organism_uniqueness(gbk_lst):
    """
    As BioSQL uses the organism entry of the records to assign
    a taxid when it is not allowed to access the ncbi online,
    we need to ensure that the same organism name is not used
    twice in the genomes to prevent bioentries from different genbank
    files to be entered under the same taxon_id (as chlamdb uses the taxon_id to 
    differentiate between genomes).

    Display an error message and stop the pipeline if this arises
    """
    organisms = defaultdict(list)
    contigs = defaultdict(list)
    locuses = defaultdict(list)

    for gbk_file in gbk_lst:
        curr_organism = None
        for record in SeqIO.parse(gbk_file, "genbank"):
            if not "organism" in record.annotations:
                raise Exception(f"No organism for file {gbk_file}")

            for feature in record.features:
                if feature.type != "CDS":
                    continue
                if not "locus_tag" in feature.qualifiers:
                    raise Exception(f"Missing locus tag in {gbk_file}")
                locus_tag = feature.qualifiers["locus_tag"][0]
                locuses[locus_tag].append(record.name)
            contigs[record.name].append(gbk_file)
            organism = record.annotations["organism"]
            if curr_organism is None:
                curr_organism = organism
            elif curr_organism != organism:
                raise Exception(f"Two different organism in {gbk_file}: {curr_organism}/{organism}")
        organisms[curr_organism].append(gbk_file)

    for locus, contig_list in locuses.items():
        if len(contig_list)>1:
            err = ",".join(contig_list)
            raise Exception(f"Duplicate locus tag {locus} in contigs {err}")

    for contig, file_lst in contigs.items():
        if len(file_lst)>1:
            genbanks = ",".join(file_lst)
            raise Exception(
                f"Duplicate contig name: {contig} in "
                f"{genbanks}")

    for organism, file_lst in organisms.items():
        if len(file_lst) > 1:
            genbanks = ",".join(file_lst)
            raise Exception(f"Genbank files {genbanks} have the same organism: {organism}")


def check_names(record, common_names):
    common_name = record.annotations.get("source", None)
    if (common_name is None) or (common_name in common_names):
        # such a common name has already been 
        # encountered: change it to the unique scientific
        # name
        # NOTE: We know that record.annotations contains "organism"
        record.annotations["source"] = record.annotations["organism"]
    return record


def check_organism_names(gbk_files):
    """
    NOTE: assumes that the check_organism_uniqueness has been called before (i.e.
    unique scientific organism name and only one organism per file).

    This function does two things:
    ** filter out unannotated contigs (TODO: check if this is really necessary)
    ** for genbank files that have the same organism common name, change it to their
       scientific name (BioSQL would ortherwise assign them the same taxon_id,
       I know, this is stupid)
    """
    common_names = set()
    for gbk_file in gbk_files:
        records = SeqIO.parse(gbk_file, "genbank")
        result_file = gbk_file.replace(".gbk", "_filtered.gbk")
        to_keep = (check_names(rec, common_names) for rec in records)
        first_rec = next(to_keep)
        common_names.add(first_rec.annotations["source"])
        SeqIO.write(itertools.chain([first_rec], to_keep), result_file, "genbank")


def convert_gbk_to_fasta(gbf_file, edited_gbf, output_fmt="faa", keep_pseudo=False):
    records = SeqIO.parse(gbf_file, 'genbank')
    edited_records = open(edited_gbf, 'w')

    for record in records:
        for feature in record.features:
            if feature.type == 'CDS':
                if not keep_pseudo and ('pseudo' in feature.qualifiers 
                    or 'pseudogene' in feature.qualifiers):
                    continue

                if "locus_tag" not in feature.qualifiers:
                    raise Exception(f"Feature without locus tag in record {record.name}")

                locus_tag = feature.qualifiers["locus_tag"][0]
                if output_fmt=="faa":
                    if "translation" not in feature.qualifiers:
                        continue
                    data = feature.qualifiers["translation"][0]
                elif output_fmt=="fna":
                    data = feature.location.extract(record).seq
                else:
                    raise Exception(f"Unsupported option: {output_fmt}, must be either faa or fna")

                edited_records.write(">%s %s\n%s\n" % (locus_tag,
                                                 record.description, data))


def convert_gbk_to_fna(gbf_file, fna_contigs):
    records = SeqIO.parse(gbf_file, 'genbank')  #from BioPython, object of class SeqRecord
    edited_records = open(fna_contigs, 'w')

    for record in records:
       edited_records.write(">%s %s\n%s\n" % (record.name, record.description, record.seq))


# all faa files are merged into fasta_file
def get_nr_sequences(fasta_file, genomes_list):
    locus2genome = {}
    for fasta in genomes_list:
        genome = os.path.basename(fasta).split('.')[0]
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome
    nr_fasta = open('nr.faa', 'w')
    nr_mapping = open('nr_mapping.tab', 'w')

    hsh_checksum_list = {}

    records = SeqIO.parse(fasta_file, "fasta")
    updated_records = []

    for record in records:

        # NOTE: the case is important for crc64, need to check whether it
        # is necessary to make all entries lower/upper case to ensure consistency.
        checksum = CheckSum.crc64(record.seq)
        nr_mapping.write("%s\t%s\t%s\n" % (record.id,
                                          checksum,
                                          locus2genome[record.id]))
        if checksum not in hsh_checksum_list:
            hsh_checksum_list[checksum] = [record]
            record.id = checksum
            record.name = ""
            updated_records.append(record)
        else:
            # NOTE: having same hash does not mean that the sequences are identical: as
            # the hash space is smaller than the sequence space, it means that collision
            # are unavoidable (but not probable) and record with same hashes should be compared
            # https://www.uniprot.org/help/uniparc (sequence comparison)
            #
            # the list of records having the same checksum, but potentially, 
            # different sequences -> compare them: python does so 
            # comparing the sequences as strings, assuming a similar alphabet
            lst_records = hsh_checksum_list[checksum]
            sequence = record.seq
            has_identical = False
            for prev_record in lst_records:
                if prev_record.seq == sequence:
                    has_identical = True
                    break
            if not has_identical:
                lst_records.append(record)
                record.id = checksum + "-" + len(lst_records)
                record.name = ""
                updated_records.append(record)

    SeqIO.write(updated_records, nr_fasta, "fasta")



# This function parses the result of orthofinder/orthoMCL, 
# retrieves the groups of orthologs detected by the tools
# and returns the core single copy orthologs (i.e. the 
# ortholog present in all samples)
#
# Orthofinder output file : Orthogroup_ID: locus1 locus2... locusN
def orthofinder2core_groups(fasta_list,
                          mcl_file,
                          n_missing=0):
    orthogroup2locus_list = {}
    with open(mcl_file, 'r') as f:
        for line in f:
            groups = line.rstrip().split(' ')
            group_id = groups[0][:-1] # remove lagging ':'
            groups = groups[1:]
            orthogroup2locus_list[group_id] = groups

    locus2genome = {}
    for fasta in fasta_list:
        tokens = os.path.basename(fasta).split('.')
        # the filename may contain additional "." that are to be included
        genome = ".".join(tokens[:-1])
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome

    df = pd.DataFrame(index=orthogroup2locus_list.keys(), columns=set(locus2genome.values()))
    df = df.fillna(0)

    # NOTE: to be replace by a map-reduce algo, to avoid
    # a slow for loop in python and to make the code more compact
    for group_id,loci_list in orthogroup2locus_list.items():
        for locus in loci_list:
            genome = locus2genome[locus]
            df.loc[group_id, genome] += 1

    df =df.apply(pd.to_numeric, args=('coerce',))

    n_genomes = len(set(locus2genome.values()))
    n_minimum_genomes = n_genomes-n_missing
    freq_missing = (n_genomes-float(n_missing))/n_genomes
    limit = freq_missing*n_genomes

    groups_with_paralogs = df[(df > 1).sum(axis=1) > 0].index
    df = df.drop(groups_with_paralogs)

    core_groups = df[(df == 1).sum(axis=1) >= limit].index.tolist()
    return core_groups, orthogroup2locus_list, locus2genome


def get_core_orthogroups(genomes_list, int_core_missing):
    core_groups, orthogroup2locus_list, locus2genome = orthofinder2core_groups(genomes_list,
          'Orthogroups.txt', int_core_missing)

    for group_id in core_groups:
        sequence_data = SeqIO.to_dict(SeqIO.parse(group_id + "_mafft.faa", "fasta"))
        dest = group_id + '_taxon_ids.faa'
        new_fasta = []
        for locus in orthogroup2locus_list[group_id]:
            tmp_seq = sequence_data[locus]
            tmp_seq.name = locus2genome[locus]
            tmp_seq.id = locus2genome[locus]
            tmp_seq.description = locus2genome[locus]
            new_fasta.append(tmp_seq)

        out_handle = open(dest, 'w')
        SeqIO.write(new_fasta, out_handle, "fasta")
        out_handle.close()


def calculate_og_identities(input_fasta, output_file):
    alignment = AlignIO.read(input_fasta, "fasta")
    values = []
    for i in range(len(alignment)):
        for j in range(i+1, len(alignment)):
            alignment_1 = alignment[i]
            alignment_2 = alignment[j]
            locus_tag_1 = alignment_1.name
            locus_tag_2 = alignment_2.name
            alignment_length = 0
            identical = 0
            for ch1, ch2 in zip(alignment_1, alignment_2):
                if ch1=="-" or ch2=="-":
                    continue
                if ch1==ch2:
                    identical += 1
                alignment_length += 1

            if alignment_length>0:
                per_identity = 100*(identical/float(alignment_length))
            else:
                per_identity = 0
            values.append((locus_tag_1, locus_tag_2, per_identity, alignment_length))
    output_fh = open(output_file, "w")
    for lt1, lt2, ident, le in values:
        print(lt1, lt2, ident, le, sep=",", file=output_fh)
    output_fh.close()


def concatenate_core_orthogroups(fasta_files):
    out_name = 'msa.faa'

    taxons = []
    all_seq_data = {}
    fasta_to_algn_length = {}
    for one_fasta in fasta_files:
        all_seq_data[one_fasta] = {}
        for record in AlignIO.read(one_fasta, "fasta"):
            if one_fasta not in fasta_to_algn_length:
                fasta_to_algn_length[one_fasta] = len(record)

            if record.id not in taxons:
                taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record

    concat_data = defaultdict(str)
    for one_fasta in fasta_files:
        for taxon in taxons:
            if taxon not in all_seq_data[one_fasta]:
                seq = "-"*fasta_to_algn_length[one_fasta]
            else:
                seq = all_seq_data[one_fasta][taxon]
            concat_data[taxon] += seq

    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")

