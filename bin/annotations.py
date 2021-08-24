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
import gzip
import re
import os



def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def merge_gbk(gbk_records, filter_size=0, gi=False, plasmids=False):
    '''
    merge multiple contigs into a single DNA molecule with 200*N between contigs
    keep source description from the first record
    remove contigs smaller than <filter_size>

    For plasmids: keep the first accession, as having several entries with the 
    same accession would profoundly hurt BioSQL.

    :param gbk_records:
    :param filter_size:
    :param gi:
    :return:
    '''

    n=0
    if len(gbk_records) == 1:
        merged_rec = gbk_records[0]
    else:
        for i, rec in enumerate(gbk_records):
            # remove source feature of all records except the first one
            if rec.features[0].type == 'source' and i != 0:
                rec.features.pop(0)
            # filter small contigs
            if len(rec) > filter_size:
                if n == 0:
                    n+=1
                    merged_rec = rec
                else:
                    merged_rec+=rec
                # you could insert a spacer if needed
                # do not add spacer after the last contig
                if i != len(gbk_records)-1:
                    merged_rec += "N" * 200

                    my_start_pos = ExactPosition(len(merged_rec)-200)
                    my_end_pos = ExactPosition(len(merged_rec))
                    my_feature_location = FeatureLocation(my_start_pos, my_end_pos)
                    my_feature = SeqFeature(my_feature_location, type="assembly_gap")
                    merged_rec.features.append(my_feature)

    accession_index = 0 if plasmids else -1
    try:
        merged_rec.id = gbk_records[0].annotations["accessions"][accession_index]
    except KeyError:
        merged_rec.id = gbk_records[0].id

    if gi:
        merged_rec.annotations["gi"] = gi

    merged_rec.description = "%s" % gbk_records[0].annotations["organism"]
    merged_rec.annotations = gbk_records[0].annotations
    try:
        merged_rec.name = gbk_records[0].annotations["accessions"][accession_index]
    except KeyError:
        merged_rec.name = gbk_records[0].id
    my_start_pos = ExactPosition(0)
    my_end_pos = ExactPosition(len(merged_rec))
    merged_rec.features[0].location = FeatureLocation(my_start_pos, my_end_pos)
    return merged_rec


def filter_plasmid(record_list):
    plasmid_record_list = []
    chromosome_record_list = []

    for record in record_list:
        if setup_chlamdb.is_plasmid(record):
            plasmid_record_list.append(record)
        else:
            chromosome_record_list.append(record)
    return chromosome_record_list, plasmid_record_list


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


def update_record_taxon_id(record, n):
    if record.features[0].type == 'source':
        if 'db_xref' in record.features[0].qualifiers:
            for item in record.features[0].qualifiers['db_xref']:
                if 'taxon' in item:
                    index = record.features[0].qualifiers['db_xref'].index(item)
                    record.features[0].qualifiers['db_xref'][index] = "taxon:%s" % n
    else:
        print('ACHRTUNG\t no source for record \t%s' % record.name)
    return record


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


def check_gbk(gbff_files):
    reannotation_list = []

    # count the number of identical source names
    source2count = {}
    accession2count = {}
    for genome_number, gbff_file in enumerate(gbff_files):

        records = list(SeqIO.parse(gbff_file, "genbank"))
        for record in records:
            n_missing, total = count_missing_locus_tags(record)
            if n_missing > 0:
                print ('Warning: %s/%s missing locus tag for record %s' % (n_missing, total, record.name))
        chromosome, plasmids = filter_plasmid(records)
        cleaned_records = []
        plasmid_reannot = False
        chromosome_reannot = False

        for n_plasmid, plasmid in enumerate(plasmids):
            annot_check = is_annotated(plasmid)
            if annot_check:

                plasmid.description = clean_description(plasmid.description)

                plasmid = update_record_taxon_id(plasmid, 1000 + genome_number)
                strain, new_source = rename_source(plasmid)
                print("plasmid:", strain, new_source )
                if new_source:
                    if not 'plasmid' in new_source:
                        new_source = "%s plasmid %s" % (new_source, n_plasmid+1)
                    if strain.lower() not in plasmid.annotations['source'].lower():
                        plasmid.description = new_source
                    if strain.lower() not in plasmid.annotations['organism'].lower():
                        plasmid.annotations['organism'] = new_source
                    if strain.lower() not in plasmid.annotations['source'].lower():
                        plasmid.annotations['source'] = new_source
                else:
                    print ('ACHTUNG\t no strain name for \t%s\t, SOUCE uniqueness should be checked manually' % merged_record.id)
                # check if accession is meaningful
                if 'NODE_' in plasmid.id or 'NODE_' in plasmid.name:
                    print ('ACHTUNG\t accession probably not unique (%s) for \t%s\t --> should be checked manually' % (merged_record.id))
                cleaned_records.append(plasmid)
            else:
                plasmid_reannot = True
                print("Warrning: unannotated genome: %s" % plasmid)

        if len(chromosome) > 0:
            '''
            Assume single chromosome bacteria.
            If multiple record founds, consider those as contigs.
            '''

            annot_check = is_annotated(chromosome[0])
            if annot_check:
                if len(chromosome) > 1:
                    merged_record = merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]
                merged_record.description = clean_description(merged_record.description)
                # rename source with strain name
                merged_record = update_record_taxon_id(merged_record, 1000 + genome_number)
                strain, new_source = rename_source(merged_record)
                if new_source:
                    if strain.lower() not in merged_record.annotations['source'].lower():
                        merged_record.description = new_source
                    if strain.lower() not in merged_record.annotations['organism'].lower():
                        merged_record.annotations['organism'] = new_source
                    if strain.lower() not in merged_record.annotations['source'].lower():
                        merged_record.annotations['source'] = new_source
                else:
                    print('ACHTUNG\t no strain name for\t%s' % gbff_file)
                # check if accession is meaningful
                if 'NODE_' in merged_record.id or 'NODE_' in merged_record.name:
                    print('ACHTUNG\t accession probably not unique (%s) for\t%s' % (merged_record.id, gbff_file))
                cleaned_records.append(merged_record)
            else:
                chromosome_reannot = True
                print("Warrning: unannotated genome: %s" % chromosome)


        if plasmid_reannot and not chromosome_reannot and len(chromosome) > 0:
            print("plasmid", plasmid_reannot)
            print("chr", chromosome_reannot)
            raise TypeError('Combination of unannotated plasmid(s) and annotated chromosome')
        elif not plasmid_reannot and chromosome_reannot and len(plasmids) > 0:
            print("plasmid", plasmid_reannot)
            print("chr", chromosome_reannot)
            raise TypeError('Combination of annotated plasmid(s) and unannotated chromosome')
        elif plasmid_reannot or chromosome_reannot:
            print("plasmid", plasmid_reannot)
            print("chr", chromosome_reannot)
            raise TypeError('Some genomes are not annotated!')
        else:
            out_name = gbff_file.split('.')[0] + '_merged.gbk'
            with open(out_name, 'w') as f:
                SeqIO.write(cleaned_records, f, 'genbank')


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
    organisms = dict()

    for gbk_file in gbk_lst:
        curr_organism = None
        for record in SeqIO.parse(gbk_file, "genbank"):
            if not "organism" in record.annotations:
                raise Exception(f"No organism for file {gbk_file}")
            
            organism = record.annotations["organism"]
            if curr_organism is None:
                curr_organism = organism
            elif curr_organism != organism:
                raise Exception(f"Two different organism in {gbk_file}: {curr_organism}/{organism}")

        gbk_lst = organisms.setdefault(organism, [])
        gbk_lst.append(gbk_file)

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


def filter_out_unannotated(gbk_files):
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
        to_keep = (check_names(rec, common_names) for rec in records if is_annotated(rec))
        first_rec = next(to_keep)
        common_names.add(first_rec.annotations["source"])
        SeqIO.write(itertools.chain([first_rec], to_keep), result_file, "genbank")


def convert_gbk_to_faa(gbf_file, edited_gbf, output_fmt="faa", keep_pseudo=False):
    records = SeqIO.parse(gbf_file, 'genbank')
    edited_records = open(edited_gbf, 'w')

    for record in records:
        for feature in record.features:
            if (feature.type == 'CDS'
                    and 'pseudo' not in feature.qualifiers
                    and 'pseudogene' not in feature.qualifiers):

                if "locus_tag" not in feature.qualifier:
                    raise Exception(f"Feature without locus tag in record {record.name}")

                locus_tag = feature.qualifiers["locus_tag"][0]
                if output_fmt=="faa":
                    data = feature.seq
                elif output_fmt=="fna":
                    data = feature.qualifiers["translation"][0]
                else:
                    raise Exception(f"Unsupported option: {output_fmt}, must be either faa or fna")

                edited_records.write(">%s %s\n%s\n" % (locus_tag,
                                                 record.description, data))


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
                          n_missing=0,
                          orthomcl=False):
    orthogroup2locus_list = {}
    with open(mcl_file, 'r') as f:
        for line in f:
            if orthomcl:
                # not sure this code will work with orthoMCL, to be checked
                groups = line.rstrip().split('\t')
                groups = [i.split('|')[1] for i in groups]
            else:
                groups = line.rstrip().split(' ')
                group_id = groups[0][:-1] # remove lagging ':'
                groups = groups[1:]
            orthogroup2locus_list[group_id] = groups

    locus2genome = {}
    for fasta in fasta_list:
        genome = os.path.basename(fasta).split('.')[0]
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
          'Orthogroups.txt', int_core_missing, False)

    for group_id in core_groups:
    # sequence_data = SeqIO.to_dict(SeqIO.parse("OG{0:07d}_mafft.faa".format(int(one_group.split('_')[1])), "fasta"))
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

def concatenate_core_orthogroups(fasta_files):
    out_name = 'msa.faa'

    # identification of all distinct fasta headers id (all unique taxons ids) in all fasta
    # storing records in all_seq_data (dico)
    taxons = []
    all_seq_data = {}
    for one_fasta in fasta_files:
        all_seq_data[one_fasta] = {}
        for record in AlignIO.read(one_fasta, "fasta"):
            if record.id not in taxons:
                taxons.append(record.id)
            all_seq_data[one_fasta][record.id] = record

    # building dictionnary of the form: dico[one_fasta][one_taxon] = sequence
    concat_data = {}
    for one_fasta in fasta_files:
        for taxon in taxons:
            # since core orthogroups, all taxids should be present everywhere
            assert(taxon in all_seq_data[one_fasta])
            if taxon not in concat_data:
                concat_data[taxon] = all_seq_data[one_fasta][taxon]
            else:
                concat_data[taxon] += all_seq_data[one_fasta][taxon]

    # concatenating the alignments, writing to fasta file
    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")

