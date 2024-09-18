import os
import re
from collections import defaultdict, namedtuple

import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import CheckSum


def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


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
            if len(groups) > 1:
                new_fasta = [sequence_data[i] for i in groups]
                out_path = "%s.faa" % group_name
                with open(out_path, "w") as out_handle:
                    SeqIO.write(new_fasta, out_handle, "fasta")


class MissingReferenceDatabase(Exception):
    pass


def check_reference_databases(params):
    db_dir = os.path.split(params.get("base_db"))[-1]
    if params.get("cog"):
        if not os.path.isdir(os.path.join(db_dir, "cog")):
            raise MissingReferenceDatabase(
                'COG database could not be found. '
                'Please set it up with "zdb setup --cog"')
    if params.get("ko"):
        if not os.path.isdir(os.path.join(db_dir, "kegg")):
            raise MissingReferenceDatabase(
                'KEGG database could not be found. '
                'Please set it up with "zdb setup --ko"')
    if params.get("pfam"):
        if not os.path.isdir(os.path.join(db_dir, "pfam")):
            raise MissingReferenceDatabase(
                'Pfam database could not be found. '
                'Please set it up with "zdb setup --pfam"')
    if params.get("swissprot"):
        if not os.path.isdir(os.path.join(db_dir, "uniprot", "swissprot")):
            raise MissingReferenceDatabase(
                'Swissprot database could not be found. '
                'Please set it up with "zdb setup --swissprot"')
    if params.get("vfdb"):
        if not os.path.isdir(os.path.join(db_dir, "vfdb")):
            raise MissingReferenceDatabase(
                'VFDB database could not be found. '
                'Please set it up with "zdb setup --vfdb"')


class InvalidInput(Exception):
    pass


class InputHandler():
    """
    This class is in charge of parsing the input csv, checking
    and if necessary revising the GBK files.
    """

    CsvEntry = namedtuple("CsvEntry", "name file groups")
    allowed_headers_expr = re.compile(r"name|file|group:[\w\s]*")
    group_header_expr = re.compile(r"group:([\w\s]*)")

    def __init__(self, csv_file):
        self.csv_entries, self.group_names = self.parse_csv(csv_file)

        # NOTE: biosql uses source to assign taxid. One must ensure
        # that the source are unique to each gbk file to avoid conflicts
        # with taxids.
        self.organisms = defaultdict(lambda: 0)
        self.locuses = defaultdict(lambda: 0)
        self.accessions = defaultdict(lambda: 0)

        # contig name also need to be unique as they are used by blast
        self.contigs = defaultdict(lambda: 0)

    def gen_new_locus_tag(self, hsh_prev_values):
        curr_value = len(hsh_prev_values)
        prefix = "ZDB_"
        locus_tag = f"{prefix}{curr_value}"
        while locus_tag in hsh_prev_values:
            curr_value += 1
            locus_tag = f"{prefix}{curr_value}"
        hsh_prev_values[locus_tag] += 1
        return locus_tag

    def gen_new_organism(self, organism):
        cnt = 1
        prefix = organism + "_"
        new_name = f"{prefix}{cnt}"
        while new_name in self.organisms:
            cnt += 1
            new_name = f"{prefix}{cnt}"
        self.organisms[new_name] += 1
        self.organisms[organism] -= 1
        return new_name

    @classmethod
    def is_header_allowed(cls, header):
        return bool(cls.allowed_headers_expr.fullmatch(header))

    @classmethod
    def header_to_group_label(cls, header):
        return cls.group_header_expr.fullmatch(header).groups()[0]

    @classmethod
    def header_is_group(cls, header):
        return bool(cls.group_header_expr.fullmatch(header))

    @staticmethod
    def is_in_group(entry, groupname):
        value = entry.get(groupname)
        if value.lower() in ["no", "n", "false", "0", ""]:
            return False
        if value.lower() in ["yes", "y", "true", "1", "x"]:
            return True
        raise InvalidInput(
            f'Invalid entry "{value}" in group column. Should be one of '
            '"yes", "true", "x", "1", "no", "false", "0" or empty')

    @classmethod
    def parse_csv(cls, csv_file):
        csv = pd.read_csv(csv_file, dtype=str,
                          skipinitialspace=True).fillna('')
        entries = []
        names = set()

        group_names = {}
        for colname in csv.columns:
            if not cls.is_header_allowed(colname):
                raise InvalidInput(
                    f'Invalid column header "{colname}" in input file.')
            if cls.header_is_group(colname):
                group_names[colname] = cls.header_to_group_label(colname)

        filenames = set()
        entries = []
        for index, entry in csv.iterrows():
            name = entry.get("name", None) or None
            if name:
                if name in names:
                    raise InvalidInput(
                        f'Name "{name}" appears twice in the input file.')
                names.add(name)

            # only get the filename, as nextflow will symlink it
            # in the current work directory
            filename = os.path.basename(entry.file)
            if not os.path.isfile(filename):
                raise InvalidInput(
                    f'File "{filename}" cannot be accessed. Please check '
                    f'that {entry.file} exists and is accessible.')

            if filename in filenames:
                raise InvalidInput(
                    f'File "{filename}" appears twice in the input file.')
            filenames.add(filename)

            groups = []
            for group_name, group_label in group_names.items():
                if cls.is_in_group(entry, group_name):
                    groups.append(group_label)
            entries.append(cls.CsvEntry(name, filename, groups))

        if len(entries) == 0:
            raise Exception(csv_file + " is empty")

        return entries, group_names

    def check_and_revise_gbks(self):
        """
        This method ensures that we do not run into name conflicts in zDB.
        For this we check that the taxid, the contig name, accession
        and locus_tags are unique and if not we assign them new names.

        As BioSQL uses the organism/source entries of the records to assign
        a taxid when it is not allowed to access the ncbi online,
        we need to ensure that the same organism name is not used
        twice in the genomes to prevent bioentries from different genbank
        files to be entered under the same taxon_id (as zdb uses the taxon_id
        to differentiate between genomes).
        """

        gbk_passed = []
        gbk_to_revise = []

        custom_names = {}

        for entry in self.csv_entries:
            gbk_file = entry.file
            curr_organism = None
            n_cds = 0
            needs_revision = False
            for record in SeqIO.parse(gbk_file, "genbank"):
                sci_name = record.annotations.get("organism", None)
                common_name = record.annotations.get("source", None)

                if entry.name is not None:
                    custom_names[entry.file] = entry.name
                    needs_revision = True
                    sci_name = entry.name
                    common_name = entry.name

                if sci_name is None:
                    raise Exception(f"No scientific name for record {record.id} "
                                    f"in {gbk_file}.")

                if record.name in self.contigs:
                    needs_revision = True
                self.contigs[record.name] += 1

                if "accessions" not in record.annotations:
                    needs_revision = True
                else:
                    acc = record.annotations["accessions"][0]
                    if acc in self.accessions:
                        needs_revision = True
                    self.accessions[acc] += 1

                if curr_organism is None:
                    curr_organism = sci_name
                    self.organisms[sci_name] += 1
                    if self.organisms[sci_name] > 1:
                        needs_revision = True
                elif curr_organism != sci_name:
                    raise Exception(f"Two different organisms in {gbk_file}: {curr_organism}/{sci_name}")

                # necessary, as BioSQL will use whatever matches in the
                # common or scientific names to assign taxid. So if the common name
                # of two assemblies is foo, but the scientific are different, BioSQL
                # will still assign them the same taxid.
                if common_name != sci_name:
                    needs_revision = True

                for feature in record.features:
                    if feature.type == "CDS":
                        n_cds += 1
                    elif feature.type not in ["tmRNA", "rRNA", "ncRNA", "tRNA"]:
                        continue

                    if "locus_tag" not in feature.qualifiers:
                        needs_revision = True
                        continue

                    locus_tag = feature.qualifiers["locus_tag"][0]
                    if locus_tag in self.locuses:
                        needs_revision = True
                    self.locuses[locus_tag] += 1

            if n_cds == 0:
                raise Exception(
                    f"No CDS in {gbk_file}, has it been correctly annotated?")

            if needs_revision:
                gbk_to_revise.append(gbk_file)
            else:
                gbk_passed.append(gbk_file)

        for filename, name in custom_names.items():
            if self.organisms[name] > 1:
                raise Exception(
                    f"The custom name {name} is already used in another file")

        # at one point, will have to rewrite this to avoid
        # re-parsing the genbank files that failed the check
        for failed_gbk in gbk_to_revise:
            records = []
            sci_name = custom_names.get(failed_gbk, None)
            for record in SeqIO.parse(failed_gbk, "genbank"):
                if sci_name is None:
                    sci_name = record.annotations["organism"]
                    if self.organisms[sci_name] > 1:
                        sci_name = self.gen_new_organism(sci_name)
                record.annotations["source"] = sci_name
                record.annotations["organism"] = sci_name

                if self.contigs[record.name] > 1:
                    record.name = self.gen_new_locus_tag(self.contigs)
                    self.contigs[record.name] -= 1

                # NOTE: for a reason I do not understand, Biopython will
                # store the accession into the "accessions" entry of a record when parsing
                # a genbank file, but won't recognize the "accessions" entry
                # when writing the same record. It instead recognizes the "accession" entry.
                if "accessions" not in record.annotations:
                    record.annotations["accession"] = [
                        self.gen_new_locus_tag(self.accessions)]
                elif self.accessions[record.annotations["accessions"][0]] > 1:
                    new_acc = self.gen_new_locus_tag(self.accessions)
                    if "accession" not in record.annotations:
                        record.annotations["accession"] = [new_acc]
                    else:
                        record.annotations["accession"][0] = new_acc
                    self.accessions[new_acc] -= 1

                for feature in record.features:
                    if feature.type == "source":
                        if "db_xref" in feature.qualifiers:
                            # need to remove this as BioSQL uses it to assign
                            # taxids and this messes with the rest of the DB
                            # XXX Ugly hack, will have to go
                            del feature.qualifiers["db_xref"]
                    if feature.type not in ["CDS", "tmRNA", "rRNA", "ncRNA", "tRNA"]:
                        continue
                    curr_locus = feature.qualifiers.get("locus_tag", None)
                    if curr_locus is None:
                        feature.qualifiers["locus_tag"] = self.gen_new_locus_tag(
                            self.locuses)
                    elif self.locuses[curr_locus[0]] > 1:
                        feature.qualifiers["locus_tag"] = self.gen_new_locus_tag(
                            self.locuses)
                        self.locuses[curr_locus[0]] -= 1
                records.append(record)
            SeqIO.write(records, "filtered/" + failed_gbk, "genbank")

        for passed in gbk_passed:
            os.symlink(os.readlink(passed), "filtered/" + passed)


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
                    raise Exception(
                        f"Feature without locus tag in record {record.name}")

                locus_tag = feature.qualifiers["locus_tag"][0]
                if output_fmt == "faa":
                    if "translation" not in feature.qualifiers:
                        continue
                    data = feature.qualifiers["translation"][0]
                elif output_fmt == "fna":
                    data = feature.location.extract(record).seq
                else:
                    raise Exception(f"Unsupported option: {output_fmt}, must be either faa or fna")

                edited_records.write(">%s %s\n%s\n" % (locus_tag,
                                                       record.description, data))
    edited_records.close()


def convert_gbk_to_fna(gbf_file, fna_contigs):
    # from BioPython, object of class SeqRecord
    records = SeqIO.parse(gbf_file, 'genbank')
    edited_records = open(fna_contigs, 'w')

    for record in records:
        edited_records.write(">%s %s\n%s\n" %
                             (record.name, record.description, record.seq))
    edited_records.close()


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
    nr_fasta.close()
    nr_mapping.close()


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
            group_id = groups[0][:-1]  # remove lagging ':'
            groups = groups[1:]
            orthogroup2locus_list[group_id] = groups

    locus2genome = {}
    for fasta in fasta_list:
        tokens = os.path.basename(fasta).split('.')
        # the filename may contain additional "." that are to be included
        genome = ".".join(tokens[:-1])
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome

    n_genomes = len(set(locus2genome.values()))
    hsh_ogs = defaultdict(lambda: defaultdict(lambda: 0))
    for group_id, loci_list in orthogroup2locus_list.items():
        # we know there will be paralogs in this case
        if len(loci_list) > n_genomes:
            continue

        # cannot possibly be a core orthogroup
        if len(loci_list) < n_genomes - n_missing:
            continue

        for locus in loci_list:
            genome = locus2genome[locus]
            hsh_to_count = hsh_ogs[genome]
            hsh_to_count[group_id] += 1

    df = pd.DataFrame(hsh_ogs).fillna(0)
    df_sum = df[df < 2].sum(axis=1)
    core_groups = df_sum[df_sum >= n_genomes - n_missing].index.tolist()
    return core_groups, orthogroup2locus_list, locus2genome


def get_core_orthogroups(genomes_list, int_core_missing):
    core_groups, orthogroup2locus_list, locus2genome = orthofinder2core_groups(genomes_list,
                                                                               'Orthogroups.txt', int_core_missing)

    og_resume_file = open("orthogroups_summary_info.tsv", "w")
    core_groups_set = set(core_groups)
    for og, locus_list in orthogroup2locus_list.items():
        og_size = len(locus_list)
        is_core = 1 if og in core_groups_set else 0
        genome_set = set()
        for locus in locus_list:
            genome_set.add(locus2genome[locus])
        num_genomes = len(genome_set)
        print(og, is_core, og_size, num_genomes, sep="\t", file=og_resume_file)
    og_resume_file.close()

    if len(core_groups) <= 0:
        raise Exception(
            "No core orthogroups, maybe try to rerun the pipeline with num_missing to a higher value")

    for group_id in core_groups:
        sequence_data = SeqIO.to_dict(
            SeqIO.parse(group_id + "_mafft.faa", "fasta"))
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
        for j in range(i + 1, len(alignment)):
            alignment_1 = alignment[i]
            alignment_2 = alignment[j]
            locus_tag_1 = alignment_1.name
            locus_tag_2 = alignment_2.name
            alignment_length = 0
            identical = 0
            for ch1, ch2 in zip(alignment_1, alignment_2):
                if ch1 == "-" or ch2 == "-":
                    continue
                if ch1 == ch2:
                    identical += 1
                alignment_length += 1

            if alignment_length > 0:
                per_identity = 100 * (identical / float(alignment_length))
            else:
                per_identity = 0
            values.append((locus_tag_1, locus_tag_2,
                          per_identity, alignment_length))
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
                seq = "-" * fasta_to_algn_length[one_fasta]
            else:
                seq = all_seq_data[one_fasta][taxon]
            concat_data[taxon] += seq

    MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
    with open(out_name, "w") as handle:
        AlignIO.write(MSA, handle, "fasta")
