#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import concat_gbk


def filter_plasmid(record_list):
    '''
    Return 2 lists: one with plasmid record(s), one without plasmids
    :param record_list:
    :return:
    '''

    plasmid_record_list = []
    chromosome_record_list = []

    for record in record_list:
        # plasmid.annotations['organism']
        if record.features[0].type == 'source':
            if 'plasmid' in record.description or "plasmid" in record.features[0].qualifiers:
                plasmid_record_list.append(record)
            else:
                chromosome_record_list.append(record)
        else:
            if 'plasmid' in record.description:
                plasmid_record_list.append(record)
            else:
                chromosome_record_list.append(record)
    return (chromosome_record_list, plasmid_record_list)


def count_missing_locus_tags(gbk_record):
    count_CDS = 0
    count_no_locus = 0
    for feature in gbk_record.features:
        if feature.type == 'CDS':
            count_CDS += 1
            try:
                test = feature.qualifiers['locus_tag']
            except:
                # missing locus case
                count_no_locus += 1
        pass
    return count_no_locus, count_CDS


def is_annotated(gbk_record):
    if len(gbk_record.features) == 1 and gbk_record.features[0].type == 'source':
        return False
    else:
        return True


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


def check_gbk(gbff_files,
              minimal_contig_length=1000):

    from Bio import SeqIO
    import gzip

    reannotation_list = []

    # count the number of identical source names
    source2count = {}
    accession2count = {}
    for genome_number, gbff_file in enumerate(gbff_files):

        records = list(SeqIO.parse(gzip.open(gbff_file, "rt"), "genbank"))

        for record in records:
            n_missing, total = count_missing_locus_tags(record)
            if n_missing > 0:
                print ('Warrning: %s/%s missing locus tag for record %s' % (n_missing, total, record.name))

        chromosome, plasmids = filter_plasmid(records)

        cleaned_records = []
        plasmid_reannot = False
        chromosome_reannot = False

        if len(plasmids) > 0:

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
            Contigs contatenation with 200 N between each (labelled as assembly_gap feature)
            '''

            ########## chromosome ###########
            if chromosome[0].seq == 'N'*len(chromosome[0].seq):
                print('Warning: No sequences for %s, skipping' % gbff_file)
                continue
            annot_check = is_annotated(chromosome[0])
            if annot_check:
                if len(chromosome) > 1:
                    merged_record = concat_gbk.merge_gbk(chromosome)
                else:
                    merged_record = chromosome[0]
                print(merged_record.description )
                merged_record.description = clean_description(merged_record.description)
                print(merged_record.description)
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



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--gbff_files', type=str, help="input gbff file(s)", nargs="+")
    parser.add_argument("-l", '--minimal_contig_length', type=int, help="input gbff file(s)")

    args = parser.parse_args()
    check_gbk(args.gbff_files,
              args.minimal_contig_length)
