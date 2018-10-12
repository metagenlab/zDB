#!/usr/bin/env python

from Bio import Entrez
Entrez.email = "trestan.pillonel@unil.ch"

def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def download_refseq_assemblies(id_list, complete=True):
    import taxid2genomes
    import os

    #handle = Entrez.esearch(db="assembly", term='txid%s[Organism:exp] AND ("reference genome"[filter])' % ncbi_taxon, retmax=100000)

    #record = Entrez.read(handle)
    #id_list = record['IdList']
    #print len(id_list)
    local_dir = os.getcwd()

    for n, id in enumerate(id_list):
        print ('download %s out of %s: %s ...' % (n+1, len(id_list), id))
        taxid2genomes.get_complete_genomes_data(id, complete=True)

def make_random_genome_selection(taxon_id2classification,
                                 rank_name,
                                 n_representative=1,
                                 exclude_unclassified=True):
    '''
    randomly chose n representative taxon for each clade of a defined rank
    i.e: get one representative genome of each Archae Order

    :param taxon_id2classification:
    :param rank_name:
    :return:
    '''

    import random
    clade2taxons = {}
    # for each clade of the rank, get list of taxons
    print('numnber of taxons: %s' % (len(taxon_id2classification)))
    unclassified_count=1
    for taxon in taxon_id2classification:
        # skip unclassified
        if exclude_unclassified:
            if 'unclassified' in taxon_id2classification[taxon]['no rank'][0]:
                print('unclassified, skipping!')
                continue
        # unnamed ranks
        if not rank_name in taxon_id2classification[taxon]:
            print ('not %s for %s!-----' % (rank_name, taxon_id2classification[taxon]))
            clade2taxons["unkown_%s" % unclassified_count] = [taxon]
            unclassified_count+=1

        else:
            ## named ranks
            if taxon_id2classification[taxon][rank_name][0] not in clade2taxons:
                clade2taxons[taxon_id2classification[taxon][rank_name][0]] = [taxon]
            else:
                clade2taxons[taxon_id2classification[taxon][rank_name][0]].append(taxon)
    keep = []
    for clade in clade2taxons:
        #print 'initial', clade2taxons[clade]
        random.shuffle(clade2taxons[clade])
        #print 'ranomized', clade2taxons[clade]
        keep.append([clade, clade2taxons[clade][0:n_representative]])
    print('number after random selection: %s' % len(keep))
    print('clade\tclade_taxid\tn_genomes_available\ttaxid_selected\tno_rank\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies')
    for clade, taxons in keep:

        for taxon in taxons:
            try:
                species = taxon_id2classification[taxon]['species'][0]
            except:
                species = '-'
            try:
                genus = taxon_id2classification[taxon]['genus'][0]
            except:
                genus = '-'
            try:
                family = taxon_id2classification[taxon]['family'][0]
            except:
                family = '-'
            try:
                order = taxon_id2classification[taxon]['order'][0]
            except:
                order = '-'
            try:
                tclass = taxon_id2classification[taxon]['class'][0]
            except:
                tclass = '-'
            try:
                phylum = taxon_id2classification[taxon]['phylum'][0]
            except:
                phylum = '-'
            if clade =='[Clostridium] clostridioforme':
                print (clade2taxons[clade])
            if rank_name in taxon_id2classification[taxon]:
                rank_taxid = taxon_id2classification[taxon][rank_name][1]
            else:
                rank_taxid = '-'
            print ('%s (%s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                                                      clade,
                                                      rank_name,
                                                      rank_taxid,
                                                      len(clade2taxons[clade]),
                                                      taxon,
                                                      taxon_id2classification[taxon]['no rank'][0],
                                                      taxon_id2classification[taxon]['superkingdom'][0],
                                                      phylum,
                                                      tclass,
                                                      order,
                                                      family,
                                                      genus,
                                                      species))
    return keep


def get_complete_genomes_taxonmy(ncbi_taxon,
                                 complete=True,
                                 representative=True,
                                 reference=False,
                                 exclude_metagenome_derived=True):

    import sequence_id2scientific_classification

    # txid588605[Organism:exp] AND (latest[filter]) AND (all[filter] NOT anomalous[filter]) AND (all[filter] NOT "derived from metagenome"[filter])

    seach_term = 'txid%s[Organism:exp] AND (latest[filter]) AND (all[filter] NOT anomalous[filter]) AND ("latest refseq"[filter])'
    if complete:
        seach_term+='AND ("complete genome"[filter]) '
    if representative and reference:
        seach_term+='AND ("representative genome"[filter] OR "reference genome"[filter]) '
    if representative and not reference:
        seach_term+='AND ("representative genome"[filter]) '
    if reference and not representative:
        seach_term+='AND ("reference genome"[filter]) '
    if exclude_metagenome_derived:
        seach_term+='AND (all[filter] NOT "derived from metagenome"[filter]) '
    print (seach_term)
    handle = Entrez.esearch(db="assembly", term=seach_term % ncbi_taxon, retmax=100000)
    #handle = Entrez.esearch(db="assembly", term='txid%s[Organism:exp] AND ("representative genome"[filter] OR "reference genome"[filter])' % ncbi_taxon, retmax=100000)
    #handle = Entrez.esearch(db="assembly", term='txid%s[Organism:exp] AND ("reference genome"[filter])' % ncbi_taxon, retmax=100000)

    record = Entrez.read(handle)
    id_list = record['IdList']
    #print 'id_list', id_list[0:10]

    id_lists = _chunks(id_list, 300)
    taxo_data = {}
    for one_list in id_lists:
        handle2 = Entrez.elink(dbfrom="assembly", db="taxonomy", id=','.join(one_list))
        record2 = Entrez.read(handle2)

        taxon_list = [i['Id'] for i in record2[0]['LinkSetDb'][0]['Link']]

        taxo_data.update(sequence_id2scientific_classification.taxon_id2scientific_classification_and_taxids(taxon_list))
    return taxo_data

def taxon2genome_subset(ncbi_taxon,
                        taxon_rank='order',
                        complete=True,
                        reference=False,
                        representative=True,
                        show_all=True,
                        print_only=True,
                        download_all=False,
                        exclude_metagenome_derived=True,
                        exclude_unclassified=True):

    taxonomy_dico = get_complete_genomes_taxonmy(ncbi_taxon,
                                                 complete=complete,
                                                 reference=reference,
                                                 representative=representative,
                                                 exclude_metagenome_derived=exclude_metagenome_derived)

    if print_only:
        print('print only')
        if not show_all:
            print('not show all')
            make_random_genome_selection(taxonomy_dico, taxon_rank)
        else:
            print('taxo_data')
            print_taxo_data(taxonomy_dico,
                            exclude_unclassified=exclude_unclassified)
    else:
        print('not print only')
        if not download_all:
            dw_taxonomy = make_random_genome_selection(taxonomy_dico,
                                                       taxon_rank,
                                                       exclude_unclassified=exclude_unclassified)
            taxon_list = []
            # get all_taxon_list
            for i in dw_taxonomy:
                taxon_list+=i[1]
            # download
            print ('downloading %s taxon...' % (len(taxon_list)))
            download_refseq_assembliesdownload_refseq_assemblies(taxon_list)
        else:
            print ('downloading %s taxon...' % (len(taxonomy_dico.keys())))
            download_refseq_assemblies(taxonomy_dico.keys())

def print_taxo_data(taxo_data, exclude_unclassified=True):
        # {'superkingdom': 'Archaea', 'no rank': 'cellular organisms', 'family': 'Natrialbaceae', 'order': 'Natrialbales', 'phylum': 'Euryarchaeota', 'species': 'Halobiforma lacisalsi', 'genus': 'Halobiforma', 'class': 'Halobacteria'}
        for taxon in taxo_data:
            try:
                species = taxo_data[taxon]['species']
            except:
                species = '-'
            try:
                genus = taxo_data[taxon]['genus']
            except:
                genus = '-'
            try:
                family = taxo_data[taxon]['family']
            except:
                family = '-'
            try:
                order = taxo_data[taxon]['order']
            except:
                order = '-'
            try:
                tclass = taxo_data[taxon]['class']
            except:
                tclass = '-'
            try:
                phylum = taxo_data[taxon]['phylum']
            except:
                phylum = '-'

            #print len(taxo_data[i]), taxo_data[i]
            if exclude_unclassified:
                if 'unclassified' in taxo_data[taxon]['no rank']:
                    continue
            print ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (taxon,
                                                      taxo_data[taxon]['no rank'],
                                                      taxo_data[taxon]['superkingdom'],
                                                      phylum,
                                                      tclass,
                                                      order,
                                                      family,
                                                      genus,
                                                      species))




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--taxon_id',type=str,help="taxon_idd")

    parser.add_argument("-c", '--complete', action='store_true', help="keep complete genomes")
    parser.add_argument("-rf", '--reference', action='store_true', help="keep reference genomes")
    parser.add_argument("-rp", '--representative', action='store_true', help="keep representative genomes")
    parser.add_argument("-d", '--download_all', action='store_true', help="download all hit assemblies")
    parser.add_argument("-r", '--download_representative', action='store_true', help="download representative hit assemblies")
    parser.add_argument("-ra", '--rank', default='order', help="download one prepresentative/clade of specified taxonomic <rank>")
    parser.add_argument("-p", '--print_only', action='store_false', help="only print data, no download (default=True)")

    args = parser.parse_args()
    #get_complete_genomes_list(args.taxon_id)
    #get_reference_genomes(args.taxon_id)
    taxon2genome_subset(args.taxon_id,
                        args.rank,
                        args.complete,
                        args.reference,
                        args.representative,
                        print_only=args.print_only,
                        download_all=args.download_all)
