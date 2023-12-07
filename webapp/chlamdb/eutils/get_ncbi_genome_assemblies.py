#!/usr/bin/env python

from Bio import Entrez
Entrez.email = "trestan.pillonel@unil.ch"


def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]


def download_refseq_assemblies(id_list, complete=True):
    import taxid2genomes
    import os

    # handle = Entrez.esearch(db="assembly", term='txid%s[Organism:exp] AND ("reference genome"[filter])' % ncbi_taxon, retmax=100000)

    # record = Entrez.read(handle)
    # id_list = record['IdList']
    # print len(id_list)
    local_dir = os.getcwd()

    for n, id in enumerate(id_list):
        print('download %s out of %s: %s ...' % (n+1, len(id_list), id))
        taxid2genomes.get_complete_genomes_data(id, complete=True)


def write_assembly_table(taxid2assemblies,
                         taxon_id2classification,
                         clade2taxons=False,
                         rank_name=False,
                         output_name='out.txt'):
    o = open(output_name, 'w')
    o.write('target_rank\trank_taxid\trepresentative_species\ttaxid_representative_species\tassembly_name\tRefSeq\tGenbank\tcontig_count\tstatus\trepres\tlength\tclade\tn_subtaxid_available\tno_rank\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\n')

    for clade in taxid2assemblies:
        assembly_list = taxid2assemblies[clade]
        for assembly in assembly_list:
            ncbi_taxon = assembly[0]
            AssemblyName = assembly[1]
            RefSeq = assembly[2]
            Genbank = assembly[3]
            contig_count = assembly[4]
            status = assembly[5]
            repres = assembly[6]
            length = assembly[7]

            try:
                species = taxon_id2classification[ncbi_taxon]['species'][0]
            except Exception:
                species = '-'
            try:
                genus = taxon_id2classification[ncbi_taxon]['genus'][0]
            except Exception:
                genus = '-'
            try:
                family = taxon_id2classification[ncbi_taxon]['family'][0]
            except Exception:
                family = '-'
            try:
                order = taxon_id2classification[ncbi_taxon]['order'][0]
            except Exception:
                order = '-'
            try:
                tclass = taxon_id2classification[ncbi_taxon]['class'][0]
            except Exception:
                tclass = '-'
            try:
                phylum = taxon_id2classification[ncbi_taxon]['phylum'][0]
            except Exception:
                phylum = '-'
            if clade2taxons:
                no_rank = taxon_id2classification[ncbi_taxon]['no rank'][0]
                superkingdom = taxon_id2classification[ncbi_taxon]['superkingdom'][0]
                if rank_name in taxon_id2classification[ncbi_taxon]:
                    rank_taxid = taxon_id2classification[ncbi_taxon][rank_name][1]
                else:
                    rank_taxid = '-'
            else:
                rank_taxid = '-'
                superkingdom = ''
                no_rank = ''

            if clade2taxons:
                n_assemblies_in_clade = len(clade2taxons[clade])
            else:
                n_assemblies_in_clade = 1

            o.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (rank_name,
                                                                                                          rank_taxid,
                                                                                                          species,
                                                                                                          ncbi_taxon,
                                                                                                          AssemblyName,
                                                                                                          Genbank,
                                                                                                          RefSeq,
                                                                                                          contig_count,
                                                                                                          status,
                                                                                                          repres,
                                                                                                          length,
                                                                                                          clade,
                                                                                                          n_assemblies_in_clade,
                                                                                                          no_rank,
                                                                                                          superkingdom,
                                                                                                          phylum,
                                                                                                          tclass,
                                                                                                          order,
                                                                                                          family,
                                                                                                          genus))


def make_random_genome_selection(taxon_id2classification,
                                 rank_name,
                                 n_representative=1,
                                 exclude_unclassified=True,
                                 output_name='out.txt'):
    import taxid2genomes

    '''
    randomly chose n representative taxon for each clade of a defined rank
    i.e: get one representative genome of each Archae Order
    Selection is random but priority is given to
    1. complete genomes
    2. reference genomes
    3. representative genomes
    4. if none of the previsously mentionned label can be found, a random assembly is selected

    :param taxon_id2classification:
    :param rank_name:
    :return:
    '''

    import random
    clade2taxons = {}
    # for each clade of the rank, get list of taxons
    print('number of unique taxids: %s' % (len(taxon_id2classification)))
    unclassified_count = 1
    for taxon in taxon_id2classification:
        # skip unclassified
        if exclude_unclassified:
            if 'unclassified' in taxon_id2classification[taxon]['no rank'][0]:
                print('Incomplete classification, skipping: %s' %
                      taxon_id2classification[taxon]['no rank'][0])
                continue
            elif ' sp. ' in taxon_id2classification[taxon][rank_name][0]:
                print('Incomplete classification, skipping: %s' %
                      taxon_id2classification[taxon][rank_name][0])
                continue

        # get all assemblies associated to each taxon_id
        # unnamed ranks
        if rank_name not in taxon_id2classification[taxon]:
            print('not %s for %s!-----' %
                  (rank_name, taxon_id2classification[taxon]))
            clade2taxons["unkown_%s" %
                         unclassified_count] = taxid2genomes.get_taxi2assembly_accession(taxon)
            unclassified_count += 1

        else:
            # named ranks
            if taxon_id2classification[taxon][rank_name][0] not in clade2taxons:
                clade2taxons[taxon_id2classification[taxon][rank_name]
                             [0]] = taxid2genomes.get_taxi2assembly_accession(taxon)
            else:
                clade2taxons[taxon_id2classification[taxon][rank_name][0]
                             ] += taxid2genomes.get_taxi2assembly_accession(taxon)

    # from assembly list
    # keep first complete genomes

    clade2assemblies = {}
    for clade in clade2taxons:
        # print 'initial', clade2taxons[clade]
        # random.shuffle(clade2taxons[clade])
        # print 'ranomized', clade2taxons[clade]
        # keep.append([clade, clade2taxons[clade][0:n_representative]])
        assembly_list = clade2taxons[clade]
        sorted_assemblies = taxid2genomes.sort_assembly_list(assembly_list)
        clade2assemblies[clade] = sorted_assemblies[0:n_representative]

    print('number after random selection: %s' % len(clade2assemblies.keys()))

    write_assembly_table(clade2assemblies,
                         taxon_id2classification,
                         clade2taxons,
                         rank_name,
                         output_name=output_name)

    return clade2assemblies


def search_assembly_database(ncbi_taxon,
                             complete=True,
                             representative=True,
                             reference=False,
                             exclude_metagenome_derived=True,
                             refseq=True,
                             genbank=True):

    # txid588605[Organism:exp] AND (latest[filter]) AND (all[filter] NOT anomalous[filter]) AND (all[filter] NOT "derived from metagenome"[filter])

    seach_term = 'txid%s[Organism:exp] AND (latest[filter]) AND (all[filter] NOT anomalous[filter]) '

    if refseq and not genbank:
        seach_term += ' AND("latest refseq"[filter]) '
    if genbank and not refseq:
        seach_term += ' AND("latest genbank"[filter]) '
    if genbank and refseq:
        seach_term += ' AND ("latest genbank"[filter] OR "latest refseq"[filter]) '
    if complete:
        seach_term += 'AND ("complete genome"[filter]) '
    if representative and reference:
        seach_term += 'AND ("representative genome"[filter] OR "reference genome"[filter]) '
    if representative and not reference:
        seach_term += 'AND ("representative genome"[filter]) '
    if reference and not representative:
        seach_term += 'AND ("reference genome"[filter]) '
    if exclude_metagenome_derived:
        seach_term += 'AND (all[filter] NOT "derived from metagenome"[filter])'
    print("search term for the NCBI assembly database:", seach_term % ncbi_taxon)
    handle = Entrez.esearch(db="assembly", term=seach_term %
                            ncbi_taxon, retmax=100000)
    # handle = Entrez.esearch(db="assembly", term='txid%s[Organism:exp] AND ("representative genome"[filter] OR "reference genome"[filter])' % ncbi_taxon, retmax=100000)
    # handle = Entrez.esearch(db="assembly", term='txid%s[Organism:exp] AND ("reference genome"[filter])' % ncbi_taxon, retmax=100000)

    record = Entrez.read(handle)
    id_list = record['IdList']
    return id_list


def get_complete_genomes_taxonomy(ncbi_taxon,
                                  complete=True,
                                  representative=True,
                                  reference=False,
                                  exclude_metagenome_derived=True,
                                  refseq=True,
                                  genbank=False):

    import sequence_id2scientific_classification

    assembly_id_list = search_assembly_database(ncbi_taxon,
                                                complete=complete,
                                                representative=representative,
                                                reference=reference,
                                                exclude_metagenome_derived=exclude_metagenome_derived,
                                                refseq=refseq,
                                                genbank=genbank)

    print('Number of assemblies: %s' % len(assembly_id_list))

    id_lists = _chunks(assembly_id_list, 300)
    taxo_data = {}
    for one_list in id_lists:
        handle2 = Entrez.elink(
            dbfrom="assembly", db="taxonomy", id=','.join(one_list))
        record2 = Entrez.read(handle2)

        taxon_list = [i['Id'] for i in record2[0]['LinkSetDb'][0]['Link']]

        taxo_data.update(
            sequence_id2scientific_classification.taxon_id2scientific_classification_and_taxids(taxon_list))
    return taxo_data


def taxon2genome_subset(ncbi_taxon,
                        taxon_rank='order',
                        complete=True,
                        reference=False,
                        representative=True,
                        genbank=False,
                        refseq=True,
                        exclude_metagenome_derived=True,
                        exclude_unclassified=True,
                        get_complete_assembly_list=False,
                        output_name='out.tab'):

    if get_complete_assembly_list:
        import taxid2genomes
        # get complete assembly list
        print_assembly_table = taxid2genomes.get_taxi2assembly_accession(ncbi_taxon,
                                                                         complete=complete,
                                                                         representative=representative,
                                                                         reference=reference,
                                                                         exclude_metagenome_derived=exclude_metagenome_derived,
                                                                         refseq=refseq,
                                                                         genbank=genbank)

        taxonomy_dico = get_complete_genomes_taxonomy(ncbi_taxon,
                                                      complete=complete,
                                                      reference=reference,
                                                      representative=representative,
                                                      exclude_metagenome_derived=exclude_metagenome_derived,
                                                      genbank=genbank,
                                                      refseq=refseq)
        taxid2assemblies = []
        taxid2assemblies[ncbi_taxon] = print_assembly_table
        write_assembly_table(taxid2assemblies,
                             taxonomy_dico,
                             output_name=output_name)

    else:
        # retrieve detailed taxonomy of each assembly
        taxonomy_dico = get_complete_genomes_taxonomy(ncbi_taxon,
                                                      complete=complete,
                                                      reference=reference,
                                                      representative=representative,
                                                      exclude_metagenome_derived=exclude_metagenome_derived,
                                                      genbank=genbank,
                                                      refseq=refseq)

        print('Performing random selection at the %s level...' % taxon_rank)
        dw_taxonomy = make_random_genome_selection(taxonomy_dico,
                                                   taxon_rank,
                                                   exclude_unclassified=exclude_unclassified,
                                                   output_name=output_name)
        # taxon_list = []
        # get all_taxon_list
        # for i in dw_taxonomy:
        #    taxon_list+=i[1]
        # download
        # print ('downloading %s taxon...' % (len(taxon_list)))
        # download_refseq_assembliesdownload_refseq_assemblies(taxon_list)
        # download_refseq_assemblies(taxonomy_dico.keys())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--taxon_id', type=str, help="taxon_idd")
    parser.add_argument("-c", '--complete', action='store_true',
                        help="Filter complete genomes (default: False)")
    parser.add_argument("-rf", '--reference', action='store_true',
                        help="Filter reference genomes (default: False)")
    parser.add_argument("-rp", '--representative', action='store_true',
                        help="Filter representative genomes (default: False)")
    parser.add_argument("-rs", '--refseq', action='store_true',
                        help="Filter refseq assemblies (default: False)")
    parser.add_argument("-gb", '--genbank', action='store_true',
                        help="Filter genbank assemblies (default: False)")
    parser.add_argument("-md", '--exclude_metagenome_derived', action='store_true',
                        help="Exclude metagenome derived assemblies (default: False)")
    parser.add_argument("-e", '--exclude_unclassified', action='store_true',
                        help='Exclude species labelled as "unclassified xxx" or as "xxxx sp. yyy" (e.g. Blautia sp. TF11-31AT) (default=False)')
    parser.add_argument("-d", '--complete_list', action='store_true',
                        help="Report complete list of assemblies (no taxonomy clustering) (default: False) ")
    parser.add_argument("-ra", '--rank', default='species',
                        help="Keep one prepresentative/clade of specified taxonomic <rank>")
    parser.add_argument("-o", '--output', default=None,
                        help="Output name (default: assemblies_<taxid>.tab )")

    args = parser.parse_args()
    # get_complete_genomes_list(args.taxon_id)
    # get_reference_genomes(args.taxon_id)

    if not args.output:
        output_name = 'assemblies_%s.tab' % args.taxon_id
    else:
        output_name = args.output

    taxon2genome_subset(args.taxon_id,
                        args.rank,
                        args.complete,
                        args.reference,
                        args.representative,
                        get_complete_assembly_list=args.complete_list,
                        exclude_unclassified=args.exclude_unclassified,
                        refseq=args.refseq,
                        genbank=args.genbank,
                        exclude_metagenome_derived=args.exclude_metagenome_derived,
                        output_name=output_name)
