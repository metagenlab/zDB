#!/usr/bin/env python

def protein_id2locus_tag(gbk_list, corresp_table, mcl_file):

    '''
    Resolution of a problem encounted with orthomcl: protein from different genomes with the same protein_id. My subsequent scripts/tools need a single accession/protein/genome
    Solution: replace all protein id by locus tag in mcl file
    - all gbk files
    - correspondance table genome accession - mcl genome id (tab separeted file): AE015925	genome_1
    It allow to distinguish proteins from different genomes with the same protein ID
    - mcl file
    '''
    
    mcl_genome2record_accession = {}
    with open(corresp_table, 'r') as f:
        for row in f:
            data = row.rstrip().split('\t')
            mcl_genome2record_accession[data[1]] = data[0]

    print mcl_genome2record_accession
    protein_id2locus_tag = {}

    from Bio import SeqIO
    for one_file in gbk_list:
        handle = open(one_file, "rU")
        for record in SeqIO.parse(handle, "genbank"):
            protein_id2locus_tag[record.name] = {}
            for i in range(0, len(record.features)):
                if record.features[i].type == "CDS":
                    try:
                        len(record.features[i].qualifiers['translation'])
                    except:
                        #print record.features[i]
                        #sys.stderr.write(seq_feature.location.start)
                        #sys.stderr.write("pseudogene?")
                        continue
                    try:
                        protein_id2locus_tag[record.name][record.features[i].qualifiers['protein_id'][0]] = record.features[i].qualifiers['locus_tag'][0]
                    except:
                        protein_id2locus_tag[record.name][record.features[i].qualifiers['locus_tag'][0]] = record.features[i].qualifiers['locus_tag'][0]

    out_file=open("merged_mcl_renamed.txt", 'w')
    with open(mcl_file, 'r') as f:
        for row in f:
            orthogroup_data = row.rstrip().split('\t')
            data_renamed = []
            for i in orthogroup_data:
                protein_data = i.split('|')
                accession = mcl_genome2record_accession[protein_data[0]]
                print accession
                try:
                    locus_tag = protein_id2locus_tag[accession][protein_data[1]]
                except:
                    print protein_id2locus_tag[accession].keys()
                    import sys
                    sys.exit()
                data_renamed.append("%s|%s" % (protein_data[0], locus_tag))
            out_file.write('\t'.join(data_renamed)+'\n')
    out_file.close()
    #print protein_id2locus_tag



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-g",'--input_gbk',type=str,help="input genbank", nargs="+")
    parser.add_argument("-c",'--corresp_table',type=str,help="corresp table: gbk locus 2 genome id used in orthomcl")
    parser.add_argument("-m",'--mcl_file', type=str, help="mcl file")



    args = parser.parse_args()
    protein_id2locus_tag(args.input_gbk, args.corresp_table, args.mcl_file)
