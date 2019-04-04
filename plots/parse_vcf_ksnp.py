#!/usr/bin/env python
import heatmap
from biosqldb import manipulate_biosqldb

def parse_vcf(vcf_file_name):
    import vcf
    import numpy as np
    vcf_reader = vcf.VCFReader(open(vcf_file_name, 'rb'))
    vcf_sample_matrix = []
    i = 0
    samples_list = []
    for record in vcf_reader:
        #print i
        i+=1
        #print record
        import os
        samples = [record.POS]
        for call in record.samples:
            #print call, call.gt_type
            #print call.called, call.gt_type, call.gt_bases
            if i == 1:
                samples_list.append(call.sample)
            if call.gt_type == None:
                
                samples.append(100)
            else:
                samples.append(call.gt_type)
        #print samples
        #os.quit()
        vcf_sample_matrix.append(samples)
    return (np.asarray(vcf_sample_matrix, dtype=np.float32), samples_list)


"""
def density_track(vcf_sorted_matrix, genome_size_bp, range = 1000):
    start = 1
    stop = range
    while stop < genome_size_bp:
        for column in range(0, len(vcf_sorted_matrix[1, :])):
            for one_genome
"""

def vcf2heatmap(vcf_file, biodb ="saureus_01_15", output="heat_ksnp_SaC.pdf"):
    my_array, sample_names = parse_vcf(vcf_file)
    server, db = manipulate_biosqldb.load_db(biodb)
    accession2description = manipulate_biosqldb.accession2description(server, biodb)

    for i in range(0, len(sample_names)):
        try:
            sample_names[i] = accession2description[sample_names[i]]
        except:
            pass

    M = my_array[:, 0:]
    #print M[0:10,:]
    # tri en fonction de la localisation
    #M = M[M[:, 0].argsort()]

    core = 0
    core_locations = []
    non_core = 0
    for index, value in enumerate(M[:,1:].sum(axis=1)):
        if value > 50:
            core+=1
            core_locations.append(int(M[index,0]))
        else:
            non_core+=1
    print "core", core, len(core_locations)
    print "non core", non_core
    print core/float((core+non_core))

    #M = heatmap.randomize_table(M[:,1:])

    #heatmap.heatmap_ksnp(M, format='png', output=output, breaks="-0.5, 0.5, 5, 99", rows=None, columns = sample_names,orderRows = False)
    return core_locations

def parse_vcf_get_pairwise_n_snp(vcf_file_name, core_list = False):
    import vcf
    print core_list
    import numpy as np
    vcf_reader = vcf.VCFReader(open(vcf_file_name, 'rb'))
    vcf_sample_matrix = []
    i = 0
    samples_list = []
    position2index = {}
    nsnp = {}
    for record in vcf_reader:
        print i, type(record.POS), record
        i+=1
        position2index[record.POS] = record.ID
        samples = [record.POS]
        for call in record.samples:

            #print call, call.gt_type
            #print call.called, call.gt_type, call.gt_bases

            if i == 1:
                samples_list.append(call.sample)
                nsnp[call.sample] = {}
                nsnp[call.sample]["identical"] = 0
                nsnp[call.sample]["diff"] = 0
                nsnp[call.sample]["NA"] = 0
                nsnp[call.sample]["snp"] = []
            if core_list:
                if call.gt_type == None and record.POS in core_list:
                    samples.append(100)
                    nsnp[call.sample]["NA"] += 1
                elif call.gt_type == 0 and record.POS in core_list:
                    nsnp[call.sample]["identical"] += 1
                elif record.POS in core_list:
                    nsnp[call.sample]["diff"] += 1
                    nsnp[call.sample]["snp"].append(record.POS)
                else:
                    continue
            else:
                if call.gt_type == None:
                    samples.append(100)
                    nsnp[call.sample]["NA"] += 1
                elif call.gt_type == 0:
                    nsnp[call.sample]["identical"] += 1
                else:
                    nsnp[call.sample]["diff"] += 1
                    nsnp[call.sample]["snp"].append(record.POS)



                #import sys
        #print samples
        #os.quit()
        vcf_sample_matrix.append(samples)
    for i in nsnp:
        print i, "NA:", nsnp[i]["NA"], "idem:", nsnp[i]["identical"],"diff:", nsnp[i]["diff"], "TOTAL:", nsnp[i]["NA"] + nsnp[i]["identical"] + nsnp[i]["diff"]
    #return (np.asarray(vcf_sample_matrix, dtype=np.float32), samples_list)


