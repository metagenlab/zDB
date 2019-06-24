#! /usr/bin/env python

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqRecord import SeqRecord
import json
import pickle
import os
from time import time
from chlamdb.biosqldb import manipulate_biosqldb
from multiprocessing import Process, Queue, JoinableQueue
import math
from Bio.Seq import Seq

def get_feature_neighborhood(feature_start, feature_end, contig_or_genome_record, neighborhood_size_bp, record_name):

    start = feature_start - neighborhood_size_bp
    end = feature_end + neighborhood_size_bp
    if contig_or_genome_record.features[0].location.start > start:
            start = 0
    if contig_or_genome_record.features[0].location.end < end:
            end = contig_or_genome_record.features[0].location.end
    record = contig_or_genome_record[int(start):int(end)]
    return record

def plot_multiple_regions_crosslink2(target_protein_list, region_record_list, plasmid_list, out_name):
    gd_diagram = GenomeDiagram.Diagram("geomic_region")
    feature_sets = []
    max_len = 0
    records = dict((rec.name, rec) for rec in region_record_list)
    n_records = len(region_record_list)
    
    record_length = [len(record) for record in region_record_list] 

    for i, record in enumerate(region_record_list):
        max_len = max(max_len, len(record))
        print "i", i
        #Allocate tracks 3 (top), 1 (bottom) for region 1 and 2
        #(empty tracks 2 useful white space to emphasise the cross links
        #and also serve to make the tracks vertically more compressed)
        gd_track_for_features = gd_diagram.new_track((2*n_records-1)-2*i, name=record.name, greytrack=True, height=0.5, start=0, end=len(record))
        if record.name not in feature_sets:
                feature_sets.append(gd_track_for_features.new_set())
        else:
           print "already in feature_sets!"
           print record
           quit
        
    for x in range(0,len(region_record_list)-1):
        print "x", x
        features_X = region_record_list[x].features
        features_Y = region_record_list[x+1].features
        set_X = feature_sets[x]
        set_Y = feature_sets[x+1]
        for feature_1 in features_X:
            if feature_1.type != "CDS":
                continue
            for feature_2 in features_Y:
                if feature_2.type != "CDS":
                    continue
                try:

                    group1 = feature_1.qualifiers["orthogroup"][0]
                    group2 = feature_2.qualifiers["orthogroup"][0]

                except:
                        group1 = "one_singleton"
                        group2 = "two_singleton"
                print "group1, group2", group1, group2
                if group1 == group2:
                    border = colors.lightgrey
                    color = colors.lightgrey
                    F_x = set_X.add_feature(SeqFeature(FeatureLocation(feature_1.location.start, feature_1.location.end, strand=0)),
                                    color=color, border=border)
                    F_y = set_Y.add_feature(SeqFeature(FeatureLocation(feature_2.location.start, feature_2.location.end, strand=0)),
                                    color=color, border=border)
                    gd_diagram.cross_track_links.append(CrossLink(F_x, F_y, color, border))

    #for x in range(0,len(region_record_list)-1):
    x = 0
    for n, record in enumerate(region_record_list):
        gd_feature_set = feature_sets[n]
        i = 0

        if plasmid_list[x]:
            #print "PLASMID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            color1 = colors.HexColor('#2837B7')
            color2 = colors.blue
        else:
            color1 = colors.HexColor('#40F13A')
            color2 = colors.HexColor('#0F600C')



        for feature in record.features:
            if feature.type != "CDS":
                continue
            try:
                a = feature.qualifiers["locus_tag"]
            except:
                # cas des pseudogenes qui sont des CDS mais n'ont pas de protein ID
                continue

            if len(gd_feature_set) % 2 == 0:
                color = color1
            else:
                color = color2

            #try:
            #    try:
            #            group = protein_id2group[feature.qualifiers["protein_id"][0]]
            #    except:
            #            group = protein_id2group[feature.qualifiers["protein_id"][1]]
            #except:
            #    # no group attributed: singleton => special color
            #    color = colors.HexColor('#E104C0')


            for target_protein in target_protein_list:
                    if target_protein in feature.qualifiers["locus_tag"]:
                        print "target prot!"
                        color = colors.red

            gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_position="middle",label_strand=1, label_size=12, label_angle=45)
            i += 1
        x += 1

    print "max", max_len
    hauteur = 250*len(region_record_list)
    largeur = max(record_length)/30
    print "hauteur", hauteur
    print "largeur", largeur
    #gd_diagram.set_page_size(, orientation)
    if hauteur > largeur:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='portrait', fragments=1,start=0, end=max_len)
    else:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='landscape', fragments=1,start=0, end=max_len)
    print "writing diagram", out_name
    
    gd_diagram.write(out_name, "SVG")

def plot_multiple_regions_crosslink(target_protein_list, region_record_list, plasmid_list, out_name):
    gd_diagram = GenomeDiagram.Diagram("geomic_region")
    feature_sets = []
    max_len = 0
    records = dict((rec.name, rec) for rec in region_record_list)
    n_records = len(region_record_list)
    
    record_length = [len(record) for record in region_record_list] 

    for i, record in enumerate(region_record_list):
        max_len = max(max_len, len(record))
        print "i", i
        #Allocate tracks 3 (top), 1 (bottom) for region 1 and 2
        #(empty tracks 2 useful white space to emphasise the cross links
        #and also serve to make the tracks vertically more compressed)
        gd_track_for_features = gd_diagram.new_track((1*n_records-1)-1*i, name=record.name, greytrack=True, height=0.5, start=0, end=len(record))
        if record.name not in feature_sets:
                feature_sets.append(gd_track_for_features.new_set())
        else:
           print "already in feature_sets!"
           print record
           quit
        
    for x in range(0,len(region_record_list)-1):
        print "x", x
        features_X = region_record_list[x].features
        features_Y = region_record_list[x+1].features
        set_X = feature_sets[x]
        set_Y = feature_sets[x+1]
        for feature_1 in features_X:
            if feature_1.type != "CDS":
                continue
            for feature_2 in features_Y:
                if feature_2.type != "CDS":
                    continue
                try:

                    group1 = feature_1.qualifiers["orthogroup"][0]
                    group2 = feature_2.qualifiers["orthogroup"][0]

                except:
                        group1 = "one_singleton"
                        group2 = "two_singleton"
                print "group1, group2", group1, group2
                if group1 == group2:
                    border = colors.lightgrey
                    color = colors.lightgrey
                    F_x = set_X.add_feature(SeqFeature(FeatureLocation(feature_1.location.start, feature_1.location.end, strand=0)),
                                    color=color, border=border)
                    F_y = set_Y.add_feature(SeqFeature(FeatureLocation(feature_2.location.start, feature_2.location.end, strand=0)),
                                    color=color, border=border)
                    gd_diagram.cross_track_links.append(CrossLink(F_x, F_y, color, border))

    #for x in range(0,len(region_record_list)-1):
    x = 0
    for n, record in enumerate(region_record_list):
        gd_feature_set = feature_sets[n]
        i = 0

        if plasmid_list[x]:
            #print "PLASMID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            color1 = colors.HexColor('#2837B7')
            color2 = colors.blue
        else:
            color1 = colors.HexColor('#40F13A')
            color2 = colors.HexColor('#0F600C')



        for feature in record.features:
            if feature.type != "CDS":
                continue
            try:
                a = feature.qualifiers["locus_tag"]
            except:
                # cas des pseudogenes qui sont des CDS mais n'ont pas de protein ID
                continue

            if len(gd_feature_set) % 2 == 0:
                color = color1
            else:
                color = color2

            #try:
            #    try:
            #            group = protein_id2group[feature.qualifiers["protein_id"][0]]
            #    except:
            #            group = protein_id2group[feature.qualifiers["protein_id"][1]]
            #except:
            #    # no group attributed: singleton => special color
            #    color = colors.HexColor('#E104C0')


            for target_protein in target_protein_list:
                    if target_protein in feature.qualifiers["locus_tag"]:
                        print "target prot!"
                        color = colors.red

            gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_position="middle",label_strand=1, label_size=10, label_angle=40)
            i += 1
        x += 1

    print "max", max_len
    hauteur = 100*len(region_record_list)
    largeur = max(record_length)/30
    print "hauteur", hauteur
    print "largeur", largeur
    #gd_diagram.set_page_size(, orientation)
    if hauteur > largeur:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='portrait', fragments=1,start=0, end=max_len)
    else:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='landscape', fragments=1,start=0, end=max_len)
    print "writing diagram", out_name
    
    gd_diagram.write(out_name, "SVG")


    
def chunks(l, n):
    "return sublists of l of minimum length n (work subdivision for the subprocesing module"
    return [l[i:i+n] for i in range(0, len(l), n)]



#def reformat_record_list(biodb_name, bioentry_id_list, out_q):

#    server, db = manipulate_biosqldb.load_db(biodb_name)
#    records = [db.lookup(accession=one_bioentry) for one_bioentry in bioentry_id_list]

#    record_result = [_record_to_dict(temp_record) for temp_record in records]
    #record_dict = [i.to_dict() for i in record_result]
    #print "n records", len(result)
    #print record_dict[0]
    #return record_result

out_q = Queue()

def proteins_id2cossplot(server, biodb, biodb_name, locus_tag_list, out_name, region_size_bp, loaded_records):
    plasmid_list = []
    sub_record_list = []
    
    bioentry_id_list, seqfeature_id_list = manipulate_biosqldb.get_bioentry_and_seqfeature_id_from_locus_tag_list(server, locus_tag_list, biodb_name)
    #print "bioentry_id_list", bioentry_id_list
    #print "seqfeature_id_list", seqfeature_id_list
    #print "getting records..."
    
    #print "formatting records..."
   
    #n_cpu = 8
    #n_poc_per_list = math.ceil(len(bioentry_id_list)/float(n_cpu))
    #query_lists = chunks(range(0, len(bioentry_id_list)), int(n_poc_per_list))
    
    #procs = []
    #for one_list in query_lists:
    #    bioentry_list = [bioentry_id_list[i] for i in one_list]
        
        #print len(comb_list)
    #    proc = Process(target=reformat_record_list, args=(biodb_name, bioentry_list, out_q))
    #    procs.append(proc)
    #    proc.start()

    
    # format: list of tuples: (seq1, seq2, identity)
    #out_q.join()
    #print "getting results"
    #reformat_records = []
    #for i in range(len(procs)):
    #    reformat_records += out_q.get()
    #print reformat_records
    
#[SeqRecord(temp_record.seq, id=temp_record.id, name=temp_record.name, description=temp_record.description, dbxrefs =temp_record.dbxrefs, features=temp_record.features, annotations=temp_record.annotations) for temp_record in records]
    
    #import time
    #time.sleep(5)

    #for proc in procs:
    #    proc.join()

    unloaded_bioentry = []
    loaded_bioentry = []

    print "bientry id list", bioentry_id_list
    print "loaded bioentry:"
    for i in loaded_records.keys():
        print i
    for bioentry in bioentry_id_list:
        if bioentry in loaded_records.keys():
            loaded_bioentry.append(loaded_records[bioentry])
        else:
            unloaded_bioentry.append(bioentry) 
            

    all_records = [biodb.lookup(accession=one_bioentry) for one_bioentry in bioentry_id_list]
    unloaded_records = [biodb.lookup(accession=one_bioentry) for one_bioentry in unloaded_bioentry]
    reformat_records = [SeqRecord(Seq(temp_record.seq.data, temp_record.seq.alphabet), id=temp_record.id, name=temp_record.name, description=temp_record.description, dbxrefs =temp_record.dbxrefs, features=temp_record.features, annotations=temp_record.annotations) for temp_record in unloaded_records]
    reformat_records = reformat_records + loaded_bioentry

    print "all_record :", reformat_records
    
    #print "creating seqfeature_id2seqfeature..."
    #seqfeature_id2seqfeature = manipulate_biosqldb.seqfeature_id2seqfeature_object_dict(*all_records)
    #print "done"
    #print "region size", region_size_bp
   
    for seqfeature_id, record in zip(seqfeature_id_list, reformat_records):
        #print "seqfeature_id", seqfeature_id
        #print "record", record
        target_feature_start,  target_feature_end, strand = manipulate_biosqldb.seqfeature_id2feature_location(server, seqfeature_id)
        
        #target = seqfeature_id2seqfeature[seqfeature_id]
        try:
            plasmid = record.features[0].qualifiers["plasmid"]
            plasmid_list.append(True)
            plas = True
        except:
            plas = False
            plasmid_list.append(False)
        if plas:
            genome_list.append(record)
        else:
            sub_record = get_feature_neighborhood(target_feature_start,  target_feature_end, record, region_size_bp, "rec")
            sub_record_list.append(sub_record)
            
    plot_multiple_regions_crosslink(locus_tag_list, sub_record_list, plasmid_list, out_name)

    
    for record in reformat_records:
        print "record id", record.id
        print "record name", record.name
        if record.id not in loaded_records.keys():
            
            loaded_records[record.id] = record
            
    return loaded_records


#server, db = manipulate_biosqldb.load_db("chlamydiales")
#proteins_id2cossplot(server, db, "chlamydiales", ["CAQ48442.1"], "first_test.svg", 16000)






def proteins_id2sub_record_list(server, biodb, biodb_name, locus_tag_list, region_size_bp):
    plasmid_list = []
    sub_record_list = []
    
    bioentry_id_list, seqfeature_id_list = manipulate_biosqldb.get_bioentry_and_seqfeature_id_from_locus_tag_list(server, locus_tag_list, biodb_name)
    #print "bioentry_id_list", bioentry_id_list
    #print "seqfeature_id_list", seqfeature_id_list
    #print "getting records..."
    all_records = [biodb.lookup(accession=one_bioentry) for one_bioentry in bioentry_id_list]
    #print "formatting records..."
    
    reformat_records = [SeqRecord(Seq(temp_record.seq.data, temp_record.seq.alphabet), id=temp_record.id, name=temp_record.name, description=temp_record.description, dbxrefs =temp_record.dbxrefs, features=temp_record.features, annotations=temp_record.annotations) for temp_record in all_records]
    #print "creating seqfeature_id2seqfeature..."
    #seqfeature_id2seqfeature = manipulate_biosqldb.seqfeature_id2seqfeature_object_dict(*all_records)
    #print reformat_records, reformat_records
   
    for seqfeature_id, record in zip(seqfeature_id_list, reformat_records):
        #print "seqfeature_id", seqfeature_id
        #print "record", record

        
        #target = seqfeature_id2seqfeature[seqfeature_id]
        
        target_feature_start,  target_feature_end, strand = manipulate_biosqldb.seqfeature_id2feature_location(server, seqfeature_id)

        print "target_feature_start,  target_feature_end, strand:", target_feature_start,  target_feature_end, strand

        
        try:
            plasmid = record.features[0].qualifiers["plasmid"]
            plasmid_list.append(True)
            plas = True
        except:
            plas = False
            plasmid_list.append(False)
        if plas:
            genome_list.append(record)
        else:
            sub_record = get_feature_neighborhood(target_feature_start, target_feature_end, record, region_size_bp, "rec")
            sub_record_list.append(sub_record)
            
    return sub_record_list

