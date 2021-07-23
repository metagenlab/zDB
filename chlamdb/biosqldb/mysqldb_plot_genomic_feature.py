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
from chlamdb.biosqldb import orthogroup_identity_db

def get_feature_neighborhood(feature_start, feature_end, contig_or_genome_record, neighborhood_size_bp, record_name):

    start = feature_start - neighborhood_size_bp
    end = feature_end + neighborhood_size_bp

    if start < 0:
        diff = abs(0-start)
        end += diff
        start = 0

    #if contig_or_genome_record.features[0].location.start > start:
    #        start = 0
    #if contig_or_genome_record.features[0].location.end < end:
    #        end = contig_or_genome_record.features[0].location.end
    if end > len(contig_or_genome_record.seq):
        diff = abs(len(contig_or_genome_record.seq) - end)
        start = start-diff
        if start < 0:
            start = 0
        record = contig_or_genome_record[int(start):]
    else:
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
        #print "i", i
        #Allocate tracks 3 (top), 1 (bottom) for region 1 and 2
        #(empty tracks 2 useful white space to emphasise the cross links
        #and also serve to make the tracks vertically more compressed)
        gd_track_for_features = gd_diagram.new_track((2*n_records-1)-2*i, name=record.name, greytrack=True, height=0.5, start=0, end=len(record))
        if record.name not in feature_sets:
                feature_sets.append(gd_track_for_features.new_set())
        else:
           print ("already in feature_sets!")
           print (record)
           quit

    for x in range(0,len(region_record_list)-1):
        #print "x", x
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
            #print "PLASMID!!!"
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
                        #print "target prot!"
                        color = colors.red

            gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_position="middle",label_strand=1, label_size=12, label_angle=45)
            i += 1
        x += 1

    #print "max", max_len
    #print "n records", len(region_record_list)
    if len(region_record_list) == 2:
        hauteur = 700
    else:
        hauteur = 250*len(region_record_list)
    largeur = max(record_length)/30
    #print "hauteur", hauteur
    #print "largeur", largeur
    #gd_diagram.set_page_size(, orientation)
    if hauteur > largeur:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='portrait', fragments=1,start=0, end=max_len)
    else:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='landscape', fragments=1,start=0, end=max_len)
    #print "writing diagram", out_name

    gd_diagram.write(out_name, "SVG")





def clamp(val, minimum=0, maximum=255):
    if val < minimum:
        return minimum
    if val > maximum:
        return maximum
    return val



def colorscale(hexstr, scalefactor):
    """
    Scales a hex string by ``scalefactor``. Returns scaled hex string.

    To darken the color, use a float value between 0 and 1.
    To brighten the color, use a float value greater than 1.

    >>> colorscale("#DF3C3C", .5)
    #6F1E1E
    >>> colorscale("#52D24F", 1.6)
    #83FF7E
    >>> colorscale("#4F75D2", 1)
    #4F75D2
    """

    hexstr = hexstr.strip('#')

    if scalefactor < 0 or len(hexstr) != 6:
        return hexstr

    r, g, b = int(hexstr[:2], 16), int(hexstr[2:4], 16), int(hexstr[4:], 16)

    r = clamp(r * scalefactor)
    g = clamp(g * scalefactor)
    b = clamp(b * scalefactor)

    return "#%02x%02x%02x" % (r, g, b)



def plot_multiple_regions_crosslink(target_protein_list,
                                    region_record_list,
                                    plasmid_list,
                                    out_name,
                                    biodb_name="chlamydia_03_15",
                                    color_locus_list = [],
                                    flip_record_based_on_first=True):


    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl
    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']

    norm = mpl.colors.Normalize(vmin=-30, vmax=100)
    cmap = cm.Blues
    m = cm.ScalarMappable(norm=norm, cmap=cmap)

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb_name)

    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    gd_diagram = GenomeDiagram.Diagram("geomic_region")
    feature_sets = []
    max_len = 0
    records = dict((rec.name, rec) for rec in region_record_list)

    n_records = len(region_record_list)

    record_length = [len(record) for record in region_record_list]

    if flip_record_based_on_first:
        region_record_list_flip = [region_record_list[0]]
        region_record_list_flip[0].name = region_record_list_flip[0].description
        for x in range(0,len(region_record_list)-1):
            same_strand_count = 0
            different_strand_count = 0
            features_X = region_record_list[x].features
            features_Y = region_record_list[x+1].features
            for feature_1 in features_X:

                if feature_1.type != "CDS":
                    continue

                for feature_2 in features_Y:
                    if feature_2.type != "CDS":
                        continue
                    try:

                        group1 = feature_1.qualifiers["orthogroup"][0]
                        group2 = feature_2.qualifiers["orthogroup"][0]
                        if group1 == group2:
                            strand1 = feature_1.location.strand
                            strand2 = feature_2.location.strand
                            if strand1 == strand2:
                                same_strand_count+=1
                            else:
                                different_strand_count+=1

                    except:
                            pass

            if different_strand_count>same_strand_count:
                region_record_list[x+1] = region_record_list[x+1].reverse_complement(id=region_record_list[x+1].id,
                                                                                     name=region_record_list[x+1].description)
            else:
                region_record_list[x+1].name = region_record_list[x+1].description


        #region_record_list = region_record_list_flip
    for i, record in enumerate(region_record_list):
        max_len = max(max_len, len(record))
        #Allocate tracks 3 (top), 1 (bottom) for region 1 and 2
        #(empty tracks 2 useful white space to emphasise the cross links
        #and also serve to make the tracks vertically more compressed)
        gd_track_for_features = gd_diagram.new_track((1*n_records-1)-1*i, name=record.name, greytrack=True, height=0.4, start=0, end=len(record))
        if record.name not in feature_sets:
                feature_sets.append(gd_track_for_features.new_set())
        else:
           print ("already in feature_sets!")
           print (record)
           quit



    #print 'looping....'
    for x in range(0,len(region_record_list)-1):
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

                if group1 == group2:
                    border = colors.lightgrey
                    color = colors.lightgrey
                    try:
                        identity = orthogroup_identity_db.check_identity(cursor, 
                                                                        feature_1.qualifiers["orthogroup"][0],
                                                                        feature_1.qualifiers["locus_tag"][0],
                                                                        feature_2.qualifiers["locus_tag"][0])
                    except:
                        identity = 0
                        print ("problem with identity table %s and locus %s %s" % (group1,
                                                                                  feature_1.qualifiers["locus_tag"][0],
                                                                                  feature_2.qualifiers["locus_tag"][0]))


                    color2 = colors.HexColor(rgb2hex(m.to_rgba(float(identity))))
                    border2 = colors.HexColor(rgb2hex(m.to_rgba(float(identity))))

                    F_x = set_X.add_feature(SeqFeature(FeatureLocation(feature_1.location.start, feature_1.location.end, strand=0)),
                                    color=color, border=border, set_id=feature_1.qualifiers["locus_tag"])
                    F_y = set_Y.add_feature(SeqFeature(FeatureLocation(feature_2.location.start, feature_2.location.end, strand=0)),
                                    color=color, border=border)
                    gd_diagram.cross_track_links.append(CrossLink(F_x, F_y, color2, border2))


    #for x in range(0,len(region_record_list)-1):
    x = 0
    all_locus = []

    for n, record in enumerate(region_record_list):
        gd_feature_set = feature_sets[n]
        i = 0


        if plasmid_list[x]:
            #print "PLASMID!!"
            color1 = colors.HexColor('#2837B7')
            color2 = colors.blue
        else:
            color1 = colors.HexColor('#40F13A')
            color2 = colors.HexColor('#0F600C')


        one_row_locus = []
        for feature in record.features:
            if feature.type == "tblast_target":
                feature.name = 'match'
                gd_feature_set.add_feature(feature, sigil="BOX", color="#ff4a0c86", label=False, label_position="middle", label_size=25, label_angle=0)


            if feature.type == "assembly_gap":
                #print "gap", feature
                feature.location.strand = None
                gd_feature_set.add_feature(feature, sigil="BOX", color="red", label=True, label_position="middle", label_strand=1, label_size=14, label_angle=40)

            if feature.type == "rRNA":

                gd_feature_set.add_feature(feature, sigil="ARROW", color="orange", label=True, label_position="middle", label_strand=1, label_size=10, label_angle=40)
                try:
                    one_row_locus.append(feature.qualifiers["locus_tag"][0])
                except:
                    pass
            if feature.type == "tRNA":

                gd_feature_set.add_feature(feature, sigil="ARROW", color="orange", label=True, label_position="middle", label_strand=1, label_size=10, label_angle=40)
                try:
                    one_row_locus.append(feature.qualifiers["locus_tag"][0])
                except:
                    print ('no locus tag for:')
                    print (feature)


            if feature.type == "repeat_region":

                gd_feature_set.add_feature(feature, sigil="BOX", color="blue", label=True, label_position="middle", label_strand=1, label_size=14, label_angle=40)

            if 'pseudo' in feature.qualifiers:

                gd_feature_set.add_feature(feature, sigil="OCTO", color="#6E6E6E", label=True, label_position="middle", label_strand=1, label_size=10, label_angle=40)


            elif feature.type != "CDS":
                continue
            else:

                try:
                    a = feature.qualifiers["locus_tag"][0]
                except:
                    # cas des pseudogenes qui sont des CDS mais n'ont pas de protein ID
                    continue


                if a in color_locus_list:
                    #print '###########################', a, color_locus_list
                    if len(gd_feature_set) % 2 == 0:
                        color = colors.HexColor('#ca4700')
                    else:
                        color = colors.HexColor('#fd7a32')
                else:
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
                            #print "target prot!"
                            color = colors.red

                gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_position="middle",label_strand=1, label_size=10, label_angle=40)
                i += 1
                try:
                    one_row_locus.append(feature.qualifiers["locus_tag"][0])
                except:
                    print ('no locus tag for:')
                    print (feature)
        all_locus = one_row_locus + all_locus


        x += 1

    #print "max", max_len
    #print "n record", len(region_record_list)

    if len(region_record_list) == 2:
        hauteur = 300
    else:
        hauteur = 150*len(region_record_list)
    largeur = max(record_length)/30
    #print "hauteur", hauteur
    #print "largeur", largeur
    #gd_diagram.set_page_size(, orientation)
    if hauteur > largeur:
            gd_diagram.draw(format="linear", pagesize=(hauteur, largeur), orientation='portrait', fragments=1,start=0, end=max_len)
    else:
            gd_diagram.draw(format="linear", pagesize=(hauteur, largeur), orientation='landscape', fragments=1,start=0, end=max_len)
    #print "writing diagram", out_name

    #gd_diagram.write(out_name, "SVG")



    import io
    from chlamdb.plots import edit_svg

    svg_diagram = io.StringIO()
    gd_diagram.write(svg_diagram, "SVG")
    svg_diagram.flush()
    #gd_diagram

    with_links = edit_svg.edit_svg(svg_diagram.getvalue(), all_locus, biodb_name)

    with_links.write(out_name)

    png_name = out_name.split('.')[0] + '.png'

    #png_handle = open(png_name, 'w')
    #gd_diagram.write(png_handle, "PNG")
    #png_handle.close()

    try:
        cmd = 'chmod 444 %s' % out_name
    except:
        pass
    from chlamdb.biosqldb import shell_command
    #print cmd
    shell_command.shell_command(cmd)

    return all_locus



def plot_simple_region(region_record, out_name):
    gd_diagram = GenomeDiagram.Diagram("geomic_region")

    gd_track_for_features = gd_diagram.new_track(1, name=region_record.name, greytrack=True, height=0.5, start=0, end=len(region_record))

    gd_feature_set = gd_track_for_features.new_set()

    color1 = colors.HexColor('#40F13A')
    color2 = colors.HexColor('#0F600C')

    for feature in region_record.features:
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
        gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_position="middle",label_strand=1, label_size=10, label_angle=40)

    hauteur = 250
    largeur = len(region_record)/30
    #print "hauteur", hauteur
    #print "largeur", largeur

    if hauteur > largeur:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='portrait', fragments=1,start=0, end=len(region_record))
    else:
            gd_diagram.draw(format="linear", pagesize=(hauteur,largeur), orientation='landscape', fragments=1,start=0, end=len(region_record))
    #print "writing diagram", out_name

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

def proteins_id2cossplot(server, 
                         biodb, 
                         biodb_name, 
                         locus_tag_list, 
                         out_name, 
                         region_size_bp, 
                         cache, 
                         color_locus_list = []):
    
    plasmid_list = []
    sub_record_list = []

    bioentry_id_list, seqfeature_id_list = manipulate_biosqldb.get_bioentry_and_seqfeature_id_from_locus_tag_list(server, locus_tag_list, biodb_name)

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

    #print "bientry id list", bioentry_id_list
    #print "loaded bioentry:"
    #for i in loaded_records.keys():
    #    print i

    for bioentry in bioentry_id_list:
        key = biodb_name + "_" + bioentry
        biorecord = cache.get(key)
        if biorecord:
            #print key, "in memory"
            #loaded_bioentry.append(loaded_records[bioentry])
            continue
        else:
            #print key, "NOT in memory"
            cache_time = None
            #unloaded_bioentry.append(bioentry)
            new_record = biodb.lookup(accession=bioentry)
            new_record_reformat = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                            id=new_record.id, name=new_record.name,
                                            description=new_record.description,
                                            dbxrefs =new_record.dbxrefs,
                                            features=new_record.features,
                                            annotations=new_record.annotations)
            record_id = biodb_name + "_" + new_record_reformat.id.split(".")[0]
            cache.set(key, new_record_reformat, cache_time)
            biorecord = cache.get(key)
            if biorecord:
                print (key, "in memory")

    #all_records = [biodb.lookup(accession=one_bioentry) for one_bioentry in bioentry_id_list]
    #unloaded_records = [biodb.lookup(accession=one_bioentry) for one_bioentry in unloaded_bioentry]
    #reformat_records = [SeqRecord(Seq(temp_record.seq.data, temp_record.seq.alphabet), id=temp_record.id, name=temp_record.name, description=temp_record.description, dbxrefs =temp_record.dbxrefs, features=temp_record.features, annotations=temp_record.annotations) for temp_record in unloaded_records]
    #reformat_records = reformat_records + loaded_bioentry

    #print "all_record :", cache

    #print "creating seqfeature_id2seqfeature..."
    #seqfeature_id2seqfeature = manipulate_biosqldb.seqfeature_id2seqfeature_object_dict(*all_records)
    #print "done"
    #print "region size", region_size_bp

    orthogroup_list = []
    for seqfeature_id, record_id in zip(seqfeature_id_list, bioentry_id_list):
        key = biodb_name + "_" + record_id
        record = cache.get(key)
        if not record:
            print (key, "still not in memory")

        #print "seqfeature_id", seqfeature_id
        #print "record OKKKKK"
        target_feature_start, target_feature_end, strand = manipulate_biosqldb.seqfeature_id2feature_location(server, seqfeature_id)
        #print "target_feature_start,  target_feature_end strand", target_feature_start,  target_feature_end, strand
        #target = seqfeature_id2seqfeature[seqfeature_id]
        try:
            plasmid = record.features[0].qualifiers["plasmid"]
            plasmid_list.append(True)
            plas = True
        except:
            plas = False
            plasmid_list.append(False)
        if plas:
            sub_record = get_feature_neighborhood(target_feature_start,  target_feature_end, record, region_size_bp, "rec")
            sub_record_list.append(sub_record) # record
            for feature in sub_record.features:
             if feature.type == 'CDS':
                 try:
                     orthogroup_list.append(feature.qualifiers['orthogroup'][0])
                 except:
                     pass
        else:
            #print "Getting region from %s" % record_id
            #print record
            sub_record = get_feature_neighborhood(target_feature_start,  target_feature_end, record, region_size_bp, "rec")
            sub_record_list.append(sub_record)
            #print
            for feature in sub_record.features:
             if feature.type == 'CDS':
                 try:
                     orthogroup_list.append(feature.qualifiers['orthogroup'][0])
                 except:
                     pass

    region_locus_list = plot_multiple_regions_crosslink(locus_tag_list,
                                                        sub_record_list,
                                                        plasmid_list,
                                                        out_name,
                                                        biodb_name,
                                                        color_locus_list=color_locus_list)
    return region_locus_list, orthogroup_list

def location2plot(db,
                  accession,
                  out_name,
                  start,
                  end,
                  cache,
                  color_locus_list = [],
                  region_highlight=[]):
    import copy
    if start < 0:
        start=0
    key = accession
    biorecord = cache.get(key)
    if biorecord:
        print (key, "in memory")
    else:
        print (key, "NOT in memory")
        cache_time = None

        db=db.server[db.db_name]
        new_record = db.lookup(accession=accession)
        new_record_reformat = SeqRecord(Seq(new_record.seq.data),
                                                         id=new_record.id, name=new_record.name,
                                                         description=new_record.description,
                                                         dbxrefs =new_record.dbxrefs,
                                                         features=new_record.features,
                                                         annotations=new_record.annotations)
        cache.set(key, new_record_reformat, cache_time)
        biorecord = cache.get(key)
        if biorecord:
            print (key, "in memory")

    fake_feature = copy.copy(biorecord.features[1])
    fake_feature.type = "tblast_target"
    print(region_highlight[0], region_highlight[1])
    fake_feature.location = FeatureLocation(region_highlight[0], region_highlight[1], strand=0)
    biorecord.features.append(fake_feature)
    print("start-end",start,end)
    sub_record = biorecord[start:end]
    print(sub_record.features)
    if len(sub_record.features) > 0:
        sub_record.features = ([sub_record.features[-1]] + sub_record.features[0:-1])
    region_locus_list = plot_multiple_regions_crosslink([],
                                                        [sub_record],
                                                        [False],
                                                        out_name,
                                                        color_locus_list=color_locus_list)
    return region_locus_list




def location2simpleplot(biodb, biodb_name, bioentry, location_start, location_stop, out_name, cache):

    key = biodb_name + "_" + bioentry
    biorecord = cache.get(key)
    if biorecord:
        print (key, "in memory")

    else:
        print (key, "NOT in memory")

        new_record = biodb.lookup(accession=bioentry)
        new_record_reformat = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
        record_id = biodb_name + "_" + new_record_reformat.id.split(".")[0]
        cache.set(record_id, new_record_reformat)


    record = cache.get(biodb_name + "_" + bioentry)
    sub_record = record[location_start:location_stop]
    plot_simple_region(sub_record, out_name)


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

        #print "target_feature_start,  target_feature_end, strand:", target_feature_start,  target_feature_end, strand


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
