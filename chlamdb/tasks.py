
from __future__ import absolute_import, unicode_literals
import string
from django.contrib.auth.models import User
from django.utils.crypto import get_random_string
from celery import shared_task, current_task
import time
from chlamdb.biosqldb import manipulate_biosqldb
from django.conf import settings
from django.core.cache import cache
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from django.shortcuts import render
from django.http import HttpResponse
from django.template.loader import render_to_string
import os
from django.template import Template, Context

@shared_task
def create_random_user_accounts2(request):
    # get data
    print("get data")
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 3,
                                    'percent': 10,
                                    'description': "Getting data"})

    time.sleep(1)
    print("plotting")
    current_task.update_state(state='PROGRESS',
                              meta={'current': 2,
                                    'total': 3,
                                    'percent': 10,
                                    'description': "Plotting circos"})
    time.sleep(1)
    fertig = True
    print(locals())
    html = '<h3>Fertig!</h3>'#render_to_string('chlamdb/celery_test.html', context=locals())
    return html#HttpResponse(html)

@shared_task
def run_circos(reference_taxon, target_taxons):

    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos
    from chlamdb.biosqldb import shell_command
    import ete3
    from chlamdb.plots import gbk2circos

    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 4,
                                    'percent': 25,
                                    'description': "Get reference genome"})

    reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)

    record_list = []
    for accession in reference_accessions:

        #print "reference accession", accession
        biorecord = cache.get(biodb + "_" + accession)

        if not biorecord:
            #print biodb + "_" + accession, "NOT in memory"
            new_record = db.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                     id=new_record.id, name=new_record.name,
                                                     description=new_record.description,
                                                     dbxrefs =new_record.dbxrefs,
                                                     features=new_record.features,
                                                     annotations=new_record.annotations)
            record_id = biorecord.id.split(".")[0]
            cache.set(biodb + "_" + record_id, biorecord)
            record_list.append(biorecord)
        else:
            #print biodb + "_" + accession, "IN memory"
            record_list.append(biorecord)

    ref_name = ''
    for i in reference_accessions:
        ref_name += i
    circos_file = "circos/%s.svg" % ref_name

    current_task.update_state(state='PROGRESS',
                              meta={'current': 2,
                                    'total': 4,
                                    'percent': 50,
                                    'description': "Get query genome(s)"})

    querries = manipulate_biosqldb.get_genome_accessions(server, biodb)

    target_accessions = [manipulate_biosqldb.taxon_id2accessions(server,int(i),biodb)[0] for i in target_taxons]

    target_accessions += reference_accessions

    draft_data = []
    for biorecord in record_list:
        draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

    home_dir = os.path.dirname(__file__)
    #print "home_dir", home_dir
    temp_location = os.path.join(home_dir, "../assets/circos/")

    #sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    sql_order = 'select A.taxon_1 from (select taxon_1, median_identity from comparative_tables.shared_og_av_id_%s where taxon_2=%s ' \
                ' union select taxon_2,median_identity from comparative_tables.shared_og_av_id_%s ' \
                ' where taxon_1=%s order by median_identity DESC) A;' % (biodb, reference_taxon, biodb, reference_taxon)
    try:
        sql_order = 'select taxon_2 from comparative_tables.core_orthogroups_identity_msa_%s where taxon_1=%s order by identity desc;' % (biodb, reference_taxon)
        #print sql_order
        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
        #print 'msa core identity order!'
    except:
        print('msa orthogroups_average_identity order!')

        # median identity
        sql_order = 'select taxon from (select taxon_2 as taxon, median_identity ' \
                    ' from comparative_tables.shared_og_av_id_%s where taxon_1=%s union ' \
                    ' select taxon_1, median_identity as taxon from comparative_tables.shared_og_av_id_%s' \
                    '  where taxon_2=%s) A order by median_identity desc;' % (biodb,
                                                                              reference_taxon,
                                                                              biodb,
                                                                              reference_taxon)
        '''
        # n shared orthogroups
        sql_order = 'select taxon_2 from comparative_tables.shared_orthogroups_%s where taxon_1=%s order by n_shared_orthogroups DESC;' % (biodb,
                                                                                                                  reference_taxon)
        print
        '''
        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]

    '''
    print tree
    t1 = ete3.Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    print t1

    node_list = []
    for node in t1.iter_leaves():
            node_list.append(node.name)

    print 'original order', node_list

    reference_index = node_list.index(reference_taxon)
    ordered_taxons = node_list[reference_index:] + node_list[:reference_index][::-1]
    '''

    highlight_BBH = True
    if highlight_BBH:
          
        sql_phylum = 'select family from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                     ' inner join taxid2species_%s t3 on t2.taxon_id=t3.taxon_id ' \
                     ' inner join species_curated_taxonomy_%s t4 on t3.species_id=t4.species_id ' \
                     ' where t1.name="%s" and t2.taxon_id=%s limit 1; ' % (biodb,
                                                                           biodb,
                                                                           biodb, 
                                                                           reference_taxon)
        reference_phylum = server.adaptor.execute_and_fetchall(sql_phylum,)[0][0]
        print(reference_phylum)
        
        try:
            sql = 'select locus_tag from blastnr.blastnr_%s t1 ' \
              ' inner join biosqldb.bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
              ' inner join biosqldb.biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' inner join blastnr.blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
              ' inner join custom_tables.locus2seqfeature_id_%s t5 ' \
              ' on t1.seqfeature_id=t5.seqfeature_id ' \
              ' where t1.hit_number=1 and t3.name="%s" and t4.family!="%s" and t1.query_taxon_id=%s;' % (biodb,
                                                                                                         biodb,
                                                                                                         biodb,
                                                                                                         reference_phylum,
                                                                                                         reference_taxon)
            BBH_color = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            sql2 = 'select locus_tag from custom_tables.locus2seqfeature_id_%s t1 ' \
                   ' left join blastnr.blastnr_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                   ' where taxon_id=%s and t2.seqfeature_id is NULL;' % (biodb, biodb, reference_taxon)
                   
            no_BBH_hit_color = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]
            
            if len(BBH_color) < 20:
                sql = 'select locus_tag from blastnr.blastnr_%s t1 ' \
                  ' inner join biosqldb.bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
                  ' inner join biosqldb.biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                  ' inner join blastnr.blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
                  ' inner join custom_tables.locus2seqfeature_id_%s t5 ' \
                  ' on t1.seqfeature_id=t5.seqfeature_id ' \
                  ' where t1.hit_number=2 and t3.name="%s" and t4.family!="%s" and t1.query_taxon_id=%s;' % (biodb,
                                                                                                             biodb,
                                                                                                             biodb,
                                                                                                             reference_phylum,
                                                                                                             reference_taxon)
                BBH_color = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

        except:
            BBH_color = []
            no_BBH_hit_color = []

    else:
        BBH_color = []
        no_BBH_hit_color = []

    current_task.update_state(state='PROGRESS',
                              meta={'current': 3,
                                    'total': 4,
                                    'percent': 75,
                                    'description': "Plotting with circos..."})

    myplot = circos.CircosAccession2multiplot(server,
                                              db,
                                              biodb,
                                              record_list,
                                              target_accessions,
                                              locus_highlight=BBH_color,
                                              locus_highlight2=no_BBH_hit_color,
                                              out_directory=temp_location,
                                              draft_fasta=draft_data,
                                              href="/chlamdb/locusx/",
                                              ordered_taxons=ordered_taxons)

    original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % ref_name
    with open(original_map_file, "r") as f:
        map_string = ''.join([line for line in f.readlines()])

    circos_html = '<!DOCTYPE html>\n' \
                  ' <html>\n' \
                  ' <body>\n' \
                  ' %s\n' \
                  ' <img src="%s.svg" usemap="#%s">' \
                  ' </body>\n' \
                  ' </html>\n' % (map_string, ref_name, ref_name)

    circos_new_file = '/assets/circos/circos_clic.html'


    current_task.update_state(state='PROGRESS',
                              meta={'current': 4,
                                    'total': 4,
                                    'percent': 100,
                                    'description': "Rendering..."})

    with open(settings.BASE_DIR + circos_new_file, "w") as f:
        f.write(circos_html)

    circos_svg_file = "circos/%s.svg" % ref_name
    circos_png_file = "circos/%s.png" % ref_name
    original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
    map_file = "circos/%s.html" % ref_name
    svg_file = "circos/%s.svg" % ref_name
    map_name = ref_name

    template = Template('''
            {% load staticfiles %}
            {% load static %}
            <div class="row">
              <a download="circos.svg" class="btn" href="{% static circos_svg_file %}"><i class="fa fa-download"></i> Download SVG</a>
              <a download="circos.png" class="btn" href="{% static circos_png_file %}"><i class="fa fa-download"></i> Download PNG</a>
            </div>

              <div class="row">
                <div class="col-lg-12">
                  <div>
                    <a href="{{ circos_new_file }}"><img src="{% static ""%}{{circos_file}}" id="circos_plot" style="position: absolute; top: 00px; left: 00px;"></a>
                    <img src="{% static "/scales/circos_legend.png" %}" id="circos_legend"  width="160" style="position: relative; top: 0; left: 0;">

                  </div>
                </div>
              </div>
            ''')

    html = template.render(Context(locals()))#render_to_string(template, context=locals())
    return html


@shared_task
def run_circos_main(reference_taxon, target_taxons, highlight):


    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos
    from chlamdb.biosqldb import shell_command
    import ete3
    from chlamdb.plots import gbk2circos

    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)


    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 4,
                                    'percent': 25,
                                    'description': "Get reference genome"})

    reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb) # ["NC_009648"] NZ_CP009208 NC_016845

    #print "reference_accessions", reference_accessions
    record_list = []
    for accession in reference_accessions:
        #print "reference accession", accession
        biorecord = cache.get(biodb + "_" + accession)

        if not biorecord:
            #print biodb + "_" + accession, "NOT in memory"
            new_record = db.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                     id=new_record.id, name=new_record.name,
                                                     description=new_record.description,
                                                     dbxrefs =new_record.dbxrefs,
                                                     features=new_record.features,
                                                     annotations=new_record.annotations)
            record_id = biorecord.id.split(".")[0]
            cache.set(biodb + "_" + record_id, biorecord)
            record_list.append(biorecord)
        else:
            #print biodb + "_" + accession, "In memory"
            record_list.append(biorecord)

    ref_name = ('').join(reference_accessions)

    circos_file = "circos/%s.svg" % ref_name


    current_task.update_state(state='PROGRESS',
                              meta={'current': 2,
                                    'total': 4,
                                    'percent': 50,
                                    'description': "Get query genome(s)"})
    

    querries = manipulate_biosqldb.get_genome_accessions(server, biodb)

    target_accessions = [manipulate_biosqldb.taxon_id2accessions(server,int(i),biodb)[0] for i in target_taxons]

    target_accessions += reference_accessions



    draft_data = []
    for biorecord in record_list:
        draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

    home_dir = os.path.dirname(__file__)

    temp_location = os.path.join(home_dir, "../assets/circos/")

    #sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    sql_order1 = 'select A.taxon_1 from (select taxon_1,median_identity from comparative_tables.shared_og_av_id_%s where taxon_2=%s ' \
                ' union select taxon_2,median_identity from comparative_tables.shared_og_av_id_%s ' \
                ' where taxon_1=%s order by median_identity DESC) A;' % (biodb, reference_taxon, biodb, reference_taxon)
    try:
        sql_order = 'select taxon_2 from comparative_tables.core_orthogroups_identity_msa_%s where taxon_1=%s order by identity desc;' % (biodb, reference_taxon)

        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
    except:
        sql_order2 = 'select taxon_2 from comparative_tables.shared_og_av_id_%s where taxon_1=%s order by median_identity desc;' % (biodb, reference_taxon)

        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order1)]

    current_task.update_state(state='PROGRESS',
                              meta={'current': 3,
                                    'total': 4,
                                    'percent': 75,
                                    'description': "Plotting with circos..."})

    myplot = circos.CircosAccession2multiplot(server,
                              db,
                              biodb,
                              record_list,
                              target_accessions,
                              locus_highlight=highlight,
                              out_directory=temp_location,
                              draft_fasta=draft_data,
                              href="/chlamdb/locusx/",
                              ordered_taxons = ordered_taxons)



    original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % ref_name
    with open(original_map_file, "r") as f:
        map_string = ''.join([line for line in f.readlines()])

    circos_html = '<!DOCTYPE html>\n' \
                  ' <html>\n' \
                  ' <body>\n' \
                  ' %s\n' \
                  ' <img src="%s.svg" usemap="#%s">' \
                  ' </body>\n' \
                  ' </html>\n' % (map_string, ref_name, ref_name)


    circos_new_file = '/assets/circos/circos_clic.html'
    #print settings.BASE_DIR + circos_new_file
    with open(settings.BASE_DIR + circos_new_file, "w") as f:
        f.write(circos_html)

    #target_map_file = settings.BASE_DIR + "/templates/circos/%s.html" % ref_name
    original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
    #target_map_file_svg = settings.BASE_DIR + "/templates/circos/%s.svg" % ref_name
    map_file = "circos/%s.html" % ref_name
    svg_file = "circos/%s.svg" % ref_name
    #a, b, c = shell_command.shell_command("mv %s %s" % (original_map_file, target_map_file))
    #a, b, c = shell_command.shell_command("cp %s %s" % (original_map_file_svg, target_map_file_svg))
    #print a,b,c
    map_name = ref_name

    template = Template('''
            {% load staticfiles %}
            {% load static %}
            <div class="row">
              <a download="circos.svg" class="btn" href="{% static circos_svg_file %}"><i class="fa fa-download"></i> Download SVG</a>
              <a download="circos.png" class="btn" href="{% static circos_png_file %}"><i class="fa fa-download"></i> Download PNG</a>
            </div>

              <div class="row">
                <div class="col-lg-12">
                  <div>
                    <a href="{{ circos_new_file }}"><img src="{% static ""%}{{circos_file}}" id="circos_plot" style="position: absolute; top: 00px; left: 00px;"></a>
                    <img src="{% static "/scales/circos_legend.png" %}" id="circos_legend"  width="160" style="position: relative; top: 0; left: 0;">

                  </div>
                </div>
              </div>
            ''')

    html = template.render(Context(locals()))#render_to_string(template, context=locals())
    return html