
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
from django.core.cache import cache


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
def extract_orthogroup_task(biodb, 
                            include,
                            exclude,
                            freq_missing,
                            single_copy,
                            accessions,
                            reference_taxon,
                            fasta_url,
                            n_missing):
  
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.views import get_locus_annotations
    
    server = manipulate_biosqldb.load_db()
    
    
    current_task.update_state(state='PROGRESS',
    meta={'current': 1,
          'total': 3,
          'percent': 25,
          'description': "Get orthogroup list"})
    
    if int(n_missing)>=len(include):
        template = Template('''
        {% load staticfiles %}
        {% load static %}
        <div class="panel panel-warning" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
            <div class="panel-heading" style="width:100%">
                <h3 class="panel-title">Help</h3>
            </div>
            <p style="margin: 10px 10px 10px 10px">You cannot set a number of missing data bigger than the number of included genomes</p>
        </div>
        ''')

        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html   
    
    print("get matrix")
    if not accessions:
        # get sub matrix and complete matrix
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                  "orthology",
                                                                  "orthogroup",
                                                                  include,
                                                                  exclude,
                                                                  freq_missing,
                                                                  single_copy=single_copy,
                                                                  accessions=accessions,
                                                                  cache=cache)
    else:
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                  "orthology",
                                                                  "id",
                                                                  include,
                                                                  exclude,
                                                                  freq_missing,
                                                                  single_copy=single_copy,
                                                                  accessions=accessions,
                                                                  cache=cache)
    print("matrix OK")
    match_groups = mat.index.tolist()

    if len(match_groups) == 0:
        template = Template('''
        {% load staticfiles %}
        {% load static %}
        <div class="panel panel-warning" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
        <div class="panel-heading" style="width:100%">
            <h3 class="panel-title"></h3>
        </div>
        <p style="margin: 10px 10px 10px 10px">No match</p>
        </div>
        ''')

        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html
    else:

        # get count in subgroup
        orthogroup2count = dict((mat > 0).sum(axis=1))
        # get count in complete database
        orthogroup2count_all = dict((mat_all > 0).sum(axis=1))

        #print cog2count_all
        max_n = max(orthogroup2count_all.values())

        # GET max frequency for template
        sum_group = len(match_groups)

        current_task.update_state(state='PROGRESS',
        meta={'current': 2,
              'total': 3,
              'percent': 50,
              'description': "Retrieve annotations"})

        match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(match_groups,
                                                                                                      biodb,
                                                                                                      taxon_filter=include,
                                                                                                      accessions=accessions)
        print("done")

        if len(include) == 1:
            # get url to get single include taxon fasta
            locus_list = [i[2] for i in extract_result]
            fasta_ref_url = '?l=' + '&l='.join(locus_list)

        columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                  'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'

        envoi_extract = True

        circos_url = '?ref=%s&' % reference_taxon
        circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(match_groups)
        fasta_url+= "&i="+('&i=').join(include)
        fasta_url+= "&e="+('&e=').join(exclude)
        fasta_url+= "&f=%s" % freq_missing
        fasta_url+= "&s=%s" % single_copy
        fasta_url_ref = fasta_url +'&ref=%s' % reference_taxon
        fasta_url_noref =fasta_url + '&ref=F'

        current_task.update_state(state='PROGRESS',
        meta={'current': 3,
              'total': 3,
              'percent': 90,
              'description': "Retrieve reference genome annotation"})

        if not accessions:
            sql = 'select locus_tag from orthology_detail_%s where orthogroup in (%s) and taxon_id=%s' % (biodb,
                                                                                                        '"' + '","'.join(match_groups) + '"',
                                                                                                        reference_taxon)
        else:
            sql = 'select locus_tag from orthology_detail_%s where orthogroup in (%s) and accession="%s"' % (biodb,
                                                                                                        '"' + '","'.join(match_groups) + '"',
                                                                                                        reference_taxon)
        locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


        fasta_ref_url = '?l=' + '&l='.join(locus_list)

        locus2annot, \
        locus_tag2cog_catego, \
        locus_tag2cog_name, \
        locus_tag2ko, \
        pathway2category, \
        module2category, \
        ko2ko_pathways, \
        ko2ko_modules,\
        locus2interpro = get_locus_annotations(biodb, locus_list)    

        template = Template('''
        {% load staticfiles %}
        {% load static %}
        {% load custom_tags %}

        <div class="row" style="padding-top:20px">
        <div class="col-lg-12">
          <ul id="tabs" class="nav nav-tabs" data-tabs="tabs">
              <li class="active"><a href="#red" data-toggle="tab">Orthogroups</a></li>
              <li><a href="#orange" data-toggle="tab">Table detail</a></li>
              <li><a href="#tab2" data-toggle="tab">Locus list reference</a></li>
          </ul>

          <div id="my-tab-content" class="tab-content">
              <div class="tab-pane active" id="red">
                <div style="padding-top:10px;">
                    <div id="export_bouttons_groups">
                        <a href="{% url 'get_fasta' %}{{fasta_url_noref}}" class="btn btn-success">Download fasta</a>
                        <a href="{% url 'circos_main' %}{{circos_url}}" class="btn btn-success">Show on circular map (takes some time)</a>
                        <a href="{% url 'orthogroup_list_cog_barchart' accessions%}{{circos_url}}" class="btn btn-success">Show cog categories of orthogroup subset vs all orthogroups</a>
                        <br/>
                    </div>

                    <h3>Number of orthogroups: {{ sum_group }} </h3>


                      <div class="panel panel-success" style="width:100%; top: 200px; margin: 10px 10px 10px 10px">
                          <div class="panel-heading" style="width:100%">
                              <h3 class="panel-title">Help</h3>
                          </div>
                          <p style="margin: 10px 10px 10px 10px; line-height: 180%">The annotation(s) of orthologous groups is a consensus of the annotation of all members of the group. For each annotation (gene name, product, COG,...), 
                          only the two most frequent annotations are reported. Number in parentheses indicates the number of occurences of the annotation (how many proteins of the orthologous group have this exact annotation). </br>
                          The second tab reports the complete list of protein of reference genome(s).                               
                          </p>
                      </div>


                    <table class="display" id="table_groups">

                        <thead>
                            <tr>
                            <th></th>
                            <th id="entete_locus">Orthogroup</th>
                            <th>Genes</th>
                            <th>Products</th>
                            <th>Cogs</th>
                            <th>Pfam</th>
                            <th>present in (/{{include|length}})</th>
                            <th>freq complete database {{max_n}}</th>
                            </tr>
                        </thead>

                        {% for value in match_groups_data %}
                            <tr>
                                <td>{{value.0}}</td>
                                <td><a href="{% url 'locusx' value.1 True%}">{{value.1}}</a></td>
                                <td>{{value.2|safe}}</td>
                                <td>{{value.3|safe}}</td>
                                <td>{{value.4|safe}}</td>
                                <td>{{value.5|safe}}</td>
                                <td>{{orthogroup2count|keyvalue:value.1}}</td>
                                <td>{{orthogroup2count_all|keyvalue:value.1}}</td>
                            </tr>
                        {% endfor %}


                    </table>
                </div>
              </div>
              <div class="tab-pane" id="orange">
                      <div id="export_bouttons_groups">
                          {% if show_reference_annot %}
                          <button type="button" class="btn btn-primary btn-xs" onclick="location.href='/chlamdb/fasta/{{fasta_ref_url|safe}}'">fasta aa</button>
                          {% endif %}
                            <br/>
                      </div>


                      <h3>Number of orthogroups: {{ sum_group }} </h3>

                      <table class="display" id="cog_table">
                          <thead>
                              <tr>
                                  <th></th>
                                  <th>Orthogroup</th>
                                  <th>Locus</th>
                                  <th>Protein id</th>
                                  <th>Start</th>
                                  <th>Stop</th>
                                  <th>S.</th>
                                  <th>Gene</th>
                                  <th>nH</th>
                                  <th>nG</th>
                                  <th>TM</th>
                                  <th>SP</th>
                                  <th>Product</th>
                                  <th>Organism</th>
                                  <th>Translation</th>
                              </tr>
                          </thead>
                          <tbody>

                              {% for values in extract_result %}
                                  <tr>
                                      <td>{{values.0}}</td>
                                      <td><a href="{% url 'locusx' values.1 True %}"  target="_top">{{values.1}}</a></td>
                                      <td>{{values.2}}</td>
                                      <td>{{values.3}}</td>
                                      <td>{{values.4}}</td>
                                      <td>{{values.5}}</td>
                                      <td>{{values.6}}</td>
                                      <td>{{values.7}}</td>
                                      <td>{{values.8}}</td>
                                      <td>{{values.9}}</td>
                                      <td>{{values.10}}</td>
                                      <td>{{values.11}}</td>
                                      <td>{{values.12}}</td>
                                      <td>{{values.13}}</td>
                                      <td>{{values.14}}</td>
                                  </tr>
                              {% endfor %}
                          </tbody>
                      </table>
              </div>
              <div class="tab-pane" id="tab2">

                  <h3> Locus annotation </h3>
                  <div id="export_bouttons_groups" style="width:100%; top: 200px; margin: 10px 10px 10px 10px">
                      <a href="{% url "get_fasta" %}{{fasta_url_ref}}" class="btn btn-success">Download fasta</a>
                      <a href="{% url 'circos_main' %}{{circos_url}}" class="btn btn-success">Show on circular map</a>
                  <br/>
                  </div>
                  <table id="reference_table" class="table">
                      <thead>
                          <tr>
                              <th></th>
                              <th id="entete_locus">Orthogroup</th>
                              <th id="entete_locus">Locus</th>
                              <th style="width:10px">C</th>
                              <th style="width:70px">COGn</th>
                              <th style="width:55px">KO</th>
                              <th style="width:260px">Pathways</th>
                              <th style="width:230px">Modules</th>
                              <th style="width:230px">Interpro</th>
                              <th id="entete_gene">Gene</th>
                              <th id="entete_n">nH</th>
                              <th id="entete_n">nG</th>
                              <th id="entete_n">TM</th>
                              <th id="entete_n">SP</th>
                              <th id="entete_product">Product</th>


                          </tr>
                      </thead>
                      <tbody>

                          {%for values in locus2annot%}
                          <tr>
                              <td>{{values.0}}</td>
                              <td>{{values.1}}</td>
                              <td><a href="{% url 'locusx' values.2 True %}" target="_top">{{values.2}}</a></td>
                              <td>{{locus_tag2cog_catego|keyvalue:values.2}}</td>
                              {% with locus_tag2cog_name|keyvalue:values.2 as name %}
                                  {% if name == '-' %}
                                      <td>{{locus_tag2cog_name|keyvalue:values.2}}</td>
                                  {% else %}
                                      <td><a href="{% url "fam" name 'cog'%}" target="_top">{{locus_tag2cog_name|keyvalue:values.2}}</a></td>
                                  {% endif %}
                              {%endwith%}
                              {% with locus_tag2ko|keyvalue:values.2 as oneko %}
                                  {% if oneko == '-' %}
                                      <td>{{locus_tag2ko|keyvalue:values.2}}</td>
                                      <td>-</td>
                                      <td>-</td>
                                  {% else %}
                                      <td><a href="{% url "fam" oneko 'ko'%}" target="_top">{{locus_tag2ko|keyvalue:values.2}}</a></td>
                                      {% with ko2ko_pathways|keyvalue:oneko as name %}
                                          {% if name == '-' %}
                                              <td>-</td>
                                          {% else %}
                                              <td>{{name|safe}}</td>
                                          {% endif %}
                                      {%endwith%}
                                      {% with ko2ko_modules|keyvalue:oneko as name %}
                                          {% if name == '-' %}
                                              <td>-</td>
                                          {% else %}
                                              <td>{{name|safe}}</td>
                                          {% endif %}
                                      {%endwith%}
                                  {% endif %}


                              {%endwith%}
                              {% with locus2interpro|keyvalue:values.2 as interpro_data %}
                              <td>
                                  {% for one_interpro in interpro_data %}
                                  {% if one_interpro.0 != '-' %}
                                      <a href="{% url 'fam' one_interpro.0 'interpro' %}" target="_top">{{one_interpro.0}}</a> /
                                  {%else%}
                                      {{one_interpro.0}}
                                  {%endif%}
                                      {{one_interpro.1}} </br>
                                  {% endfor %}
                                  {%endwith%}
                              </td>
                              <td>{{values.7}}</td>
                              <td>{{values.8}}</td>
                              <td>{{values.9}}</td>
                              <td>{{values.10}}</td>
                              <td>{{values.11}}</td>
                              <td>{{values.12}}</td>

                          </tr>

                          {%endfor%}
                      </tbody>
                  </table>



              </div>

          </div>
        </div>
        </div>




        ''')

        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html




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
          
        sql_phylum = 'select phylum from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
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
              ' where t1.hit_number=1 and t3.name="%s" and t4.phylum!="%s" and t1.query_taxon_id=%s;' % (biodb,
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
                  ' where t1.hit_number=2 and t3.name="%s" and t4.phylum!="%s" and t1.query_taxon_id=%s;' % (biodb,
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



@shared_task
def plot_neighborhood_task(biodb, target_locus, region_size):

    from chlamdb.phylo_tree_display import ete_motifs
    from tempfile import NamedTemporaryFile
    from chlamdb.biosqldb import mysqldb_plot_genomic_feature

    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 3,
                                    'percent': 25,
                                    'description': "Plotting region"})

 
    server, db = manipulate_biosqldb.load_db(biodb)

    sql2 = 'select orthogroup, taxon_id from orthology_detail_%s where locus_tag = "%s"' % (biodb, target_locus)

    reference_orthogroup = server.adaptor.execute_and_fetchall(sql2, )[0]
    reference_taxid = reference_orthogroup[1]

    if reference_orthogroup:
        
        orthogroup = reference_orthogroup[0]
        reference_taxon_id = reference_orthogroup[1]
        

        operon_locus = []
        operon_ofs = False
        temp_location = os.path.join(settings.BASE_DIR, "assets/temp/")
        temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
        name = 'temp/' + os.path.basename(temp_file.name)

        sql10 = 'select operon_id from custom_tables.locus2seqfeature_id_%s t1 ' \
                ' inner join custom_tables.DOOR2_operons_%s t2 on t1.seqfeature_id=t2.seqfeature_id' \
                ' where t1.locus_tag="%s"' % (biodb,
                                                biodb,
                                                target_locus)
        try:
            operon_id = server.adaptor.execute_and_fetchall(sql10, )[0][0]
            sqlo = 'select operon_id,gi,locus_tag,old_locus_tag,COG_number,product from custom_tables.DOOR2_operons_%s t1 ' \
                    ' left join custom_tables.locus2seqfeature_id_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                    ' where operon_id=%s;' % (biodb, biodb, operon_id)
            operon = server.adaptor.execute_and_fetchall(sqlo, )
            operon_locus = [i[2] for i in operon]
        except IndexError:
            try:
                sqlo = 'select C.locus_tag' \
                        ' from (select operon_id from custom_tables.locus2seqfeature_id_%s t1 ' \
                        ' inner join custom_tables.ofs_operons_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                        ' where t1.locus_tag="%s") A ' \
                        ' inner join custom_tables.ofs_operons_%s B on A.operon_id=B.operon_id ' \
                        ' inner join custom_tables.locus2seqfeature_id_%s C on B.seqfeature_id=C.seqfeature_id' % (biodb,
                                                                                            biodb,
                                                                                            target_locus,
                                                                                            biodb,
                                                                                            biodb)



                operon_ofs = server.adaptor.execute_and_fetchall(sqlo, )
                operon_locus = [i[0] for i in operon_ofs]

            except:
                operon_locus = []
                operon_ofs = False

            locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, 
                                                                                            db, 
                                                                                            biodb, 
                                                                                            [target_locus],
                                                                                            temp_file.name, 
                                                                                            int(region_size),
                                                                                            cache, 
                                                                                            operon_locus)


    current_task.update_state(state='PROGRESS',
    meta={'current': 2,
            'total': 3,
            'percent': 50,
            'description': "Preparing profile data"})
    
    taxon2orthogroup2count_all = ete_motifs.get_taxon2orthogroup2count(biodb, orthogroup_list)

    labels = orthogroup_list

    n_orthogroup = orthogroup_list.index(orthogroup)

    taxon2locus2identity_closest = ete_motifs.get_locus2taxon2identity(biodb, locus_tags)

    locus2taxon = {}
    for i in locus_tags:
        locus2taxon[i] = reference_taxid

    # remove potential pseudogenes from the list
    locus_tags_labels = []
    for i in locus_tags:
        if i in locus2taxon.keys():
            locus_tags_labels.append(i)

    current_task.update_state(state='PROGRESS',
    meta={'current': 3,
            'total': 3,
            'percent': 90,
            'description': "Plotting profile"})

    tree2, style = ete_motifs.multiple_profiles_heatmap(biodb,
                                                locus_tags_labels,
                                                taxon2locus2identity_closest,
                                                identity_scale=True,
                                                show_labels=False,
                                                reference_taxon=locus2taxon)


    path = settings.BASE_DIR + '/assets/temp/region.svg'
    asset_path = '/temp/region.svg'
    tree2.render(path, dpi=800, h=600, tree_style=style)
 
 
    template = Template('''
            {% load staticfiles %}
            {% load static %}



            <div class="row">
            <h3>Region</h3>
              <a download="genomic_region.svg" class="btn" href="{% static name %}"><i class="fa fa-download"></i> Download SVG</a>
              <a onclick='exportPNG("genomic_region", "genomic_region");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>

                <object width="105%"  type="image/svg+xml" style="margin-left: 10px" data="{% static name %}" id="genomic_region"></object>
            </div>

            <div class="row">
              <h3>Conservation profile</h3>
              <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
              <a onclick='exportPNG("profile", "profile");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>

                 <object type="image/svg+xml" data="{% static asset_path %}" id="profile"></object>
                 <!--<object width="98%" height="80%" align="left" style="margin-left: 20px" type="image/svg+xml" data="{% static asset_path %}" id="profile"></object>-->
            </div>

            ''')

    html = template.render(Context(locals()))
    return html



@shared_task
def TM_tree_task(biodb, 
                 orthogroup):

    from chlamdb.phylo_tree_display import ete_motifs
    from tempfile import NamedTemporaryFile
    from chlamdb.biosqldb import mysqldb_plot_genomic_feature
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.phylo_tree_display import ete_motifs
    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 50,
                                    'description': "Plotting TM tree"})





    server, db = manipulate_biosqldb.load_db(biodb)

    home_dir = os.path.dirname(__file__)
    alignment_fasta = "../assets/%s_fasta/%s.fa" % (biodb, orthogroup)
    alignment_path = os.path.join(home_dir, alignment_fasta)
    if os.path.exists(alignment_path):
        locus2TM_data = ete_motifs.get_TM_data(biodb, orthogroup, aa_alignment=False, signal_peptide=True)
    else:
        locus2TM_data = ete_motifs.get_TM_data(biodb, orthogroup, aa_alignment=False, signal_peptide=True)
    sql_tree = 'select phylogeny from biosqldb_phylogenies.%s where orthogroup="%s"' % (biodb, orthogroup)
    tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]

    t, ts, leaf_number = ete_motifs.draw_TM_tree(tree, locus2TM_data)
    path = settings.BASE_DIR + '/assets/temp/TM_tree.svg'
    asset_path = '/temp/TM_tree.svg'
    if leaf_number < 10:
        leaf_number = 10
    
    t.render(path, dpi=500, tree_style=ts)


    template = Template('''
            {% load staticfiles %}
            {% load static %}

            {% if not no_tree %}
            <h3>Phylogeny:</h3>
        
        
            <div class="panel panel-success" style="width:80%; top: 200px; margin: 10px 10px 10px 10px">
                <div class="panel-heading" style="width:100%">
                    <h3 class="panel-title">Help</h3>
                </div>
                <p style="margin: 10px 10px 10px 10px">This phylogeny includes all orthologs identified with <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a>. 
                The first column at the right of the phylogeny reports the locus tag of each protein sequence (it can be used in the search bar to search for the corresponding locus). 
                The right part is a representation of the amino acid sequence. The length of the line reflects the length of the sequence. 
                <font color="green">Green blocks</font> are predicted <font color="green">transmembrane domains</font> and <font color="red">red blocks</font> are predicted <font color="red">signal peptides</font>. Both were identified with <a href="http://phobius.sbc.su.se/">Phobius</a>.<br>
                </p>
            </div>
        
            <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
            <a onclick='exportPNG("pfam_tree", "pfam_tree");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>
        
            <div id="pfam_tree_div">
                <object type="image/svg+xml" data="{% static asset_path %}" id="pfam_tree" style="width:90%"></object>
            </div>

            {% else %}
            No tree for {{orthogroup}}
            {% endif %}

            ''')

    html = template.render(Context(locals()))
    return html


@shared_task
def pfam_tree_task(biodb, 
                 orthogroup):

    from chlamdb.phylo_tree_display import ete_motifs
    from tempfile import NamedTemporaryFile
    from chlamdb.biosqldb import mysqldb_plot_genomic_feature
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.phylo_tree_display import ete_motifs
    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 50,
                                    'description': "Plotting TM tree"})

    print ('pfam tree %s -- %s' % (biodb, orthogroup))
    server, db = manipulate_biosqldb.load_db(biodb)

    #sql_locus2protein_id = 'select locus_tag, protein_id from orthology_detail_%s where orthogroup="%s"' % (biodb, orthogroup)

    #locus2protein_id= manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_locus2protein_id,))

    alignment_fasta = "../assets/%s_fasta/%s.fa" % (biodb, orthogroup)

    home_dir = os.path.dirname(__file__)

    alignment_path = os.path.join(home_dir, alignment_fasta)
    print("get data")
    if os.path.exists(alignment_path):
        #pass
        locus2pfam_data = ete_motifs.get_pfam_data(orthogroup, biodb, aa_alignment=False) # alignment_path
    else:

        locus2pfam_data = ete_motifs.get_pfam_data(orthogroup, biodb, aa_alignment=False)
    print("done")
    motif_count = {}
    for data in locus2pfam_data.values():
        for motif in data:
            try:
                if motif[4] not in motif_count:
                    motif_count[motif[4]] = [1, motif[5]]
                else:
                    motif_count[motif[4]][0]+=1
            except:
                print ("motif", motif)

    print("get tree")
    sql_tree = 'select phylogeny from biosqldb_phylogenies.%s where orthogroup="%s"' % (biodb, orthogroup)

    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        no_tree = True
        return render(request, 'chlamdb/pfam_tree.html', locals())


    #sql = 'select taxon_id, family from genomes_classification;'

    #taxon_id2family = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    print("draw tree")
    t, ts, leaf_number = ete_motifs.draw_pfam_tree(tree, locus2pfam_data, False, taxon_id2family=False)
    path = settings.BASE_DIR + '/assets/temp/pfam_tree.svg'
    asset_path = '/temp/pfam_tree.svg'
    if leaf_number < 10:
        leaf_number = 10
    ts.show_branch_support = False
    t.render(path, dpi=500, tree_style=ts)


    template = Template('''
            {% load staticfiles %}
            {% load static %}

            {% if not no_tree %}
            <h3>Phylogeny:</h3>
        
        
            <div class="panel panel-success" style="width:80%; top: 200px; margin: 10px 10px 10px 10px">
                <div class="panel-heading" style="width:100%">
                    <h3 class="panel-title">Help</h3>
                </div>
                <p style="margin: 10px 10px 10px 10px">This phylogeny includes all orthologs identified with <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a>. 
                The first column at the right of the phylogeny reports the locus tag of each protein sequence (that can be used in the search bar to search for the corresponding locus). 
                The right part is a representation of the amino acid sequence. The length of the line reflects the length of the sequence, and <font color="green">green blocks</font> are identified <a href="https://pfam.xfam.org/">PFAM domains</a>.<br>
                </p>
            </div>
        
            <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
            <a onclick='exportPNG("pfam_tree", "pfam_tree");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>
        
            <div id="pfam_tree_div">
                <object type="image/svg+xml" data="{% static asset_path %}" id="pfam_tree" style="width:90%"></object>
            </div>
            <div  id="pfam_table_div">
                <h3>Motifs:</h3>
        
        
                <table class="table">
                    <thead>
                    <tr>
                        <th>Pfam domain</th>
                        <th>Count</th>
                        <th>Description</th>
                    </tr>
                    </thead>
                    <tbody>
                    {% for key, values in motif_count.items %}
                    <tr>
                        <td><a href="{% url 'fam' key "pfam" %}" target="_top">{{key}}</a></td>
                        <td>{{values.0}}</td>
                        <td>{{values.1}}</td>
                    </tr>
                    {% endfor %}
                    </tbody>
                </table>
            </div>
            {% else %}
            No tree for {{orthogroup}}
            {% endif %}
            ''')

    html = template.render(Context(locals()))
    return html



@shared_task
def phylogeny_task(biodb, 
                   orthogroup):


    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 50,
                                    'description': "Plotting TM tree"})


    from chlamdb.biosqldb import manipulate_biosqldb
    from ete3 import Tree
    import os
    from chlamdb.plots import orthogroup2phylogeny_best_refseq_uniprot_hity

    sqlpsw = os.environ['SQLPSW']

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select phylogeny from biosqldb_phylogenies.BBH_%s where orthogroup="%s";' % (biodb, orthogroup)

    ete3_tree = Tree(server.adaptor.execute_and_fetchall(sql,)[0][0])

    t, ts = orthogroup2phylogeny_best_refseq_uniprot_hity.plot_tree(ete3_tree,
                                                                    orthogroup,
                                                                    biodb,
                                                                    mysql_pwd=sqlpsw)
    path = settings.BASE_DIR + '/assets/temp/BBH_tree.svg'
    asset_path = '/temp/BBH_tree.svg'

    t.render(path, tree_style=ts)


    template = Template('''
            {% load staticfiles %}
            {% load static %}
            {% if not no_tree %}
            <h3>Phylogeny including all orthologs and their closest homologs in RefSeq and SwissProt database</h3>
            <p>
              The alignment was made with mafft and the phylogeny was reconstructed with FastTree with default parameters.
            </p>
            <div id="pfam_tree_div">
                <object type="image/svg+xml" data="{% static asset_path %}" id="pfam_tree" style="width:80%"></object>
            </div>
            {% else %}
            No tree for {{orthogroup}}
            {% endif %}
            ''')

    html = template.render(Context(locals()))
    return html


@shared_task
def plot_heatmap_task(biodb,
                      targets, 
                      accessions,
                      type):

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.plots import heatmap
    import numpy as np
    import time
        
    server, db = manipulate_biosqldb.load_db(biodb)
    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 2,
                                    'percent': 50,
                                    'description': "Get data"})

    # particularity of orthology table
    if type == 'orthology':
        col_id = 'orthogroup'
    else:
        col_id = 'id'

    if not accessions:
        # get sub matrix and complete matrix
        taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                    type,
                                                                    col_id,
                                                                    targets,
                                                                    [],
                                                                    ratio=1/float(len(targets)),
                                                                    single_copy=False,
                                                                    accessions=accessions,
                                                                        cache=cache)
        taxon_list = [i.split("_")[1] for i in list(mat.columns.values)]
        labels = [taxon2description[i] for i in taxon_list]
        print(labels)
    else:
        accession2description = manipulate_biosqldb.accession2description(server,biodb)
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                    type,
                                                                    'id',
                                                                    targets,
                                                                    [],
                                                                    ratio=1/float(len(targets)),
                                                                    single_copy=False,
                                                                    accessions=accessions,
                                                                        cache=cache)
        accession_list = list(mat.columns.values)
        labels = [accession2description[i] for i in accession_list]

    m = np.array(mat.transpose())
    m = m.astype(float)
    collabels = [""]*len(m[1,:])
    for i in range(0,len(m[1,:]), 100):
        collabels[i] = i

    path = settings.BASE_DIR + '/assets/temp/heatmap_%s.png' % type
    asset_path = '/temp/heatmap_%s.png' % type

    current_task.update_state(state='PROGRESS',
                              meta={'current': 2,
                                    'total': 2,
                                    'percent': 75,
                                    'description': "Plotting"})

    heatmap.heatmap_pangenome(m,
                                output=path,
                    breaks="-0.5, 0.5, 1.5, 2.5",
                    rows=labels,
                    format="png",
                    orderCols=True)

    named_tuple = time.localtime()
    timestamp = time.strftime("%H%M%S", named_tuple)
    
    template = Template('''
            {% load staticfiles %}
            {% load static %}
            {% if not no_tree %}
            <h3>Heatmap</h3>
            <div id="heatmap">
                <img type="image/svg+xml" src="{% static asset_path %}?{{timestamp}}" alt="heatmap" style="width:90%;">
            </div>
            {% else %}
            No tree for {{orthogroup}}
            {% endif %}
            ''')

    html = template.render(Context(locals()))
    return html
