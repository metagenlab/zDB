
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
def extract_interpro_task(biodb, 
                          include,
                          exclude,
                          freq_missing,
                          reference_taxon,
                          accessions,
                          n_missing):

    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.views import get_locus_annotations
    
    server, db = manipulate_biosqldb.load_db(biodb)
    
    
    current_task.update_state(state='PROGRESS',
    meta={'current': 1,
          'total': 3,
          'percent': 25,
          'description': "Get InterPro entry list"})
    
    if int(n_missing)>=len(include):
        template = Template('''
        
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

    # get sub matrix and complete matrix
    if not accessions:
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                    "interpro",
                                                                    "id",
                                                                    include,
                                                                    exclude,
                                                                    freq_missing,
                                                                    cache=cache)
    else:
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                  "interpro",
                                                                  "id",
                                                                  include,
                                                                  exclude,
                                                                  freq_missing,
                                                                  accessions=accessions,
                                                                  cache=cache)


    match_groups = mat.index.tolist()

    if len(match_groups) == 0:
        template = Template('''
        
        {% load static %}
        <div class="panel panel-warning" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
        <div class="panel-heading" style="width:100%">
            <h3 class="panel-title"></h3>
        </div>
        <p style="margin: 10px 10px 10px 10px">No match</p>
        </div>
        ''')

        html = template.render(Context(locals()))#render_to_string(template, context=locals())
 
        current_task.update_state(state='SECCESS',
                              meta={'current': 3,
                                    'total': 3,
                                    'percent': 100,
                                    'description': "Done"})
 
        return html
    

    # get count in subgroup
    interpro2count = dict((mat > 0).sum(axis=1))
    # get count in complete database
    interpro2count_all = dict((mat_all > 0).sum(axis=1))

    max_n = max(list(interpro2count_all.values()))

    # GET max frequency for template
    sum_group = len(match_groups)

    current_task.update_state(state='PROGRESS',
    meta={'current': 2,
            'total': 3,
            'percent': 50,
            'description': "Retrieve annotations"})

    filter = '"' + '","'.join(match_groups) + '"'

    sql2 = 'select interpro_accession, interpro_description from interpro' \
    ' where interpro_accession in (%s) group by interpro_accession,interpro_description;' % (filter)

    raw_data = list(server.adaptor.execute_and_fetchall(sql2,))

    match_data = []
    for one_match in raw_data:
        match_data.append(list(one_match)+[interpro2count[one_match[0]], interpro2count_all[one_match[0]]])

    sql_include = 'taxon_id ='
    for i in range(0, len(include)-1):
        sql_include+='%s or taxon_id =' % include[i]
    sql_include+=include[-1]
    n = 1
    search_result = []

    group_filter = 'where ('

    for i, group in enumerate(match_groups):
        if i == 0:
            group_filter += 'orthogroup="%s"' % group
        else:
            group_filter += ' or orthogroup="%s"' % group
    group_filter+=')'

    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
    sql_2 = 'select %s from orthology_detail %s' % (columns, group_filter)

    raw_data = server.adaptor.execute_and_fetchall(sql_2,)

    n = 1
    extract_result = []
    for one_hit in raw_data:
        extract_result.append((n,) + one_hit)
        n+=1


    interpro_list = '"' + '","'.join(match_groups) + '"'

    locus_list_sql = 'select locus_tag from interpro where taxon_id=%s ' \
                    ' and interpro_accession in (%s)' % (reference_taxon, interpro_list)

    locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]

    circos_url = '?ref=%s&' % reference_taxon
    circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)
    
    target_circos_taxons = include + exclude

    envoi_extract = True
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
            
            {% load static %}
            {% load custom_tags %}

            <ul id="tabs" class="nav nav-tabs" data-tabs="tabs">
                <li class="active"><a href="#red" data-toggle="tab">Orthogroups</a></li>
                <li><a href="#tab2" data-toggle="tab">Locus list reference</a></li>
            </ul>

            <div id="my-tab-content" class="tab-content">

                <div class="tab-pane active" id="red">


                    <div id="export_bouttons_groups">
                    <form name="circos_form" id="circos_form" action="{% url 'circos_main' %}" method="post">
                    {% csrf_token %}
                        <input type="hidden" name="reference_taxon" value="{{ reference_taxon }}">
                        <input type="hidden" name="target_list" value="{{ target_circos_taxons }}">
                        <input type="hidden" name="highlight" value="{{ locus_list }}">
                        <button class="btn btn-success">Show on circular map</button>
                    </form>
                        <br/>
                    </div>

                    <p>Number of Interpro ID: {{ sum_group }} </p>

                    <table class="display" id="interpro_table">

                        <thead>
                            <tr>
                                <th>Interpro accession</th>
                                <th>Description</th>
                                <th>Count (/{{include|length}})</th><th>Count all (/{{max_n}})</th>
                            </tr>
                        </thead>

                        {% for values in match_data %}
                            <tr>
                                <td><a href="{% url 'fam' values.0 'interpro' %}">{{values.0}}</a></td>
                                <td>{{values.1}}</td>
                                <td>{{values.2}}</td>
                                <td>{{values.3}}</td>
                            </tr>
                        {% endfor %}
                    </table>
                </div>
                <div class="tab-pane" id="orange">

                </div>

                <div class="tab-pane" id="tab2">

                    <h3> Locus annotation </h3>
                    <div id="export_bouttons_groups">
                        <a href="{% url 'get_fasta' %}{{fasta_url_ref}}" class="btn btn-success">Download fasta</a>
                    <form name="circos_form" id="circos_form" action="{% url 'circos_main' %}" method="post">
                    {% csrf_token %}
                        <input type="hidden" name="reference_taxon" value="{{ reference_taxon }}">
                        <input type="hidden" name="target_list" value="{{ target_circos_taxons }}">
                        <input type="hidden" name="highlight" value="{{ locus_list }}">
                        <button class="btn btn-success">Show on circular map</button>
                    </form>

                        <br/>
                    </div>
                    <table id="reference_genome" class="display">
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
                                <td><a href="http://www.ncbi.nlm.nih.gov/protein/{{values.3}}" target="_top">{{locus_tag2cog_catego|keyvalue:values.2}}</a></td>
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



            ''')

    html = template.render(Context(locals()))#render_to_string(template, context=locals())

    current_task.update_state(state='SECCESS',
                            meta={'current': 3,
                                'total': 3,
                                'percent': 100,
                                'description': "Done"})

    return html



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
    
    server, db = manipulate_biosqldb.load_db(biodb)
    
    
    current_task.update_state(state='PROGRESS',
    meta={'current': 1,
          'total': 3,
          'percent': 25,
          'description': "Get orthogroup list"})
    
    if int(n_missing)>=len(include):
        template = Template('''
        
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

    match_groups = mat.index.tolist()

    if len(match_groups) == 0:
        template = Template('''
        
        {% load static %}
        <div class="panel panel-warning" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
        <div class="panel-heading" style="width:100%">
            <h3 class="panel-title"></h3>
        </div>
        <p style="margin: 10px 10px 10px 10px">No match</p>
        </div>
        ''')

        html = template.render(Context(locals()))#render_to_string(template, context=locals())
 
        current_task.update_state(state='SECCESS',
                              meta={'current': 3,
                                    'total': 3,
                                    'percent': 100,
                                    'description': "Done"})
 
        return html
    else:

        # get count in subgroup
        orthogroup2count = dict((mat > 0).sum(axis=1))
        # get count in complete database
        orthogroup2count_all = dict((mat_all > 0).sum(axis=1))

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

        if len(include) == 1:
            # get url to get single include taxon fasta
            locus_list = [i[2] for i in extract_result]
            fasta_ref_url = '?l=' + '&l='.join(locus_list)

        columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                  'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'

        envoi_extract = True

        target_circos_taxons = include + exclude

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
            sql = 'select locus_tag from orthology_detail where orthogroup in (%s) and taxon_id=%s' % ('"' + '","'.join(match_groups) + '"',
                                                                                                        reference_taxon)
        else:
            sql = 'select locus_tag from orthology_detail where orthogroup in (%s) and accession="%s"' % ('"' + '","'.join(match_groups) + '"',
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
                        
                        <form name="circos_form" id="circos_form" action="{% url 'circos_main' %}" method="post">
                            <input type="hidden" name="reference_taxon" value="{{ reference_taxon }}">
                            <input type="hidden" name="target_list" value="{{ target_circos_taxons }}">
                            <input type="hidden" name="highlight" value="{{ locus_list }}">
                            <button class="btn btn-success">Show on circular map</button>
                        </form>
                        
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
                      <button class="btn btn-success">Show on circular map</button>

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
 
        current_task.update_state(state='SECCESS',
                              meta={'current': 3,
                                    'total': 3,
                                    'percent': 100,
                                    'description': "Done"})
 
        return html




@shared_task
def run_circos(reference_taxon, target_taxons):

    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos
    from chlamdb.biosqldb import shell_command
    import ete3
    from chlamdb.plots import gbk2circos
    from tempfile import NamedTemporaryFile

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
        biorecord = cache.get(biodb + "_" + accession)

        if not biorecord:
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
            record_list.append(biorecord)


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

    sql_order = 'select A.taxon_1 from (select taxon_1, median_identity from comparative_tables_shared_og_av_id where taxon_2=%s ' \
                ' union select taxon_2,median_identity from comparative_tables_shared_og_av_id ' \
                ' where taxon_1=%s order by median_identity DESC) A;' % (reference_taxon,
                                                                         reference_taxon)
    try:
        sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa where taxon_1=%s order by identity desc;' % (reference_taxon)
        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]

    except:

        # median identity
        sql_order = 'select taxon from (select taxon_2 as taxon, median_identity ' \
                    ' from comparative_tables_shared_og_av_id where taxon_1=%s union ' \
                    ' select taxon_1, median_identity as taxon from comparative_tables_shared_og_av_id' \
                    '  where taxon_2=%s) A order by median_identity desc;' % (reference_taxon,
                                                                              reference_taxon)
        '''
        # n shared orthogroups
        sql_order = 'select taxon_2 from comparative_tables_shared_orthogroups where taxon_1=%s order by n_shared_orthogroups DESC;' % (biodb,
                                                                                                                  reference_taxon)
        '''
        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]

    '''

    t1 = ete3.Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)

    node_list = []
    for node in t1.iter_leaves():
            node_list.append(node.name)

    reference_index = node_list.index(reference_taxon)
    ordered_taxons = node_list[reference_index:] + node_list[:reference_index][::-1]
    '''

    # TODO add this option to view
    highlight_BBH = False
    
    if highlight_BBH:
          
        sql_phylum = 'select phylum from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                     ' inner join taxid2species t3 on t2.taxon_id=t3.taxon_id ' \
                     ' inner join species_curated_taxonomy t4 on t3.species_id=t4.species_id ' \
                     ' where t1.name="%s" and t2.taxon_id=%s limit 1; ' % (biodb, 
                                                                           reference_taxon)
        reference_phylum = server.adaptor.execute_and_fetchall(sql_phylum,)[0][0]
        
        try:
            sql = 'select locus_tag from blastnr_blastnr t1 ' \
              ' inner join bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
              ' inner join custom_tables_locus2seqfeature_id t5 ' \
              ' on t1.seqfeature_id=t5.seqfeature_id ' \
              ' where t1.hit_number=1 and t3.name="%s" and t4.phylum!="%s" and t1.query_taxon_id=%s;' % (biodb,
                                                                                                         reference_phylum,
                                                                                                         reference_taxon)
            BBH_color = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            sql2 = 'select locus_tag from custom_tables_locus2seqfeature_id t1 ' \
                   ' left join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                   ' where taxon_id=%s and t2.seqfeature_id is NULL;' % (reference_taxon)
                   
            no_BBH_hit_color = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]
            
            if len(BBH_color) < 20:
                sql = 'select locus_tag from blastnr_blastnr t1 ' \
                  ' inner join bioentry t2 on t1.query_bioentry_id=t2.bioentry_id ' \
                  ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                  ' inner join blastnr_blastnr_taxonomy t4 on t1.subject_taxid=t4.taxon_id ' \
                  ' inner join custom_tables_locus2seqfeature_id t5 ' \
                  ' on t1.seqfeature_id=t5.seqfeature_id ' \
                  ' where t1.hit_number=2 and t3.name="%s" and t4.phylum!="%s" and t1.query_taxon_id=%s;' % (biodb,
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


    temp_location = os.path.join(settings.BASE_DIR, "assets/circos/")
    temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
    circos_file = 'circos/' + os.path.basename(temp_file.name)
    base_file_name = os.path.basename(temp_file.name).split(".svg")[0]

    circos_svg_file = "circos/%s.svg" % base_file_name
    circos_png_file = "circos/%s.png" % base_file_name
    original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % base_file_name
    map_file = "circos/%s.html" % base_file_name
    svg_file = "circos/%s.svg" % base_file_name
    map_name = base_file_name

    original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % base_file_name

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
                                              ordered_taxons=ordered_taxons,
                                              outfile_prefix=base_file_name)

    
    
    
    with open(original_map_file, "r") as f:
        map_string = ''.join([line for line in f.readlines()])

    circos_html = '<!DOCTYPE html>\n' \
                  ' <html>\n' \
                  ' <body>\n' \
                  ' %s\n' \
                  ' <img src="%s" usemap="#%s">' \
                  ' </body>\n' \
                  ' </html>\n' % (map_string, os.path.basename(temp_file.name), base_file_name)

    circos_new_file = '/assets/circos/circos_clic_%s.html' % base_file_name


    current_task.update_state(state='PROGRESS',
                              meta={'current': 4,
                                    'total': 4,
                                    'percent': 100,
                                    'description': "Rendering..."})

    with open(settings.BASE_DIR + circos_new_file, "w") as f:
        f.write(circos_html)


    template = Template('''
            
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

    current_task.update_state(state='SECCESS',
                              meta={'current': 4,
                                    'total': 4,
                                    'percent': 100,
                                    'description': "Done"})

    return html


@shared_task
def run_circos_main(reference_taxon, target_taxons, highlight):


    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos
    from chlamdb.biosqldb import shell_command
    import ete3
    from chlamdb.plots import gbk2circos
    from tempfile import NamedTemporaryFile

    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)


    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 4,
                                    'percent': 25,
                                    'description': "Get reference genome"})

    reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb) # ["NC_009648"] NZ_CP009208 NC_016845

    record_list = []
    for accession in reference_accessions:
        biorecord = cache.get(biodb + "_" + accession)

        if not biorecord:
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

    sql_order1 = 'select A.taxon_1 from (select taxon_1,median_identity from comparative_tables_shared_og_av_id where taxon_2=%s ' \
                ' union select taxon_2,median_identity from comparative_tables_shared_og_av_id ' \
                ' where taxon_1=%s order by median_identity DESC) A;' % (reference_taxon, 
                                                                         reference_taxon)
    try:
        sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa where taxon_1=%s order by identity desc;' % (reference_taxon)

        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
    except:
        sql_order2 = 'select taxon_2 from comparative_tables_shared_og_av_id where taxon_1=%s order by median_identity desc;' % (reference_taxon)

        ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order1)]

    current_task.update_state(state='PROGRESS',
                              meta={'current': 3,
                                    'total': 4,
                                    'percent': 75,
                                    'description': "Plotting with circos..."})


    temp_location = os.path.join(settings.BASE_DIR, "assets/circos/")
    temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
    
    circos_file = 'circos/' + os.path.basename(temp_file.name)
    base_file_name = os.path.basename(temp_file.name).split(".svg")[0]
    original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % base_file_name
    #original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % base_file_name
    #map_file = "circos/%s.html" % base_file_name
    circos_svg_file = "circos/%s.svg" % base_file_name
    circos_png_file = "circos/%s.png" % base_file_name
    #map_name = ref_name
        
    myplot = circos.CircosAccession2multiplot(server,
                              db,
                              biodb,
                              record_list,
                              target_accessions,
                              locus_highlight=highlight,
                              out_directory=temp_location,
                              draft_fasta=draft_data,
                              href="/chlamdb/locusx/",
                              ordered_taxons = ordered_taxons,
                              outfile_prefix=base_file_name)


    # read original map file generated by circos
    with open(original_map_file, "r") as f:
        map_string = ''.join([line for line in f.readlines()])

    circos_html = '<!DOCTYPE html>\n' \
                  ' <html>\n' \
                  ' <body>\n' \
                  ' %s\n' \
                  ' <img src="%s.svg" usemap="#%s">' \
                  ' </body>\n' \
                  ' </html>\n' % (map_string, base_file_name, base_file_name)

    # write a new html file including svg + mapping data
    circos_new_file = '/assets/circos/_%s.html' % base_file_name
    with open(settings.BASE_DIR + circos_new_file, "w") as f:
        f.write(circos_html)

    template = Template('''
            
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

    current_task.update_state(state='SECCESS',
                              meta={'current': 4,
                                    'total': 4,
                                    'percent': 100,
                                    'description': "Done"})

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

    sql2 = 'select orthogroup, taxon_id from orthology_detail where locus_tag = "%s"' % (target_locus)

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

        sql10 = 'select operon_id from custom_tables_locus2seqfeature_id t1 ' \
                ' inner join custom_tables_DOOR2_operons t2 on t1.seqfeature_id=t2.seqfeature_id' \
                ' where t1.locus_tag="%s"' % (target_locus)
        try:
            operon_id = server.adaptor.execute_and_fetchall(sql10, )[0][0]
            sqlo = 'select operon_id,gi,locus_tag,old_locus_tag,COG_number,product from custom_tables_DOOR2_operons t1 ' \
                    ' left join custom_tables_locus2seqfeature_id t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                    ' where operon_id=%s;' % (operon_id)
                    
            operon = server.adaptor.execute_and_fetchall(sqlo, )
            operon_locus = [i[2] for i in operon]
        except:
            try:
                sqlo = 'select C.locus_tag' \
                        ' from (select operon_id from custom_tables_locus2seqfeature_id t1 ' \
                        ' inner join custom_tables_ofs_operons t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                        ' where t1.locus_tag="%s") A ' \
                        ' inner join custom_tables_ofs_operons B on A.operon_id=B.operon_id ' \
                        ' inner join custom_tables_locus2seqfeature_id C on B.seqfeature_id=C.seqfeature_id' % (target_locus)
                        
                operon_ofs = server.adaptor.execute_and_fetchall(sqlo, )
                operon_locus = [i[0] for i in operon_ofs]
            # case when no operon
            # or no operon table
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
    tree2.render(path, dpi=800, w=800, tree_style=style)
 
 
    template = Template('''
            
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

    current_task.update_state(state='SECCESS',
                              meta={'current': 3,
                                    'total': 3,
                                    'percent': 100,
                                    'description': "Done"})

    return html


@shared_task
def basic_tree_task(biodb, 
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

    sql_tree = 'select phylogeny from phylogenies where orthogroup="%s"' % (orthogroup)
    
    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        template = Template('''
        
        {% load static %}
        <div class="panel panel-danger" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
            <div class="panel-heading" style="width:100%">
                <h3 class="panel-title">Warning</h3>
            </div>
            <p style="margin: 10px 10px 10px 10px">No tree</p>
        </div>
        ''')
        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html

    sql = f'select distinct locus_tag,t4.description from orthology_orthogroup t1 ' \
          f' inner join orthology_seqfeature_id2orthogroup t2 on t1.orthogroup_id=t2.orthogroup_id ' \
          f' inner join annotation_seqfeature_id2locus t3 on t2.seqfeature_id=t3.seqfeature_id ' \
          f' inner join bioentry t4 on t3.bioentry_id=t4.bioentry_id where orthogroup_name="{orthogroup}";'
    
    locus2organism = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql, ))

    t, ts, leaf_number = ete_motifs.draw_basic_tree(tree, locus2organism)
    path = settings.BASE_DIR + '/assets/temp/TM_tree.svg'
    asset_path = '/temp/TM_tree.svg'
    if leaf_number < 10:
        leaf_number = 10
    
    t.render(path, dpi=500, tree_style=ts)


    template = Template('''
            
            {% load static %}

            {% if not no_tree %}     
        
            <div class="panel panel-success" style="width:80%; top: 200px; margin: 10px 10px 10px 10px">
                <div class="panel-heading" style="width:100%">
                    <h3 class="panel-title">Help</h3>
                </div>
                <p style="margin: 10px 10px 10px 10px">
                Phylogeny
                </p>
            </div>
        
            <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
            <a onclick='exportPNG("pfam_tree", "pfam_tree");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>
            <a class="btn" href="{% url 'get_newick_tree' orthogroup False %}"><i class="fa fa-download"></i> Download tree (newick format)</a>

            <div id="pfam_tree_div">
                <object type="image/svg+xml" data="{% static asset_path %}" id="pfam_tree" style="width:90%"></object>
            </div>

            {% else %}
            No tree for {{orthogroup}}
            {% endif %}

            ''')

    html = template.render(Context(locals()))

    current_task.update_state(state='SECCESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 100,
                                    'description': "Done"})
   
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
    sql_tree = 'select phylogeny from phylogenies where orthogroup="%s"' % (orthogroup)
    
    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        template = Template('''
        
        {% load static %}
        <div class="panel panel-danger" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
            <div class="panel-heading" style="width:100%">
                <h3 class="panel-title">Warning</h3>
            </div>
            <p style="margin: 10px 10px 10px 10px">No tree</p>
        </div>
        ''')
        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html

    t, ts, leaf_number = ete_motifs.draw_TM_tree(tree, locus2TM_data)
    path = settings.BASE_DIR + '/assets/temp/TM_tree.svg'
    asset_path = '/temp/TM_tree.svg'
    if leaf_number < 10:
        leaf_number = 10
    
    t.render(path, dpi=500, tree_style=ts)


    template = Template('''
            
            {% load static %}

            {% if not no_tree %}    
        
            <div class="panel panel-success" style="width:80%; top: 200px; margin: 10px 10px 10px 10px">
                <div class="panel-heading" style="width:100%">
                    <h3 class="panel-title">Help</h3>
                </div>
                <p style="margin: 10px 10px 10px 10px">This phylogeny includes all orthologs identified with <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a>. 
                (see detailed method <a href="/docs/methods/annotation.html#id1" target="_top">here</a>). 
                The first column at the right of the phylogeny reports the locus tag of each protein sequence (it can be used in the search bar to search for the corresponding locus) 
                The right part is a representation of the amino acid sequence. The length of the line reflects the length of the sequence. 
                <font color="green">Green blocks</font> are predicted <font color="green">transmembrane domains</font> and <font color="red">red blocks</font> are predicted <font color="red">signal peptides</font>. Both were identified with <a href="http://phobius.sbc.su.se/">Phobius</a>.<br>
                </p>
            </div>
        
            <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
            <a onclick='exportPNG("pfam_tree", "pfam_tree");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>
            <a class="btn" href="{% url 'get_newick_tree' orthogroup False %}"><i class="fa fa-download"></i> Download tree (newick format)</a>

            <div id="pfam_tree_div">
                <object type="image/svg+xml" data="{% static asset_path %}" id="pfam_tree" style="width:90%"></object>
            </div>

            {% else %}
            No tree for {{orthogroup}}
            {% endif %}

            ''')

    html = template.render(Context(locals()))

    current_task.update_state(state='SECCESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 100,
                                    'description': "Done"})
   
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

    server, db = manipulate_biosqldb.load_db(biodb)

    #sql_locus2protein_id = 'select locus_tag, protein_id from orthology_detail where orthogroup="%s"' % (biodb, orthogroup)

    #locus2protein_id= manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_locus2protein_id,))

    alignment_fasta = "../assets/%s_fasta/%s.fa" % (biodb, orthogroup)

    home_dir = os.path.dirname(__file__)

    alignment_path = os.path.join(home_dir, alignment_fasta)

    if os.path.exists(alignment_path):
        #pass
        locus2pfam_data = ete_motifs.get_pfam_data(orthogroup, biodb, aa_alignment=False) # alignment_path
    else:

        locus2pfam_data = ete_motifs.get_pfam_data(orthogroup, biodb, aa_alignment=False)

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

    sql_tree = 'select phylogeny from phylogenies where orthogroup="%s"' % (orthogroup)

    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        no_tree = True
        template = Template('''
        
        {% load static %}
        <div class="panel panel-danger" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
            <div class="panel-heading" style="width:100%">
                <h3 class="panel-title">Warning</h3>
            </div>
            <p style="margin: 10px 10px 10px 10px">No tree</p>
        </div>
        ''')
        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html

    #sql = 'select taxon_id, family from genomes_classification;'

    #taxon_id2family = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    t, ts, leaf_number = ete_motifs.draw_pfam_tree(tree, locus2pfam_data, False, taxon_id2family=False)
    path = settings.BASE_DIR + '/assets/temp/pfam_tree.svg'
    asset_path = '/temp/pfam_tree.svg'
    if leaf_number < 10:
        leaf_number = 10
    ts.show_branch_support = False
    t.render(path, dpi=500, tree_style=ts)


    template = Template('''
            
            {% load static %}

            {% if not no_tree %}
             
            <div class="panel panel-success" style="width:80%; top: 200px; margin: 10px 10px 10px 10px">
                <div class="panel-heading" style="width:100%">
                    <h3 class="panel-title">Help</h3>
                </div>
                <p style="margin: 10px 10px 10px 10px">This phylogeny includes all orthologs identified with <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a> 
                (see detailed method <a href="/docs/methods/annotation.html#id1" target="_top">here</a>).
                The first column at the right of the phylogeny reports the locus tag of each protein sequence (that can be used in the search bar to search for the corresponding locus). 
                The right part is a representation of the amino acid sequence. The length of the line reflects the length of the sequence, and <font color="green">green blocks</font> are identified <a href="https://pfam.xfam.org/">PFAM domains</a>.<br>
                </p>
            </div>
        
            <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
            <a onclick='exportPNG("pfam_tree", "pfam_tree");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>
            <a class="btn" href="{% url 'get_newick_tree' orthogroup False %}"><i class="fa fa-download"></i> Download tree (newick format)</a>
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

    current_task.update_state(state='SECCESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 100,
                                    'description': "Done"})

    return html



@shared_task
def phylogeny_task(biodb, 
                   orthogroup):


    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 50,
                                    'description': "Plotting phylogenetic tree"})


    from chlamdb.biosqldb import manipulate_biosqldb
    from ete3 import Tree
    import os
    from chlamdb.plots import orthogroup2phylogeny_best_refseq_uniprot_hity

    sqlpsw = os.environ['SQLPSW']

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select phylogeny from phylogenies_BBH where orthogroup="%s";' % (orthogroup)

    ete3_tree = Tree(server.adaptor.execute_and_fetchall(sql,)[0][0])

    t, ts = orthogroup2phylogeny_best_refseq_uniprot_hity.plot_tree(ete3_tree,
                                                                    orthogroup,
                                                                    biodb,
                                                                    mysql_pwd=sqlpsw)
    path = settings.BASE_DIR + '/assets/temp/BBH_tree.svg'
    asset_path = '/temp/BBH_tree.svg'

    t.render(path, tree_style=ts)


    template = Template('''
            
            {% load static %}
            {% if not no_tree %}
            <h3>Phylogeny including best RefSeq hits</h3>

            <div id="pfam_tree_div">
                <div class="row">

                    <div class="panel panel-success" style="width:80%; top: 200px; margin: 10px 10px 10px 10px">
                        <div class="panel-heading" style="width:100%">
                            <h3 class="panel-title">Help</h3>
                        </div>
                        <p style="margin: 10px 10px 10px 10px">This phylogeny includes all orthologs identified with <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a> as well as the 4 best RefSeq hit of each protein 
                        (see detailed method <a href="/docs/methods/annotation.html#phlogeny-including-top-refseq-hits" target="_top">here</a>). </br>
                        PVC bacteria are coloured in red. Other species are coloured according to their phylum-level classification.
                         </p>
                    </div>

                    <a download="phylogeny.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
                    <a onclick='exportPNG("pyhlogeny", "pyhlogeny");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>
                    <a class="btn" href="{% url 'get_newick_tree' orthogroup True %}"><i class="fa fa-download"></i> Download tree (newick format)</a>
                </div>
                <div class="row">
                    <object type="image/svg+xml" data="{% static asset_path %}" id="pyhlogeny" style="width:80%"></object>
                </div>
            </div>
            {% else %}
            No tree for {{orthogroup}}
            {% endif %}
            ''')

    html = template.render(Context(locals()))
    
    current_task.update_state(state='SUCCESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 100,
                                    'description': "Plotting phylogenetic tree"})
    
    return html


@shared_task
def plot_heatmap_task(biodb, targets, type):
    from chlamdb.plots import heatmap
    import numpy as np
    import time

    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 2,
                                    'percent': 50,
                                    'description': "Get data"})

    # Not too happy with this solution. Really need a database connection
    # manager.
    db = db_utils.DB.load_db_from_name(biodb)
    # particularity of orthology table
    if type == 'orthology':
        col_id = 'orthogroup'
    else:
        col_id = 'id'

    mat = db.get_cog_hits(targets, only_best_hits=True, as_count=True)
    target2description = db.get_genomes_description(targets)
    bioentries_list = [i for i in list(mat.columns.values)]
    labels = [target2description[i] for i in bioentries_list]

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

    current_task.update_state(state='SUCCESS',
                              meta={'current': 2,
                                    'total': 2,
                                    'percent': 100,
                                    'description': "Done"})
    return html


@shared_task
def KEGG_map_ko_task(biodb, 
                     map_name):
    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 50,
                                    'description': "Preparing data"})


    from chlamdb.phylo_tree_display import ete_motifs
    from chlamdb.plots import kegg_maps

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select pathway_name,pathway_category,description,C.EC, C.ko_accession, D.definition, A.pathway_id from ' \
            ' (select * from enzyme_kegg_pathway where pathway_name="%s") A inner join enzyme_pathway2ko as B ' \
            ' on A.pathway_id=B.pathway_id inner join enzyme_ko_annotation as C on B.ko_id=C.ko_id ' \
            ' inner join enzyme_ko_annotation as D on B.ko_id=D.ko_id;' % (map_name)

    map_data = server.adaptor.execute_and_fetchall(sql,)
    try:
        pathway_id = map_data[0][-1]
    except IndexError:
        template = Template('''
        
        {% load static %}
        <div class="panel panel-danger" style="width:500px ; top: 200px; margin: 10px 10px 10px 10px">
            <div class="panel-heading" style="width:100%">
                <h3 class="panel-title">Warning</h3>
            </div>
            <p style="margin: 10px 10px 10px 10px">Invalid map accession</p>
        </div>
        ''')

        html = template.render(Context(locals()))#render_to_string(template, context=locals())
        return html
    
    ko_list = [i[4] for i in map_data]

    # fetch ko pylogenetic profiles
    sql = 'select id from comparative_tables_ko where id in (%s);' % ('"' + '","'.join(ko_list) + '"')

    ko_list_found_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    # get list of all orthogroups with corresponding KO
    sql = 'select distinct aa.ko_accession, bb.orthogroup from ' \
            ' (select C.ko_accession, E.seqfeature_id from   ' \
            ' (select * from enzyme_kegg_pathway  where pathway_id=%s) A ' \
            ' inner join enzyme_pathway2ko as B  on A.pathway_id=B.pathway_id  ' \
            ' inner join enzyme_ko_annotation as C on B.ko_id=C.ko_id  ' \
            ' inner join enzyme_seqfeature_id2ko E on B.ko_id=E.ko_id ' \
            ' group by B.ko_id, E.seqfeature_id) aa ' \
            ' inner join orthology_detail bb on aa.seqfeature_id=bb.seqfeature_id;' % (pathway_id)

    orthogroup_data = server.adaptor.execute_and_fetchall(sql,)
    ko2orthogroups = {}
    orthogroup_list = []
    for i in orthogroup_data:
        if i[0] not in ko2orthogroups:
            ko2orthogroups[i[0]] = [i[1]]
        else:
            ko2orthogroups[i[0]].append(i[1])
        orthogroup_list.append(i[1])

    taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup")
    taxon2ko2count = ete_motifs.get_taxon2name2count(biodb, ko_list_found_in_db, type="ko")

    labels = ko_list_found_in_db
    tree, style = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2ko2count)

    tree2, style2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                    labels,
                                                    taxon2orthogroup2count,
                                                    taxon2ko2count,
                                                    ko2orthogroups)


    if len(labels) > 70:
        big = True
        path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.png' % map_name
        asset_path = '/temp/KEGG_tree_%s.png' % map_name
        tree.render(path, dpi=1200, tree_style=style)

        path2 = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s_complete.png' % map_name
        asset_path2 = '/temp/KEGG_tree_%s_complete.png' % map_name

        tree2.render(path2, dpi=800, tree_style=style2)

    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.svg' % map_name
        asset_path = '/temp/KEGG_tree_%s.svg' % map_name
        tree.render(path, dpi=800, tree_style=style)

        path2 = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s_complete.svg' % map_name
        asset_path2 = '/temp/KEGG_tree_%s_complete.svg' % map_name
        tree2.render(path2, dpi=800, tree_style=style2)

    sql = 'select bb.ko_accession, count(*) from (select C.ko_accession, E.seqfeature_id from  ' \
            ' (select * from enzyme_kegg_pathway  where pathway_name="%s") A ' \
            ' inner join enzyme_pathway2ko as B  on A.pathway_id=B.pathway_id ' \
            ' inner join enzyme_ko_annotation as C on B.ko_id=C.ko_id ' \
            ' inner join enzyme_seqfeature_id2ko E on B.ko_id=E.ko_id ' \
            ' group by B.ko_id,seqfeature_id) bb group by ko_accession;' % (map_name)

    ko2freq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    path_map_freq = settings.BASE_DIR + '/assets/temp/KEGG_map_freq_%s' % map_name
    path_map_freq_svg = '/temp/KEGG_map_freq_%s.svg' % map_name

    keg_path_highlight = '+'+'+'.join(ko_list_found_in_db)

    kegg_maps.map2highlighted_map(map_name, ko_list_found_in_db, ko2freq, biodb, path_map_freq+'.pdf')

    template = Template('''
                      
                      {% load static %}
                      {% load custom_tags %}

                     <h3> {% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}: {% with map_data|first as first_doc %}{{ first_doc.2 }}{% endwith %} ({% with map_data|first as first_doc %}{{ first_doc.1 }}{% endwith %}) {% if organism%}:  {{organism}} {%endif%}</h3>

                      <ul id="tabs" class="nav nav-tabs" data-tabs="tabs">
                          <li class="active"><a href="#tab1" data-toggle="tab">General</a></li>
                          <li><a href="#tab2" data-toggle="tab">Profile</a></li>
                          <li><a href="#tab3" data-toggle="tab">Profile + homologs</a></li>
                          <li><a href="#tab4" data-toggle="tab">KEGG map</a></li>

                      </ul>

                      <div id="my-tab-content" class="tab-content">
                          <div class="tab-pane active" id="tab1">
                              </br>
                              <ol> 
                                  <li><a href='http://www.genome.jp/kegg-bin/show_pathway?{% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}{{keg_path_highlight}}'>Link to KEGG website <i class="fa fa-external-link"></i></a></li>
                              </ol>

                              <h4>List of Kegg Orthologs associated to KEGG pathways {% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}</h4>

                              <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                                    <div class="panel-heading" style="width:100%">
                                        <h3 class="panel-title">Help</h3>
                                    </div>
                                    <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                        The "KO coccurence" column report the number of proteins annotated with this Kegg Ortholog in the entire database. 
                                        Click on the KO accession to retrieve the detailed list of proteins.
                                    </p>
                                </div>

                              <table id="ko_table">
                                  <thead>
                                  <tr>
                                      <th>EC(s)</th>
                                      <th>KO</th>
                                      <th>KO occurences</th>
                                      <th>Description</th>
                                  </tr>
                                  </thead>
                                  <tbody>
                                  {% for values in map_data%}
                                      <tr>
                                          {%if not ' ' in values.3 and not '-' in values.3 and not ',' in values.3 %}
                                              <td><a href="{% url 'fam'  values.3 'EC' %}" target="_top">{{values.3}}</a></td>
                                          {% else %}
                                              <td>{{values.3}}</td>
                                          {% endif %}
                                          {%if ko2freq|keyvalue:values.4 != None %}
                                             <td><a href="{% url 'fam'  values.4 'ko' %}" target="_top">{{values.4}}</a></td>
                                             <td>{{ko2freq|keyvalue:values.4}}</td>
                                          {%else%}
                                             <td>{{values.4}}</td>
                                             <td>0</td>
                                          {%endif%}

                                          <td>{{values.5}}</td>
                                      </tr>
                                  {% endfor %}

                                  </tbody>
                              </table>
                          </div>
                          <div class="tab-pane" id="tab2">

                                <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                                        <div class="panel-heading" style="width:100%">
                                            <h3 class="panel-title">Help</h3>
                                        </div>
                                        <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                            Conservation of the pathway among species of the PVC superphylum. A blue background indicate that one or multiple 
                                            proteins were annotated with the corresponding Kegg Ortholog using KoFamScan (<a href="https://www.genome.jp/tools/kofamkoala/">https://www.genome.jp/tools/kofamkoala/</a>).
                                        </p>
                                    </div>

                                <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
                                <a onclick='exportPNG("profile", "profile");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>

                              <object type="image/svg+xml" data="{% static asset_path %}" id="profile" style="max-width:100%"></object>
                          </div>
                          <div class="tab-pane" id="tab3">

                                <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                                        <div class="panel-heading" style="width:100%">
                                            <h3 class="panel-title">Help</h3>
                                        </div>
                                        <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                            Conservation of the pathway among species of the PVC superphylum. A blue background indicate that one or multiple 
                                            proteins were annotated with the corresponding Kegg Ortholog using KoFamScan (<a href="https://www.genome.jp/tools/kofamkoala/">https://www.genome.jp/tools/kofamkoala/</a>).
                                            A green background indicates that no protein were annoated using KoFamScan but an homolog was identified using <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a>. 
                                        </p>
                                </div>

                                <a download="profile_all.svg" class="btn" href="{% static asset_path2 %}"><i class="fa fa-download"></i> Download SVG</a>
                                <a onclick='exportPNG("profile_homologs", "profile_homologs");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>

                              <object type="image/svg+xml" data="{% static asset_path2 %}" id="profile_homologs" style="max-width:100%"></object>
                          </div>
                          <div class="tab-pane" id="tab4">

                           <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                              <div class="panel-heading" style="width:100%">
                                  <h3 class="panel-title">Help</h3>
                              </div>
                              <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                  {% if organism %}
                                      <strong>Red scale:</strong> no homolog identified in the genome of {{organism}}. The color scale reflects the number of homologs in other genomes (the more the darker). <br>
                                      <strong>Green scale:</strong> homolog identified in the genome of {{organism}}. The color scale reflects the number of homologs in other genomes.
                                      <br>
                                  {% else %}
                                      <strong>Green scale:</strong> reflects the conservation of the KO among all organisms included in the database. 
                                      You can click on pathway names and KO accessions to retrieve more information about that step.
                                  {% endif %}

                              </p>
                          </div>

                              <object type="image/svg+xml" data="{% static path_map_freq_svg %}"></object>
                          </div>
                      </div>

                    <script>
                        $(document).ready(function() {
                            $('#ko_table').DataTable( {
                                dom: 'Bfrtip',
                                "order": [[2, "desc" ]],
                                "pageLength": 15,
                                "paging":   true,
                                "ordering": true,
                                "info":     false,
                                buttons: [
                                {
                                    extend: 'excel',
                                    title: '{% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}'
                                },
                                {
                                    extend: 'csv',
                                    title: '{% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}'
                                }
                                ],
                            } );
                        } );
                    </script>

            ''')

    html = template.render(Context(locals()))
    
    current_task.update_state(state='SUCCESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 100,
                                    'description': "Plotting phylogenetic tree"})
    
    return html

@shared_task
def KEGG_map_ko_organism_task(biodb, 
                              map_name,
                              taxon_id):
    
    current_task.update_state(state='PROGRESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 50,
                                    'description': "Preparing data"})

    from chlamdb.phylo_tree_display import ete_motifs
    from chlamdb.plots import kegg_maps

    print ('kegg mapp organism %s -- %s -- %s' % (biodb, map_name, taxon_id))

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select pathway_name,pathway_category,description,C.EC, C.ko_accession, D.definition, A.pathway_id from ' \
            ' (select * from enzyme_kegg_pathway where pathway_name="%s") A inner join enzyme_pathway2ko as B ' \
            ' on A.pathway_id=B.pathway_id inner join enzyme_ko_annotation as C on B.ko_id=C.ko_id ' \
            ' inner join enzyme_ko_annotation as D on B.ko_id=D.ko_id;' % (map_name)

    map_data = server.adaptor.execute_and_fetchall(sql,)

    ko_list = [i[4] for i in map_data]

    sql = 'select id from comparative_tables_ko where id in (%s);' % ('"' + '","'.join(ko_list) + '"')
    
    ko_list_found_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    # get list of all orthogroups with corresponding KO
    sql = 'select distinct ko_id,orthogroup from enzyme_locus2ko as t1 ' \
            ' where ko_id in (%s);' % ( '"' + '","'.join(ko_list_found_in_db) + '"')
            
    orthogroup_data = server.adaptor.execute_and_fetchall(sql,)
    ko2orthogroups = {}
    orthogroup_list = []
    for i in orthogroup_data:
        if i[0] not in ko2orthogroups:
            ko2orthogroups[i[0]] = [i[1]]
        else:
            ko2orthogroups[i[0]].append(i[1])
        orthogroup_list.append(i[1])

    taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup")
    taxon2ko2count = ete_motifs.get_taxon2name2count(biodb, ko_list_found_in_db, type="ko")

    labels = ko_list_found_in_db
    tree, style = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2ko2count)

    tree2, style2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                    labels,
                                                    taxon2orthogroup2count,
                                                    taxon2ko2count,
                                                    ko2orthogroups)


    if len(labels) > 70:
        big = True
        path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.png' % map_name
        asset_path = '/temp/KEGG_tree_%s.png' % map_name
        tree.render(path, dpi=1200, tree_style=style)

    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.svg' % map_name
        asset_path = '/temp/KEGG_tree_%s.svg' % map_name
        tree.render(path, dpi=800, tree_style=style)

        path2 = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s_complete.svg' % map_name
        asset_path2 = '/temp/KEGG_tree_%s_complete.svg' % map_name

        tree2.render(path2, dpi=800, tree_style=style2)

    sql = 'select bb.ko_accession, count(*) from (select C.ko_accession, E.seqfeature_id from  ' \
            ' (select * from enzyme_kegg_pathway  where pathway_name="%s") A ' \
            ' inner join enzyme_pathway2ko as B  on A.pathway_id=B.pathway_id ' \
            ' inner join enzyme_ko_annotation as C on B.ko_id=C.ko_id ' \
            ' inner join enzyme_seqfeature_id2ko E on B.ko_id=E.ko_id ' \
            ' group by B.ko_id,seqfeature_id) bb group by ko_accession;' % (map_name)

    ko2freq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select ko_accession from enzyme_seqfeature_id2ko t1 ' \
            ' inner join enzyme_pathway2ko t2 on t1.ko_id=t2.ko_id ' \
            ' inner join enzyme_kegg_pathway t3 on t2.pathway_id=t3.pathway_id ' \
            ' inner join orthology_detail t4 on t1.seqfeature_id=t4.seqfeature_id ' \
            ' inner join enzyme_ko_annotation t5 on t1.ko_id=t5.ko_id' \
            ' where taxon_id=%s and pathway_name="%s";' % (taxon_id,
                                                           map_name)

    ko_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    path_map_freq = settings.BASE_DIR + '/assets/temp/KEGG_map_freq_%s' % map_name
    path_map_freq_svg = '/temp/KEGG_map_freq_%s.svg' % map_name

    keg_path_highlight = '+'+'+'.join(ko_list_found_in_db)

    kegg_maps.map2highlighted_map(map_name, ko_list, ko2freq, biodb, path_map_freq+'.pdf', taxon_id=taxon_id)

    sql = 'select t1.description from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
            ' where t2.name="%s" and t1.taxon_id=%s' % (biodb,
                                                        taxon_id)

    organism = server.adaptor.execute_and_fetchall(sql,)[0][0]


    template = Template('''
                    
                    {% load static %}
                    {% load custom_tags %}

                    <h3> {% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}: {% with map_data|first as first_doc %}{{ first_doc.2 }}{% endwith %} ({% with map_data|first as first_doc %}{{ first_doc.1 }}{% endwith %}) {% if organism%}:  {{organism}} {%endif%}</h3>

                    <ul id="tabs" class="nav nav-tabs" data-tabs="tabs">
                        <li class="active"><a href="#tab1" data-toggle="tab">General</a></li>
                        <li><a href="#tab2" data-toggle="tab">Profile</a></li>
                        <li><a href="#tab3" data-toggle="tab">Profile + homologs</a></li>
                        <li><a href="#tab4" data-toggle="tab">KEGG map</a></li>

                    </ul>

                    <div id="my-tab-content" class="tab-content">
                        <div class="tab-pane active" id="tab1">
                            </br>
                            <ol> 
                                <li><a href='http://www.genome.jp/kegg-bin/show_pathway?{% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}{{keg_path_highlight}}'>Link to KEGG website <i class="fa fa-external-link"></i></a></li>
                            </ol>

                            <h4>List of Kegg Orthologs associated to KEGG pathways {% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}</h4>

                            <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                                <div class="panel-heading" style="width:100%">
                                    <h3 class="panel-title">Help</h3>
                                </div>
                                <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                    The "KO coccurence" column report the number of proteins annotated with this Kegg Ortholog in the entire database. 
                                    Click on the KO accession to retrieve the detailed list of proteins.
                                </p>
                            </div>

                            <table id="ko_table">
                                <thead>
                                <tr>
                                    <th>EC(s)</th>
                                    <th>KO</th>
                                    <th>KO occurences</th>
                                    <th>Description</th>
                                </tr>
                                </thead>
                                <tbody>
                                {% for values in map_data%}
                                    <tr>
                                        {%if not ' ' in values.3 and not '-' in values.3 and not ',' in values.3 %}
                                            <td><a href="{% url 'fam'  values.3 'EC' %}" target="_top">{{values.3}}</a></td>
                                        {% else %}
                                            <td>{{values.3}}</td>
                                        {% endif %}
                                        {%if ko2freq|keyvalue:values.4 != None %}
                                            <td><a href="{% url 'fam'  values.4 'ko' %}" target="_top">{{values.4}}</a></td>
                                            <td>{{ko2freq|keyvalue:values.4}}</td>
                                        {%else%}
                                            <td>{{values.4}}</td>
                                            <td>0</td>
                                        {%endif%}

                                        <td>{{values.5}}</td>
                                    </tr>
                                {% endfor %}

                                </tbody>
                            </table>
                        </div>
                        <div class="tab-pane" id="tab2">

                            <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                                    <div class="panel-heading" style="width:100%">
                                        <h3 class="panel-title">Help</h3>
                                    </div>
                                    <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                        Conservation of the pathway among species of the PVC superphylum. A blue background indicate that one or multiple 
                                        proteins were annotated with the corresponding Kegg Ortholog using KoFamScan (<a href="https://www.genome.jp/tools/kofamkoala/">https://www.genome.jp/tools/kofamkoala/</a>).
                                    </p>
                                </div>

                            <a download="profile.svg" class="btn" href="{% static asset_path %}"><i class="fa fa-download"></i> Download SVG</a>
                            <a onclick='exportPNG("profile", "profile");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>

                            <object type="image/svg+xml" data="{% static asset_path %}" id="profile" style="max-width:100%"></object>
                        </div>
                        <div class="tab-pane" id="tab3">

                            <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                                    <div class="panel-heading" style="width:100%">
                                        <h3 class="panel-title">Help</h3>
                                    </div>
                                    <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                        Conservation of the pathway among species of the PVC superphylum. A blue background indicate that one or multiple 
                                        proteins were annotated with the corresponding Kegg Ortholog using KoFamScan (<a href="https://www.genome.jp/tools/kofamkoala/">https://www.genome.jp/tools/kofamkoala/</a>).
                                        A green background indicates that no protein were annoated using KoFamScan but an homolog was identified using <a href="https://github.com/davidemms/OrthoFinder">OrthoFinder</a>. 
                                    </p>
                            </div>

                            <a download="profile_all.svg" class="btn" href="{% static asset_path2 %}"><i class="fa fa-download"></i> Download SVG</a>
                            <a onclick='exportPNG("profile_homologs", "profile_homologs");' class="btn" id="png_button"><i class="fa fa-download"></i> Download PNG</a>

                            <object type="image/svg+xml" data="{% static asset_path2 %}" id="profile_homologs" style="max-width:100%"></object>
                        </div>
                        <div class="tab-pane" id="tab4">

                        <div class="panel panel-success" style="width:80%; top: 400px; margin: 10px 10px 10px 10px">
                            <div class="panel-heading" style="width:100%">
                                <h3 class="panel-title">Help</h3>
                            </div>
                            <p style="margin: 10px 10px 10px 10px; line-height: 180%">
                                {% if organism %}
                                    <strong>Red scale:</strong> no homolog identified in the genome of {{organism}}. The color scale reflects the number of homologs in other genomes (the more the darker). <br>
                                    <strong>Green scale:</strong> homolog identified in the genome of {{organism}}. The color scale reflects the number of homologs in other genomes.
                                    <br>
                                {% else %}
                                    <strong>Green scale:</strong> reflects the conservation of the KO among all organisms included in the database. 
                                    You can click on pathway names and KO accessions to retrieve more information about that step.
                                {% endif %}

                            </p>
                        </div>

                            <object type="image/svg+xml" data="{% static path_map_freq_svg %}"></object>
                        </div>
                    </div>

                <script>
                    $(document).ready(function() {
                        $('#ko_table').DataTable( {
                            dom: 'Bfrtip',
                            "order": [[2, "desc" ]],
                            "pageLength": 15,
                            "paging":   true,
                            "ordering": true,
                            "info":     false,
                            buttons: [
                            {
                                extend: 'excel',
                                title: '{% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}'
                            },
                            {
                                extend: 'csv',
                                title: '{% with map_data|first as first_doc %}{{ first_doc.0 }}{% endwith %}'
                            }
                            ],
                        } );
                    } );
                </script>

        ''')

    html = template.render(Context(locals()))
    
    current_task.update_state(state='SUCCESS',
                              meta={'current': 1,
                                    'total': 1,
                                    'percent': 100,
                                    'description': "Plotting phylogenetic tree"})
    
    return html
