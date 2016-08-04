# -*- coding: utf-8 -*-

# todo circos gc file curently written in home directory, move it to other place
# todo save temp files in temp folder


#from django.shortcuts import render
#from datetime import datetime
from django.shortcuts import render
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import circos_orthology
# from django.core.cache import cache
#import pylibmc
#from django.core.cache import cache
import os
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from forms import make_contact_form
from forms import make_plot_form
from forms import SearchForm
from forms import BiodatabaseForm
from forms import make_circos_form
from forms import make_circos2genomes_form
from forms import make_mummer_form
from forms import make_blast_form
from forms import make_crossplot_form
from forms import ConnexionForm
from forms import DBForm
from forms import make_motif_form
from forms import PCRForm
from forms import make_extract_form
from forms import make_circos_orthology_form
from forms import make_interpro_from
from forms import make_metabo_from
from forms import make_module_overview_form
from forms import make_extract_region_form
from forms import make_venn_from
from forms import make_priam_form
from forms import AnnotForm

from django.contrib.auth import logout
from django.conf import settings
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required

import manipulate_biosqldb
import mysqldb_plot_genomic_feature
from django.core.cache import get_cache
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from gbk2table import Record
import models
import simplejson
import string
import random

def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))


def extract_alphanumeric(input_string):
    from string import ascii_letters, digits
    return "".join([ch for ch in input_string if ch in (ascii_letters + digits + '_-.')])


def choose_db(request):
    server = manipulate_biosqldb.load_db()
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = BiodatabaseForm(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides

            biodb = form.cleaned_data['biodatabase']
            envoi = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = BiodatabaseForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/choose_db.html', locals())


def create_user(username, mail, password, first_name, last_name, staff=True):
    from django.contrib.auth.models import User
    user = User.objects.create_user(username, mail, password)
    user.first_name, user.last_name = first_name, last_name
    user.is_staff = staff
    user.save()



def chlamdb_login(request):
    error = False

    if request.method == "POST":
        print "longin!"
        form = ConnexionForm(request.POST)
        if form.is_valid():
            username = form.cleaned_data["username"]
            password = form.cleaned_data["password"]
            user = authenticate(username=username, password=password)  # Nous vérifions si les données sont correctes
            if user:  # Si l'objet renvoyé n'est pas None
                login(request, user)  # nous connectons l'utilisateur
                print user, "logged"
                #return HttpResponseRedirect("/chlamdb/home")
            else: # sinon une erreur sera affichée
                error = True

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides

            biodb = form.cleaned_data['biodatabase']
            envoi = True
    else:
        form = ConnexionForm()

    return render(request, 'chlamdb/login.html', locals())


def logout_view(request):
    logout(request)
    return render(request, 'chlamdb/logout.html', locals())

@login_required
def home(request, biodb):
    from Bio.SeqUtils import GC

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select * from genomes_info_%s' % biodb

    genomes_data = server.adaptor.execute_and_fetchall(sql,)

    return render(request, 'chlamdb/home.html', locals())


def substription():
    create_user('tpillone', 'trestan.pillonel@gmail.com', 'estrella3', "Trestan", "Pillonel")

#cache = pylibmc.Client(['127.0.0.1:8000'])


@login_required
def circos_homology(request, biodb):


    cache = get_cache('default')
    print "loading db..."
    server, db = manipulate_biosqldb.load_db(biodb)
    print "db loaded..."

    circos_orthology_form_class = make_circos_orthology_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = circos_orthology_form_class(request.POST)

        if form.is_valid():

            accession = form.cleaned_data['accession']

            print "accession", accession

            sql = 'select accession from bioentry' \
                  ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
                  ' and biodatabase.name = "%s"' % biodb



            all_accession = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            columns = 'orthogroup, accession, start, stop'
            sql_ref = 'select %s from orthology_detail_%s where locus_tag = "%s" or protein_id = "%s" or orthogroup = "%s"' % (columns,
                                                                                                          biodb,
                                                                                                          accession,
                                                                                                          accession, accession)




            ref_record = server.adaptor.execute_and_fetchall(sql_ref,)[0]

            orthogroup = ref_record[0]

            columns = 'accession, start, stop'
            sql_targets = 'select %s from orthology_detail_%s where orthogroup ="%s"' % (columns,
                                                                                          biodb,
                                                                                          orthogroup)

            target_records = server.adaptor.execute_and_fetchall(sql_targets,)

            print "ref_record", ref_record
            print "target_records", target_records

            record_list = []
            for accession in all_accession:
                if accession == "CP001848" or accession == "BX119912":
                    continue
                print "accession", accession
                biorecord = cache.get(biodb + "_" + accession)


                if not biorecord:
                    print biodb + "_" + accession, "NOT in memory"
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


            path = settings.BASE_DIR + "/assets/circos"

            circos_orthology.circos_orthology(record_list, ref_record[1:], target_records, location = path)
            circos_file = "circos/circos_ortho.svg"
            envoi_circos = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = circos_orthology_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/circos_homology.html', locals())

@login_required
def extract_orthogroup(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        #form2 = ContactForm(request.POST)
        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables

            print request.POST
            print form.cleaned_data.keys()
            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]



            try:
                single_copy = request.POST['button_single_copy']
                single_copy = True
            except:
                single_copy = False
            n_missing = form.cleaned_data['frequency']

            if int(n_missing)>=len(include):
                wrong_n_missing = True
            else:
                server, db = manipulate_biosqldb.load_db(biodb)


                freq_missing = (len(include)-float(n_missing))/len(include)

                # get sub matrix and complete matrix
                mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                              "orthology",
                                                                              "orthogroup",
                                                                              include,
                                                                              exclude,
                                                                              freq_missing,
                                                                              single_copy=single_copy)

                match_groups = mat.index.tolist()

                if len(match_groups) == 0:
                    no_match = True
                else:

                    # get count in subgroup
                    orthogroup2count = dict((mat > 0).sum(axis=1))
                    # get count in complete database
                    orthogroup2count_all = dict((mat_all > 0).sum(axis=1))

                    #print cog2count_all
                    max_n = max(orthogroup2count_all.values())

                    # GET max frequency for template
                    sum_group = len(match_groups)

                    match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(match_groups, biodb)

                    envoi_extract = True

                    circos_url = '?ref=%s&' % reference_taxon
                    circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(match_groups)
                    print "circos_url", circos_url

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_orthogroup.html', locals())



@login_required
def orthogroup_annotation(request, biodb):
    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = AnnotForm(request.POST)  # Nous reprenons les données

        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables

            server, db = manipulate_biosqldb.load_db(biodb)

            match_groups = [i.rstrip() for i in form.cleaned_data['orthogroups'].rstrip().split('\n')]
            print 'match groups', match_groups

            match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(match_groups, biodb)





            envoi_annot = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form = AnnotForm()
    return render(request, 'chlamdb/orthogroup_annotation.html', locals())


@login_required
def venn_orthogroup(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_orthogroups_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select orthogroup from orthology_%s where `%s` > 0' % (biodb, target)
                print sql
                orthogroups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_orthogroups_list += orthogroups
                data = '"' + '","'.join(orthogroups) + '"'
                series+=template_serie % (taxon_id2genome[target], orthogroups) + ','
            series = series[0:-1] + ']'


            #h['Marilyn Monroe'] = 1;

            orthogroup2description = ''
            sql = 'select orthogroup, gene from orthology_detail_%s' % biodb
            data = server.adaptor.execute_and_fetchall(sql,)
            orthogroup2genes = {}
            for i in data:
                if i[0] not in orthogroup2genes:
                    orthogroup2genes[i[0]] = [i[1]]
                else:
                    if i[1] in orthogroup2genes[i[0]]:
                        pass
                    else:
                        orthogroup2genes[i[0]].append(i[1])
            sql = 'select orthogroup, product from orthology_detail_%s' % biodb
            data = server.adaptor.execute_and_fetchall(sql,)
            orthogroup2product = {}
            for i in data:
                if i[0] not in orthogroup2product:
                    orthogroup2product[i[0]] = [i[1]]
                else:
                    if i[1] in orthogroup2product[i[0]]:
                        pass
                    else:
                        orthogroup2product[i[0]].append(i[1])

            for i in orthogroup2genes:
                if i in all_orthogroups_list:
                    genes = '<br>'.join(orthogroup2genes[i])
                    products = '<br>'.join(orthogroup2product[i])
                    orthogroup2description+='h["%s"] = "%s</td><td>%s;"\n' % (i, genes, products)
                else:
                    continue
            print orthogroup2description
            #print series
            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_orthogroup.html', locals())


@login_required
def extract_pfam(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables
            print request.POST
            print form.cleaned_data.keys()
            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            n_missing = form.cleaned_data['frequency']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]


            if int(n_missing)>=len(include):
                wrong_n_missing = True
            else:
                server, db = manipulate_biosqldb.load_db(biodb)


                freq_missing = (len(include)-float(n_missing))/len(include)

                # get sub matrix and complete matrix
                mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                              "Pfam",
                                                                              "id",
                                                                              include,
                                                                              exclude,
                                                                              freq_missing)

                match_groups = mat.index.tolist()

                if len(match_groups) == 0:
                    no_match = True
                else:

                    # get count in subgroup
                    pfam2count = dict((mat > 0).sum(axis=1))
                    # get count in complete database
                    pfam2count_all = dict((mat_all > 0).sum(axis=1))

                    #print cog2count_all
                    max_n = max(pfam2count_all.values())

                    # GET max frequency for template
                    sum_group = len(match_groups)

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
                    group_filter += ')'

                    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                              'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
                    sql_2 = 'select %s from orthology_detail_%s %s' % (columns, biodb, group_filter)
                    #print sql_2
                    raw_data = server.adaptor.execute_and_fetchall(sql_2,)


                    import biosql_own_sql_tables
                    pfam2descr = biosql_own_sql_tables.pfam2description(biodb)
                    match_groups_data = []
                    for i, pfam in enumerate(match_groups):
                        match_groups_data.append([i, pfam, pfam2descr[pfam], pfam2count[pfam], pfam2count_all[pfam]])

                    n = 1
                    extract_result = []
                    for one_hit in raw_data:
                        extract_result.append((n,) + one_hit)
                        n+=1

                    envoi_extract = True
                    asset_path = '/assets/temp/profil_tree.svg'

                    print 'getting locus list'
                    motif_list = '"' + '","'.join(match_groups) + '"'

                    locus_list_sql = 'select locus_tag from interpro_%s where taxon_id=%s ' \
                                 ' and signature_accession in (%s)' % (biodb, reference_taxon, motif_list)
                    print locus_list_sql
                    locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
                    print locus_list
                    circos_url = '?ref=%s&' % reference_taxon
                    circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)
                    print "circos_url", circos_url



    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_Pfam.html', locals())


@login_required
def extract_ko(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables

            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            n_missing = form.cleaned_data['frequency']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]

            if int(n_missing)>=len(include):
                wrong_n_missing = True
            else:
                server, db = manipulate_biosqldb.load_db(biodb)


                freq_missing = (len(include)-float(n_missing))/len(include)

                # get sub matrix and complete matrix
                mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                              "ko",
                                                                              "id",
                                                                              include,
                                                                              exclude,
                                                                              freq_missing)

                match_groups = mat.index.tolist()

                if len(match_groups) == 0:
                    no_match = True
                else:

                    # get count in subgroup
                    ko2count = dict((mat > 0).sum(axis=1))
                    # get count in complete database
                    ko2count_all = dict((mat_all > 0).sum(axis=1))

                    #print cog2count_all
                    max_n = max(ko2count_all.values())

                    # GET max frequency for template
                    sum_group = len(match_groups)


            sql = 'select ec, value, pathway_name, pathway_category, description from ' \
            ' (select enzyme_id, ec,value from enzyme.enzymes as t1 inner join enzyme.enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id ' \
            ' where line="description") A left join enzyme.kegg2ec as B on A.enzyme_id=B.ec_id ' \
            ' left join enzyme.kegg_pathway on B.pathway_id=kegg_pathway.pathway_id;'

            sql = 'select ko_id, name, definition, EC, pathways, modules from enzyme.ko_annotation;'
            print sql
            ko2description_raw = server.adaptor.execute_and_fetchall(sql,)

            ko2description_dico = {}


            sql = 'select module_name,description from enzyme.kegg_module'
            sql2 = 'select pathway_name,description from enzyme.kegg_pathway'
            module2category = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            map2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

            for i in ko2description_raw:
                #if i[3] != "1.0 Global and overview maps":

                ko2description_dico[i[0]] = [list(i[1:4])]
                if i[4] != '-':
                    path_str = ''
                    for path in i[4].split(','):
                        try:
                            path_str+='<a href="http://www.genome.jp/dbget-bin/www_bget?map%s">%s</a><br>' % (path[2:], map2description['map'+path[2:]])
                        except:
                            path_str+='<a href="http://www.genome.jp/dbget-bin/www_bget?map%s">%s</a><br>' % (path[2:], path)
                    ko2description_dico[i[0]][0].append(path_str[0:-1])
                else:
                    ko2description_dico[i[0]][0].append('-')
                if i[5] != '-':
                    mod_str = ''
                    for mod in i[5].split(','):
                        mod_str+='<a href="http://www.genome.jp/dbget-bin/www_bget?md:%s">%s</a><br>' % (mod, module2category[mod])
                    ko2description_dico[i[0]][0].append(mod_str[0:-5])
                else:
                    ko2description_dico[i[0]][0].append('-')

                #print ko2description_dico

            #print 'ko2description_dico[', ko2description_dico
            #enzyme2data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            match_groups_data = []

            for i, ko in enumerate(match_groups):
                print i, ko
                for one_pathway in ko2description_dico[ko]:
                    match_groups_data.append([i, ko, one_pathway, ko2count[ko], ko2count_all[ko]])


            ko_list = '"' + '","'.join(match_groups) + '"'

            #print extract_result
            locus_list_sql = 'select locus_tag from enzyme.locus2ko_%s where taxon_id=%s and ko_id in (%s);' % (biodb,
                                                                                                                 reference_taxon,
                                                                                                                 ko_list)
            print locus_list_sql
            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
            print locus_list
            circos_url = '?ref=%s&' % reference_taxon
            circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)
            print "circos_url", circos_url

            envoi_extract = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_ko.html', locals())


@login_required
def extract_EC(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables

            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            n_missing = form.cleaned_data['frequency']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]

            if int(n_missing)>=len(include):
                wrong_n_missing = True
            else:
                server, db = manipulate_biosqldb.load_db(biodb)


                freq_missing = (len(include)-float(n_missing))/len(include)

                # get sub matrix and complete matrix
                mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                              "EC",
                                                                              "id",
                                                                              include,
                                                                              exclude,
                                                                              freq_missing)

                match_groups = mat.index.tolist()

                if len(match_groups) == 0:
                    no_match = True
                else:

                    # get count in subgroup
                    ec2count = dict((mat > 0).sum(axis=1))
                    # get count in complete database
                    ec2count_all = dict((mat_all > 0).sum(axis=1))

                    #print cog2count_all
                    max_n = max(ec2count_all.values())

                    # GET max frequency for template
                    sum_group = len(match_groups)


            sql = 'select ec, value, pathway_name, pathway_category, description from ' \
            ' (select enzyme_id, ec,value from enzyme.enzymes as t1 inner join enzyme.enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id ' \
            ' where line="description") A left join enzyme.kegg2ec as B on A.enzyme_id=B.ec_id ' \
            ' left join enzyme.kegg_pathway on B.pathway_id=kegg_pathway.pathway_id;'
            ec2description_raw = server.adaptor.execute_and_fetchall(sql,)
            print sql
            ec2description_dico = {}

            for i in ec2description_raw:
                #if i[3] != "1.0 Global and overview maps":
                if i[0] not in ec2description_dico:
                    ec2description_dico[i[0]] = [list(i[1:len(i)])]
                else:
                    ec2description_dico[i[0]].append(list(i[1:len(i)]))

            #enzyme2data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            match_groups_data = []

            for i, ec in enumerate(match_groups):
                for one_pathway in ec2description_dico[ec]:
                    match_groups_data.append([i, ec, one_pathway, ec2count[ec], ec2count_all[ec]])


            EC_list = '"' + '","'.join(match_groups) + '"'

            biodb_id = server.adaptor.execute_and_fetchall('select biodatabase_id from biodatabase where name="%s"' % biodb,)[0][0]

            #print extract_result
            locus_list_sql = 'select locus_tag from (select taxon_id,locus_tag,ec_id from enzyme.locus2ec_%s as t1  ' \
                             ' inner join biosqldb.bioentry as t2 on t1.accession=t2.accession ' \
                             ' where biodatabase_id=%s) A inner join enzyme.enzymes as B on A.ec_id=B.enzyme_id' \
                             ' where A.taxon_id=%s and B.ec in (%s);' % (biodb,
                                                                         biodb_id,
                                                                         reference_taxon,
                                                                         EC_list)
            print locus_list_sql
            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
            print locus_list
            circos_url = '?ref=%s&' % reference_taxon
            circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)
            print "circos_url", circos_url



            # get phylogenetic profile of match if not too big
            if len(match_groups) < 50:
                import ete_motifs
                sql = 'select distinct ec,orthogroup from enzyme.locus2ec_%s as t1 ' \
                      ' inner join enzyme.enzymes as t2 on t1.ec_id=t2.enzyme_id where ec in (%s);' % (biodb,
                                                                  '"' + '","'.join(match_groups) + '"')
                orthogroup_data = server.adaptor.execute_and_fetchall(sql,)
                ec2orthogroups = {}
                orthogroup_list = []
                for i in orthogroup_data:
                    if i[0] not in ec2orthogroups:
                        ec2orthogroups[i[0]] = [i[1]]
                    else:
                        ec2orthogroups[i[0]].append(i[1])
                    orthogroup_list.append(i[1])

                taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup")
                taxon2enzyme2count = ete_motifs.get_taxon2name2count(biodb, match_groups, type="EC")


                labels = match_groups
                tree2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                             labels,
                                                             taxon2orthogroup2count,
                                                             taxon2enzyme2count,
                                                             ec2orthogroups)



                if len(labels) > 40:
                    print 'BIGGGGGGGGGGG', len(labels)
                    big = True
                    path = settings.BASE_DIR + '/assets/temp/profil_tree.png'
                    asset_path = '/assets/temp/profil_tree.png'
                    tree2.render(path, dpi=1200, h=600)



                else:
                    print 'not BIGGGGGGGGGG', len(labels)
                    big = False

                    path2 = settings.BASE_DIR + '/assets/temp/profil_tree.svg'
                    asset_path = '/assets/temp/profil_tree.svg'

                    tree2.render(path2, dpi=800, h=600)





            envoi_extract = True








    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_EC.html', locals())

@login_required
def venn_pfam(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_pfam_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.Pfam_%s where `%s` > 0' % (biodb, target)
                print sql
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_pfam_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            pfam2description = ''
            sql = 'select signature_accession, signature_description, count(*) from interpro_%s where analysis="Pfam" group by signature_accession;' % biodb
            print sql
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            for i in data:
                if i in all_pfam_list:
                    #print 'ok'
                    pfam2description+='h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1])
                else:
                    pass
                    #print 'pas ok'
            #print pfam2description
            #print series

            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_Pfam.html', locals())


@login_required
def venn_EC(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_ec_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.EC_%s where `%s` > 0' % (biodb, target)
                print sql
                ec_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_ec_list += ec_list
                data = '"' + '","'.join(ec_list) + '"'
                series+=template_serie % (taxon_id2genome[target], ec_list) + ','
            series = series[0:-1] + ']'


            ec2description = ''
            sql = 'select ec, value, pathway_name, pathway_category, description from ' \
                  ' (select enzyme_id, ec,value from enzyme.enzymes as t1 inner join enzyme.enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id ' \
                  ' where line="description") A left join enzyme.kegg2ec as B on A.enzyme_id=B.ec_id ' \
                  ' left join enzyme.kegg_pathway on B.pathway_id=kegg_pathway.pathway_id;'
            ec2description_raw = server.adaptor.execute_and_fetchall(sql,)
            print sql
            ec2description_dico = {}

            for i in ec2description_raw:
                #if i[3] != "1.0 Global and overview maps":
                if i[0] not in ec2description_dico:
                    ec2description_dico[i[0]] = [list(i[1:len(i)])]
                else:
                    ec2description_dico[i[0]].append(list(i[1:len(i)]))

            for one_ec in all_ec_list:
                data = ec2description_dico[one_ec]
                tmp_str = ''
                for one_pathway in data:
                    tmp_str+= "<td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr><tr>" % (one_ec,
                                                                                                    one_pathway[0],
                                                                                                    one_pathway[1],
                                                                                                    one_pathway[2],
                                                                                                    one_pathway[3])
                ec2description+='h["%s"] = "%s";' % (one_ec, tmp_str[0:-12])


            print ec2description
            #print series

            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_EC.html', locals())


@login_required
def extract_interpro(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données

        #form2 = ContactForm(request.POST)
        if 'comparison' in request.POST and form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables

            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            n_missing = form.cleaned_data['frequency']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]

            if int(n_missing)>=len(include):
                wrong_n_missing = True
            else:
                freq_missing = (len(include)-float(n_missing))/len(include)
                print 'freq_missing', freq_missing

                # get sub matrix and complete matrix
                mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                              "interpro",
                                                                              "id",
                                                                              include,
                                                                              exclude,
                                                                              freq_missing)

                match_groups = mat.index.tolist()
                # get count in subgroup
                interpro2count = dict((mat > 0).sum(axis=1))
                # get count in complete database
                interpro2count_all = dict((mat_all > 0).sum(axis=1))
                print interpro2count_all

                #print interpro2count_all.values()

                #print cog2count_all
                max_n = max(list(interpro2count_all.values()))

                # GET max frequency for template
                sum_group = len(match_groups)

                filter = '"' + '","'.join(match_groups) + '"'

                sql2 = 'select interpro_accession, interpro_description from interpro_%s' \
                ' where interpro_accession in (%s) group by interpro_accession;' % (biodb, filter)

                raw_data = list(server.adaptor.execute_and_fetchall(sql2,))

                match_data = []
                for one_match in raw_data:
                    match_data.append(list(one_match)+[interpro2count[one_match[0]], interpro2count_all[one_match[0]]])

                print match_data[0]
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
                #print group_filter


                columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                          'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
                sql_2 = 'select %s from orthology_detail_%s %s' % (columns, biodb, group_filter)
                #print sql_2
                raw_data = server.adaptor.execute_and_fetchall(sql_2,)

                n = 1
                extract_result = []
                for one_hit in raw_data:
                    extract_result.append((n,) + one_hit)
                    n+=1
                    #print n

                interpro_list = '"' + '","'.join(match_groups) + '"'

                #print extract_result
                locus_list_sql = 'select locus_tag from interpro_%s where taxon_id=%s ' \
                             ' and interpro_accession in (%s)' % (biodb, reference_taxon, interpro_list)
                print locus_list_sql
                locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
                print locus_list
                circos_url = '?ref=%s&' % reference_taxon
                circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)
                print "circos_url", circos_url
                envoi_extract = True


    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_interpro.html', locals())




@login_required
def venn_interpro(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if 'venn' in request.POST and form_venn.is_valid():
            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_pfam_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.interpro_%s where `%s` > 0' % (biodb, target)
                print sql
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_pfam_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            interpro2description = ''
            sql = 'select interpro_accession,interpro_description, count(*) from interpro_%s group by interpro_accession;' % biodb
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for i in data:
                if i in all_pfam_list:
                    print 'ok'
                    interpro2description+='h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1])
                else:
                    print 'pas ok'
            print interpro2description
            #print series

            envoi_venn = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET  # Nous créons un formulaire vide
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_interpro.html', locals())


@login_required
def extract_cog(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    extract_form_class = make_extract_form(biodb)


    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = extract_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            import biosql_own_sql_tables

            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            n_missing = form.cleaned_data['frequency']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]

            try:
                single_copy = request.POST['button_single_copy']
                single_copy = True
            except:
                single_copy = False

            if int(n_missing)>=len(include):
                wrong_n_missing = True
            else:
                server, db = manipulate_biosqldb.load_db(biodb)


                freq_missing = (len(include)-float(n_missing))/len(include)

                # get sub matrix and complete matrix
                mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                              "COG",
                                                                              "id",
                                                                              include,
                                                                              exclude,
                                                                              freq_missing)

                match_groups = mat.index.tolist()

                # get count in subgroup
                cog2count = dict((mat > 0).sum(axis=1))
                # get count in complete database
                cog2count_all = dict((mat_all > 0).sum(axis=1))

                #print cog2count_all
                max_n = max(cog2count_all.values())

                # GET max frequency for template
                sum_group = len(match_groups)

                # get data for each matching cog
                cog_data = []
                for i in match_groups:
                    sql = 'select * from COG.cog_names_2014 where COG_id ="%s"' % i
                    tmp = list(server.adaptor.execute_and_fetchall(sql,)[0])
                    cog_data.append(tmp+[cog2count[tmp[0]], cog2count_all[tmp[0]]])
                print cog_data
                interpro_list = '"' + '","'.join(match_groups) + '"'

                biodb_id = server.adaptor.execute_and_fetchall('select biodatabase_id from biodatabase where name="%s"' % biodb,)[0][0]

                #print extract_result
                locus_list_sql = 'select locus_tag from (select taxon_id,locus_tag,COG_id from COG.locus_tag2gi_hit_%s as t1 ' \
                                 ' inner join biosqldb.bioentry as t2 on t1.accession=t2.accession ' \
                                 ' where biodatabase_id=%s) A where A.taxon_id=%s and A.COG_id in (%s);' % (biodb,
                                                                                                            biodb_id,
                                                                                                            reference_taxon,
                                                                                                            interpro_list)
                print locus_list_sql
                locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
                print locus_list
                circos_url = '?ref=%s&' % reference_taxon
                circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)
                print "circos_url", circos_url

                envoi_extract = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = extract_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/extract_cogs.html', locals())



@login_required
def venn_cog(request, biodb):

    cache = get_cache('default')
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    venn_form_class = make_venn_from(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form_venn = venn_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form_venn.is_valid():  # Nous vérifions que les données envoyées sont valides

            targets = form_venn.cleaned_data['targets']

            server, db = manipulate_biosqldb.load_db(biodb)

            all_cog_list = []
            series = '['
            taxon_id2genome = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            for target in targets:
                template_serie = '{name: "%s", data: %s}'
                sql ='select id from comparative_tables.COG_%s where `%s` > 0' % (biodb, target)
                print sql
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_cog_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            #h['Marilyn Monroe'] = 1;

            cog2description = []
            sql = 'select * from COG.cog_names_2014'
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for i in data:
                if i in all_cog_list:

                    #print 'ok'
                    cog2description.append('h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1]))
                else:
                    pass
                    #print 'pas ok'
            #cog2description = cog2description[0:120000]
            print 'COG0755' in all_cog_list
            print len(cog2description)

            #print cog2description
            #print 'all_cog_list', all_cog_list
            #print series
            envoi_venn = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form_venn = venn_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/venn_cogs.html', locals())



@login_required
def extract_region(request, biodb):
    extract_region_form_class = make_extract_region_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)
    from Bio.Alphabet import generic_protein
    import re

    cache = get_cache('default')

    if request.method == 'POST':

        form = extract_region_form_class(request.POST)

        if form.is_valid():
            genome_accession = form.cleaned_data['genome']

            region = form.cleaned_data['region']

            start_stop = re.sub(' ', '', form.cleaned_data['region']).split(",")

            extract = form.cleaned_data['extract']

            genome_description = manipulate_biosqldb.accession2description(server, biodb)[genome_accession]


            if extract == 'annotation':
                get_annotation = True
                columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                          'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'


                sql2 = 'select %s from orthology_detail_%s where start > %s and stop < %s and accession="%s"' % (columns,
                                                                                                                 biodb,
                                                                                                                 start_stop[0],
                                                                                                                 start_stop[1],
                                                                                                                 genome_accession)


                try:
                    raw_data = server.adaptor.execute_and_fetchall(sql2, )
                except IndexError:
                    valid_id = False
                    return render(request, 'chlamdb/extract_region.html', locals())
                if not raw_data:
                        valid_id = False
                else:
                    n = 1
                    search_result = []
                    for one_hit in raw_data:

                        sql = 'select contig from '


                        if one_hit[2] != '-':
                            interpro_id = one_hit[2]
                        else:
                            interpro_id = one_hit[1]
                        search_result.append((n,) + one_hit + (interpro_id,))
                        n+=1


            if extract == 'sequence' or extract == 'sequence_trans':
                get_sequence = True
                seq = manipulate_biosqldb.location2sequence(server, genome_accession, biodb, start_stop[0], int(start_stop[1])-int(start_stop[0]))

                record = SeqRecord(Seq(seq, generic_protein),
                   id="%s_%s_%s" % (genome_accession, start_stop[0], start_stop[1]), description=genome_description)

                temp_file = os.path.join(settings.BASE_DIR, "assets/temp/%s_region.fa" % genome_accession)
                temp_location = "temp/%s_region.fa" % genome_accession
                print "temp_file", temp_file
                with open(temp_file, 'w') as f:
                    SeqIO.write(record, f, 'fasta')


                if extract == 'sequence_trans':
                    extract_trans = True
                    from Bio.Seq import reverse_complement, translate
                    anti = reverse_complement(seq)
                    comp = anti[::-1]
                    length = len(seq)
                    frames = {}
                    for i in range(0, 3):
                        fragment_length = 3 * ((length-i) // 3)
                        frames[i+1] = translate(seq[i:i+fragment_length], 1)
                        frames[-(i+1)] = translate(anti[i:i+fragment_length], 1)[::-1]



                    frame_plus_1 = frames[1]
                    frame_plus_2 = frames[2]
                    frame_plus_3 = frames[3]
                    frame_minus_1 = frames[-1]
                    frame_minus_2 = frames[-2]
                    frame_minus_3 = frames[-3]



            envoi = True

    else:
        form = extract_region_form_class()
    return render(request, 'chlamdb/extract_region.html', locals())


@login_required
def locusx(request, biodb, locus=None, menu=False):

    print 'biodb', biodb
    print 'locus', locus
    print 'menu', menu

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST

        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        #sql1 = 'SELECT column_name FROM information_schema.columns WHERE table_name="orthology_detail_chlamydia_03_15"'

        sql1 =   'SELECT' \
                 ' CASE' \
                 '   WHEN locus_tag = "%s" THEN "locus_tag"' \
                 '   WHEN protein_id = "%s" THEN "protein_id"' \
                 '   WHEN orthogroup = "%s" THEN "orthogroup"'\
                 ' END AS "which_column"'\
                 ' FROM' \
                 ' orthology_detail_%s where locus_tag="%s" or protein_id like "%%%%%s%%%%" or orthogroup="%s"' % (locus,
                                                                                                                   locus,
                                                                                                                   locus,
                                                                                                                   biodb,
                                                                                                                   locus,
                                                                                                                   locus,
                                                                                                                   locus)



        try:
            print sql1
            input_type = server.adaptor.execute_and_fetchall(sql1, )[0][0]
        except IndexError:
            print 'not a valid id'
            valid_id = False
            return render(request, 'chlamdb/locus.html', locals())

        else:
            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql2 = 'select %s from orthology_detail_%s where %s="%s"' % (columns, biodb, input_type, locus)
            if input_type == 'locus_tag':
                sql4 = 'select accession from orthology_detail_%s where locus_tag="%s" limit 1' % (biodb, locus)
                genome_accession = server.adaptor.execute_and_fetchall(sql4,)[0][0]

                sql3 = 'select t2.COG_id,t2.functon,t2.name from COG.locus_tag2gi_hit_%s ' \
                       ' as t1 inner join COG.cog_names_2014 as t2 on t1.COG_id=t2.COG_id where locus_tag="%s"' % (biodb, locus)
                sql4 = 'select analysis, signature_accession, signature_description, interpro_accession, interpro_description ' \
                       ' from interpro_%s where locus_tag="%s";' % (biodb, locus)
                try:
                    cog_data = server.adaptor.execute_and_fetchall(sql3, )[0]

                except IndexError:
                    cog_data = False
                try:
                    interpro_data = server.adaptor.execute_and_fetchall(sql4, )
                except IndexError:
                    interpro_data= False

                try:
                    sql_interpro = 'select interpro_accession, interpro_description from interpro_%s' \
                                   ' where locus_tag="%s" and interpro_accession !="0"' \
                                   ' group by interpro_accession;' % (biodb, locus)
                    interpro_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_interpro, )]

                except:
                    interpro_data = False

                try:
                    sql_pfam = 'select signature_accession, signature_description' \
                               ' from interpro_%s where locus_tag="%s" ' \
                               ' and analysis="Pfam";' % (biodb, locus)
                    pfam_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pfam, )]
                    print pfam_data

                except:
                    pfam_data = False

                try:
                    sql_pathway = 'select locus_tag, pathways, interpro_description from interpro_%s where ' \
                                  ' locus_tag="%s" and pathways!="0" group by pathways;'  % (biodb, locus)

                    pathway_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pathway, )]
                    all_path = {}
                    for one_path in pathway_data:
                        data = one_path[1].split('|')
                        for one_db in data:
                            one_db_info = one_db.split(':')
                            db = one_db_info[0]
                            if db not in all_path:
                                all_path[db] = [one_db_info[1][1:]]
                            else:
                                if one_db_info[1][1:] not in all_path[db]:
                                    all_path[db].append(one_db_info[1][1:])
                except:
                    pathways_data = False



            data = list(server.adaptor.execute_and_fetchall(sql2, )[0])

            if input_type == 'locus_tag':
                seq_start = int(data[3])
                seq_end = int(data[4])
                strand = int(data[5])
                leng = (seq_end-seq_start)+100

                seq = manipulate_biosqldb.location2sequence(server, genome_accession, biodb, seq_start-50, leng)
                if strand == -1:

                    from Bio.Seq import Seq
                    seq_obj = Seq(seq)
                    print seq
                    seq = str(seq_obj.reverse_complement())
                    print seq

                    seq = seq[0:49] + '<font color="red">' + seq[49:-50] + '</font>' + seq[-50:len(seq)]

                else:
                    seq = seq[0:50] + '<font color="red">' + seq[50:-50] + '</font>' + seq[-50:len(seq)]

            if data[2] == '-':
                data[2] = data[1]

            orthogroup = data[0]

            fasta = "%s_fasta/%s.txt" % (biodb, orthogroup)
            alignment = "%s_fasta/%s.html" % (biodb, orthogroup)
            alignment_fasta = "%s_fasta/%s.fa" % (biodb, orthogroup)
            alignment_fasta_nucl = "%s_fasta_nucl/%s_nucl.txt" % (biodb, orthogroup)
            tree_unrooted = "%s_fasta/%s_tree.svg" % (biodb, orthogroup)
            tree_rooted = "%s_fasta/%s_tree_reroot.svg" % (biodb, orthogroup)
            tree_file = "%s_fasta/%s.phy_phyml_tree.txt" % (biodb, orthogroup)

            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                  'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql3 = 'select %s from orthology_detail_%s where orthogroup = "%s" ' % (columns, biodb, orthogroup)

            homologues = list(server.adaptor.execute_and_fetchall(sql3, ))
            print homologues



            if len(homologues) >1:
                orthologs = True
            else:
                orthologs = False
            import orthogroup_identity_db
            if len(homologues) > 1:
                orthogroup2identity_dico = orthogroup_identity_db.orthogroup2identity_dico(biodb, orthogroup)

                print "orthologs", orthologs, len(homologues)
                for count, value in enumerate(homologues):
                    value = list(value)
                    locus_2 = value[1]
                    if value[2] != '-':
                        interpro_id = value[2]
                    else:
                        value[2] = value[1]
                    #print value + [orthogroup2identity_dico[data[1]][locus_2]]
                    homologues[count] = [count+1] + value + [orthogroup2identity_dico[data[1]][locus_2]]
                    #print homologues[count]

            else:
                homologues[0] = (1,) + homologues[0] + (100,)

        home_dir = os.path.dirname(os.path.realpath(__file__))
        print 'home', home_dir
        local_file = "/../assets/%s/interpro/%s.html" % (biodb, data[2])
        interpro_check = home_dir + local_file
        print "interpro_check", interpro_check
        if os.path.isfile(interpro_check):
            interpro_protein = True
        else:
            interpro_protein = False

        envoi = True


    return render(request, 'chlamdb/locus.html', locals())


@login_required
def fam(request, biodb, fam, type):

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST

        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        print 'type', type, fam

        #sql1 = 'SELECT column_name FROM information_schema.columns WHERE table_name="orthology_detail_chlamydia_03_15"'
        if type =='pfam':
            sql1 =   'select locus_tag from interpro_%s where signature_accession="%s" group by locus_tag' % (biodb, fam)
            sql2 = 'select signature_description from interpro_%s where signature_accession="%s" limit 1' % (biodb, fam)
            info = server.adaptor.execute_and_fetchall(sql2, )[0]
        elif type == 'cog':
            sql1 = 'select locus_tag from COG.locus_tag2gi_hit_%s where COG_id="%s"' % (biodb, fam)
            sql2 = 'select functon, name from COG.cog_names_2014 where COG_id = "%s"' % (fam)
            print 'cog type', sql2
            info = server.adaptor.execute_and_fetchall(sql2, )[0]

        elif type == 'interpro':
            sql1 = 'select locus_tag from interpro_%s where interpro_accession="%s" group by locus_tag' % (biodb, fam)
            sql2 = 'select signature_description from interpro_%s where interpro_accession="%s" limit 1' % (biodb, fam)
            print sql1
            print sql2
            info = server.adaptor.execute_and_fetchall(sql2, )[0]

        elif type == 'EC':
            sql1 = 'select locus_tag from enzyme.locus2ec_%s as t1 ' \
                   ' inner join enzyme.enzymes as t2 on t1.ec_id=t2.enzyme_id where ec="%s" group by locus_tag;' % (biodb,
                                                                                                                    fam)
            sql2 = 'select line,value from (select * from enzyme.enzymes where ec="%s") t1 ' \
                   ' inner join enzyme.enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id;' % (fam)
            path = fam.split('.')
            print path
            external_link = 'http://www.chem.qmul.ac.uk/iubmb/enzyme/EC%s/%s/%s/%s.html' % (path[0], path[1], path[2], path[3])
            print sql1
            print sql2

            sql_pathways = 'select pathway_name,pathway_category,description ' \
                           ' from (select * from enzyme.enzymes where ec = "%s") t1 ' \
                           ' inner join enzyme.kegg2ec as t2 on t2.ec_id=t1.enzyme_id ' \
                           ' inner join enzyme.kegg_pathway as t3 on t2.pathway_id=t3.pathway_id' \
                           ' where pathway_category !="1.0 Global and overview maps";' % (fam)

            pathway_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pathways, )]
            print "pathway_data",pathway_data
            info =  server.adaptor.execute_and_fetchall(sql2, )

        elif type == 'ko':
            sql1 = 'select locus_tag from enzyme.locus2ko_%s where ko_id="%s" group by locus_tag;' % (biodb,
                                                                                                      fam)
            sql2 = 'selectselect * from enzyme.ko_annotation where ko_id="%s"' % (fam)

            external_link = 'http://www.genome.jp/dbget-bin/www_bget?%s' % (fam)
            print sql1
            print sql2

            sql_modules = 'select pathways, modules from enzyme.ko_annotation where ko_id="%s";' % (fam)
            data = server.adaptor.execute_and_fetchall(sql_modules,)[0]
            if len(data)>0:
                print data
                pathway_list = ','.split(data[0])
                module_list = ','.split(data[0])

        else:
            valid_id = False
            return render(request, 'chlamdb/fam.html', locals())
        try:
            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1, )]
            print 'locus', locus_list
            locus_list_form = '"' + '","'.join(locus_list) + '"'
        except IndexError:
            valid_id = False
            return render(request, 'chlamdb/fam.html', locals())
        else:
            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql2 = 'select %s from orthology_detail_%s where locus_tag in (%s)' % (columns, biodb, locus_list_form)

            all_locus_raw_data = server.adaptor.execute_and_fetchall(sql2, )
            orthogroup_list = [i[0] for i in all_locus_raw_data]
            all_locus_data = []
            group_count = []
            for i in range(0, len(all_locus_raw_data)):
                 all_locus_data.append([i] + list(all_locus_raw_data[i]))
                 if all_locus_raw_data[i][0] not in group_count:
                    group_count.append(all_locus_raw_data[i][0])
        envoi = True
        menu = True
        import ete_motifs
        if type =='pfam':
            taxon2orthogroup2count_reference = get_taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb, [fam], 'Pfam')

            sql3= 'select distinct taxon_id,orthogroup,signature_accession from interpro_%s ' \
                  ' where analysis="Pfam" and orthogroup in (%s);' % (biodb,'"'+'","'.join(set(orthogroup_list))+'"')




            print "taxon2orthogroup2count_reference", taxon2orthogroup2count_reference
        elif type == 'cog':
            taxon2orthogroup2count_reference = get_taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb, [fam], 'COG')

            sql3='select distinct taxon_id,orthogroup,COG_id from (select taxon_id,locus_tag,orthogroup ' \
                 ' from biosqldb.orthology_detail_%s where orthogroup in (%s)) A ' \
                 ' inner join COG.locus_tag2gi_hit_%s as B on A.locus_tag=B.locus_tag;' % (biodb,'"'+'","'.join(set(orthogroup_list))+'"', biodb)

            print "taxon2orthogroup2count_reference", taxon2orthogroup2count_reference
        elif type == 'interpro':
            taxon2orthogroup2count_reference = get_taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb, [fam], 'interpro')

            sql3 = 'select distinct taxon_id,orthogroup,interpro_accession from ' \
                   ' interpro_%s where orthogroup in (%s);' % (biodb,'"'+'","'.join(set(orthogroup_list))+'"')

        elif type == 'EC':
            taxon2orthogroup2count_reference = get_taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb, [fam], 'EC')
            sql3 = 'select distinct taxon_id,t1.orthogroup,t2.ec ' \
                   'from (select orthogroup,locus_tag,ec_id from enzyme.locus2ec_%s ' \
                   'where orthogroup in (%s)) t1 ' \
                   'left join enzyme.enzymes as t2 on t1.ec_id=t2.enzyme_id ' \
                   'left join biosqldb.orthology_detail_%s as t3 ' \
                   'on t1.locus_tag=t3.locus_tag;' % (biodb,'"'+'","'.join(set(orthogroup_list))+'"', biodb)
        else:
            taxon2orthogroup2count_reference = get_taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb, [fam], 'ko')
            sql3 = 'select distinct A.taxon_id,A.orthogroup,B.ko_id from (' \
                   ' select locus_tag,orthogroup,taxon_id from biosqldb.orthology_detail_%s ' \
                   ' where orthogroup in (%s)) A inner join enzyme.locus2ko_%s as B ' \
                   ' on A.locus_tag=B.locus_tag;' % (biodb,'"'+'","'.join(set(orthogroup_list))+'"', biodb)

        print "sql3", sql3
        data = server.adaptor.execute_and_fetchall(sql3,)
        taxon2orthogroup2ec = {}
        for one_row in data:
            taxon = one_row[0]
            group = one_row[1]
            ec = one_row[2]
            if taxon not in taxon2orthogroup2ec:
                taxon2orthogroup2ec[taxon] = {}
                taxon2orthogroup2ec[taxon][group] = [ec]
            else:
                if group not in taxon2orthogroup2ec[taxon]:
                    taxon2orthogroup2ec[taxon][group] = [ec]
                else:
                    if ec not in taxon2orthogroup2ec[taxon][group]:
                        taxon2orthogroup2ec[taxon][group].append(ec)




        taxon2orthogroup2count = ete_motifs.get_taxon2orthogroup2count(biodb, group_count)
        print "taxon2orthogroup2count", taxon2orthogroup2count
        merged_dico = taxon2orthogroup2count
        for i in taxon2orthogroup2count_reference:
            merged_dico[i] = taxon2orthogroup2count_reference[i]
        print 'merged dico', merged_dico
        labels = [fam] + group_count
        #try:



        tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, merged_dico, taxon2group2value=taxon2orthogroup2ec,highlight_first_column=True)
        #except:
        #    tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, merged_dico)
        print tree
        if len(labels) > 30:
            big = True
            path = settings.BASE_DIR + '/assets/temp/fam_tree_%s.png' % fam
            asset_path = '/assets/temp/cog_tree_%s.png' % fam
            tree.render(path, dpi=1200, h=600)
        else:
            big = False
            path = settings.BASE_DIR + '/assets/temp/fam_tree_%s.svg' % fam
            asset_path = '/assets/temp/fam_tree_%s.svg' % fam

            tree.render(path, dpi=800, h=600)

    return render(request, 'chlamdb/fam.html', locals())

@login_required
def KEGG_module_map(request, biodb, module_name):

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST
        import ete_motifs
        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select module_sub_cat,module_sub_sub_cat,description,ko_id,ko_description ' \
              ' from enzyme.module2ko as t1 inner join enzyme.kegg_module as t2 on t1.module_id=t2.module_id' \
              ' where module_name="%s";' % (module_name)
        map_data = server.adaptor.execute_and_fetchall(sql,)

        print map_data

        ko_list = [i[3] for i in map_data]

        # get list of all orthogroups with corresponding ko
        sql = 'select distinct orthogroup from enzyme.locus2ko_%s where ko_id in (%s);' % (biodb, '"' + '","'.join(ko_list) + '"')
        orthogroup_data = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
        print orthogroup_data
        ko2orthogroups = {}
        orthogroup_list = []
        for i in orthogroup_data:
            if i not in ko2orthogroups:
                ko2orthogroups[i] = [i]
            else:
                ko2orthogroups[i].append(i)
            orthogroup_list.append(i)
        print ko2orthogroups
        taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup")
        taxon2ko2count = ete_motifs.get_taxon2name2count(biodb, ko_list, type="ko")

        print taxon2ko2count
        labels = ko_list
        tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2ko2count)

        tree2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                     labels,
                                                     taxon2orthogroup2count,
                                                     taxon2ko2count,
                                                     ko2orthogroups)


        if len(labels) > 40:
            print 'BIGGGGGGGGGGG', len(labels)
            big = True
            path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.png' % module_name
            asset_path = '/assets/temp/KEGG_tree_%s.png' % module_name
            tree.render(path, dpi=1200, h=600)



        else:
            print 'not BIGGGGGGGGGG', len(labels)
            big = False
            path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.svg' % module_name
            asset_path = '/assets/temp/KEGG_tree_%s.svg' % module_name
            tree.render(path, dpi=800, h=600)

            path2 = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s_complete.svg' % module_name
            asset_path2 = '/assets/temp/KEGG_tree_%s_complete.svg' % module_name

            tree2.render(path2, dpi=800, h=600)
        envoi = True
        menu = True
        valid_id = True


    return render(request, 'chlamdb/KEGG_module_map.html', locals())



@login_required
def KEGG_mapp(request, biodb, map_name):

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST
        import ete_motifs
        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select pathway_name,pathway_category,description,ec,value from ' \
              '(select pathway_name,pathway_category,description,ec_id from enzyme.kegg_pathway as t1 ' \
              'inner join enzyme.kegg2ec as t2 on t1.pathway_id=t2.pathway_id where pathway_name="%s") A ' \
              'inner join enzyme.enzymes as B on A.ec_id=B.enzyme_id inner join enzyme.enzymes_dat on enzymes_dat.enzyme_dat_id=enzyme_id ' \
              'where line="description";' % (map_name)
        map_data = server.adaptor.execute_and_fetchall(sql,)

        print map_data

        enzyme_list = [i[3] for i in map_data]

        sql = 'select id from comparative_tables.EC_%s where id in (%s);' % (biodb,
                                                          '"' + '","'.join(enzyme_list) + '"')
        enzyme_list_found_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

        # get list of all orthogroups with corresponding EC
        sql = 'select distinct ec,orthogroup from enzyme.locus2ec_%s as t1 ' \
              ' inner join enzyme.enzymes as t2 on t1.ec_id=t2.enzyme_id where ec in (%s);' % (biodb,
                                                          '"' + '","'.join(enzyme_list_found_in_db) + '"')
        orthogroup_data = server.adaptor.execute_and_fetchall(sql,)
        ec2orthogroups = {}
        orthogroup_list = []
        for i in orthogroup_data:
            if i[0] not in ec2orthogroups:
                ec2orthogroups[i[0]] = [i[1]]
            else:
                ec2orthogroups[i[0]].append(i[1])
            orthogroup_list.append(i[1])

        taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup")
        taxon2enzyme2count = ete_motifs.get_taxon2name2count(biodb, enzyme_list_found_in_db, type="EC")

        print taxon2enzyme2count
        labels = enzyme_list_found_in_db
        tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2enzyme2count)

        tree2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                     labels,
                                                     taxon2orthogroup2count,
                                                     taxon2enzyme2count,
                                                     ec2orthogroups)


        if len(labels) > 40:
            print 'BIGGGGGGGGGGG', len(labels)
            big = True
            path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.png' % map_name
            asset_path = '/assets/temp/KEGG_tree_%s.png' % map_name
            tree.render(path, dpi=1200, h=600)



        else:
            print 'not BIGGGGGGGGGG', len(labels)
            big = False
            path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.svg' % map_name
            asset_path = '/assets/temp/KEGG_tree_%s.svg' % map_name
            tree.render(path, dpi=800, h=600)

            path2 = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s_complete.svg' % map_name
            asset_path2 = '/assets/temp/KEGG_tree_%s_complete.svg' % map_name

            tree2.render(path2, dpi=800, h=600)
        envoi = True
        menu = True
        valid_id = True


    return render(request, 'chlamdb/KEGG_map.html', locals())

@login_required
def sunburst(request, biodb, locus):

    cache = get_cache('default')

    #cache.clear()

    if request.method == 'GET':  # S'il s'agit d'une requête POST
        print 'okkkkkkkkkkkkkkk'
        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession from orthology_detail_%s where locus_tag = "%s"' % (biodb, locus)
        accession = server.adaptor.execute_and_fetchall(sql,)[0][0]
        print accession
        try:
            sql1 = 'select t3.superkingdom,  t3.phylum,  t3.order,  t3.family,  t3.genus,  t3.species  from ' \
                   ' blastnr.blastnr_hits_%s_%s as t1  inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s ' \
                   ' as t2 on t1.nr_hit_id = t2.nr_hit_id  inner join blastnr.blastnr_taxonomy as t3 on ' \
                   ' t2.subject_taxon_id = t3.taxon_id inner join blastnr.blastnr_hsps_%s_%s as t4 ' \
                   ' on t1.nr_hit_id=t4.nr_hit_id where t1.locus_tag="%s"' % (biodb, accession, biodb, accession, biodb, accession, locus)
            print sql
            raw_data = server.adaptor.execute_and_fetchall(sql1,)

        except:
            print sql1
            valid_id = False
            return render(request, 'chlamdb/sunburst.html', locals())

        print "asdffffffffffffffffffffff", sql1
        dico =  {}

        y = 1
        for data in raw_data:

            if data[0] not in dico:
                dico[data[0]] = {}
            else:
                if data[1] not in dico[data[0]]:
                    dico[data[0]][data[1]] = {}

                if data[2] not in dico[data[0]][data[1]]:
                    dico[data[0]][data[1]][data[2]] = {}

                if data[3] not in dico[data[0]][data[1]][data[2]]:
                    dico[data[0]][data[1]][data[2]][data[3]] = {}

                if data[4] not in dico[data[0]][data[1]][data[2]][data[3]]:
                    dico[data[0]][data[1]][data[2]][data[3]][data[4]] = {}

                if data[5] not in dico[data[0]][data[1]][data[2]][data[3]][data[4]]:
                    dico[data[0]][data[1]][data[2]][data[3]][data[4]][data[5]] = 1

                dico[data[0]][data[1]][data[2]][data[3]][data[4]][data[5]] += 1

                y += 1
        i = 1

        tt = settings.BASE_DIR + '/assets/out.tab'
        print 'sdsdfsdf', tt
        out = open(tt, 'w')
        for superkingdom in dico:
            for phylum in dico[superkingdom]:
                for order in dico[superkingdom][phylum]:
                    for family in dico[superkingdom][phylum][order]:
                        for genus in dico[superkingdom][phylum][order][family]:
                            for species in dico[superkingdom][phylum][order][family][genus]:
                                out.write('%s, %s, %s, %s\n' % (i, 1, superkingdom, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 2, phylum, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 3, order, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 4, family, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 5, genus, 0))
                                out.write("%s, %s, %s, %s\n" % (i, 6, species, dico[superkingdom][phylum][order][family][genus][species]))
                                i+=1

                                print '\'%s, %s, %s, %s\\n\' +' % (i, 1, superkingdom, 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 2, phylum, 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 3, order, 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 4, family, 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 5, genus, 0)
                                print "\'%s, %s, %s, %s\\n\' +" % (i, 6, species, dico[superkingdom][phylum][order][family][genus][species])
        out.close()

        envoi = True
        menu = True


    return render(request, 'chlamdb/sunburst.html', locals())



@login_required
def blastnr(request, biodb, locus_tag):


    print biodb, locus_tag

    cache = get_cache('default')
    print "cache", cache
    #cache.clear()

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    if request.method == 'GET':  # S'il s'agit d'une requête POST



        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession, organism from orthology_detail_%s where locus_tag="%s"' % (biodb, locus_tag)
        data = server.adaptor.execute_and_fetchall(sql,)[0]
        accession = data[0]
        organism = data[1]

        server, db = manipulate_biosqldb.load_db(biodb)

        columns = 't1.hit_number, t1.locus_tag, t1.query_accession, t1.subject_accession, t1.subject_taxid, t1.subject_scientific_name,' \
                  't1.subject_title, t1.subject_kingdom, t2.evalue, t2.percent_identity, t2.gaps, t2.length,' \
                  't2.query_start, t2.query_end, t2.subject_start, t2.subject_end, t1.subject_title'
        sql = 'select %s from blastnr.blastnr_hits_%s_%s as t1 inner join blastnr.blastnr_hsps_%s_%s as t2 on' \
              ' t1.nr_hit_id = t2.nr_hit_id where  t1.locus_tag="%s" ' % (columns, biodb, accession, biodb, accession, locus_tag)
        print sql
        blast_data = list(server.adaptor.execute_and_fetchall(sql))
        #blast_data = [i for i in ]


        if len(blast_data) > 0:
            valid_id = True
            blast_query_locus = blast_data[0][1]
            blast_query_protein_id = blast_data[0][2]
            if blast_query_protein_id == blast_query_locus:
                 blast_query_protein_id = ''

            for n, one_hit in enumerate(blast_data):
                blast_data[n] = [i for i in one_hit]
                subject_taxids = one_hit[4].split(';')
                subject_scientific_names = one_hit[5].split(';')
                all_taxonomy = ''
                for taxon, name in zip(subject_taxids, subject_scientific_names):
                    all_taxonomy += '<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s">%s<a> ' % (taxon, name)
                blast_data[n][5] = all_taxonomy



        return render(request, 'chlamdb/blastnr.html', locals())


    return render(request, 'chlamdb/blastnr.html', locals())


@login_required
def homology(request, biodb):
    import shell_command
    cache = get_cache('default')
    print "cache", cache
    #cache.clear()

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    contact_form_class = make_contact_form(server, biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST


        form = contact_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            accession = request.POST['accession']

            envoi = True
    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = make_contact_form(server, biodb)  # Nous créons un formulaire vide
    return render(request, 'chlamdb/homology.html', locals())


@login_required
def orthogroup_identity(request, biodb, orthogroup, group=False):


    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    #if request.method == 'POST':
    import numpy
    import pandas as pd
    sql = 'SELECT * FROM orth_%s.%s;' % (biodb, orthogroup)

    try:
        data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])
        homologs = True
    except:
        homologs = False
        return render(request, 'chlamdb/orthogroup_identity.html', locals())
    locus_list = '"' + '","'.join(data[0:,1]) + '"'

    sql2 = 'select locus_tag, organism from orthology_detail_%s where locus_tag in (%s)' % (biodb, locus_list)
    print sql2
    locus2organism = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    print locus2organism

    #print data

    print data[0:,2:]

    columns = data[0:,1]
    rows = [i + " (%s)" % locus2organism[i] for i in columns]
    frame = pd.DataFrame(data[0:,2:], index=rows, columns=columns)
    frame = frame.astype(float)
    frame = frame/100

    for i in range(0, len(frame)):

        frame.ix[i, i] = None

    path = settings.BASE_DIR + '/assets/temp/%s.json' % orthogroup
    print 'writing to', path
    print frame
    with open(path, 'w') as f:
        f.write(frame.to_json(orient="split"))

    return render(request, 'chlamdb/orthogroup_identity.html', locals())


def ortho_id_plot(request, group):
    return render(request, 'chlamdb/orthogroup_identity_plot.html', locals())




@login_required
def plot_neighborhood(request, biodb, target_locus, region_size=23000):

    cache = get_cache('default')
    print "cache", cache

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    print "request.method", request.method

    valid_id = True

    server, db = manipulate_biosqldb.load_db(biodb)

    sql2 = 'select orthogroup, taxon_id from orthology_detail_%s where locus_tag = "%s"' % (biodb, target_locus)
    print sql2
    reference_orthogroup = server.adaptor.execute_and_fetchall(sql2, )[0]
    print "orthogroup", reference_orthogroup
    if not reference_orthogroup:
            valid_id = False
    if valid_id:
        orthogroup = reference_orthogroup[0]
        reference_taxon_id = reference_orthogroup[1]
        if plot_region:
            print "plotting!!!!!!!!!!!!!!"

            home_dir = os.path.dirname(__file__)
            print "home_dir", home_dir
            temp_location = os.path.join(home_dir, "../assets")
            print "temp loc", temp_location
            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
            print "temp file", temp_file.name
            name = os.path.basename(temp_file.name)
            print name.split('.')
            name_png = name.split('.')[0] + '.png'


            locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, biodb, [target_locus],
                                                                              temp_file.name, int(region_size),
                                                                              cache)

    envoi = True
    import ete_motifs
    taxon2orthogroup2count_all =  ete_motifs.get_taxon2orthogroup2count(biodb, orthogroup_list)

    labels = orthogroup_list

    n_orthogroup = orthogroup_list.index(orthogroup)
    print "####################### reference_taxon_id", reference_taxon_id
    tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2orthogroup2count_all, reference_taxon_id, n_orthogroup)

    '''
    #print tree
    if len(labels) > 30:
        big = True
        path = settings.BASE_DIR + '/assets/temp/cog_tree.png'
        asset_path = '/assets/temp/cog_tree.png'
        print path
        tree.render(path, dpi=1200, h=1200)
    else:
    '''
    big = False
    path = settings.BASE_DIR + '/assets/temp/cog_tree.svg'
    asset_path = '/assets/temp/cog_tree.svg'

    tree.render(path, dpi=800, h=600)
    print asset_path
    return render(request, 'chlamdb/plot_region_and_profile.html', locals())

@login_required
def plot_region(request, biodb):

    cache = get_cache('default')
    print "cache", cache

    #bioentry_in_memory = cache.get("biodb")
    print "loading db..."
    server = manipulate_biosqldb.load_db()
    print "db loaded..."
    plot_region_form_class = make_plot_form(biodb)
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = plot_region_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            valid_id = True

            accession = extract_alphanumeric(form.cleaned_data['accession'])

            server, db = manipulate_biosqldb.load_db(biodb)

            region_size = form.cleaned_data['region_size']

            genomes = form.cleaned_data['genomes']



            sql2 = 'select orthogroup, locus_tag, protein_id, start, stop, strand, organism from orthology_detail_%s where locus_tag like "%%%%%s%%%%" or protein_id like "%%%%%s%%%%"' % (biodb, accession, accession)
            print sql2
            data = server.adaptor.execute_and_fetchall(sql2, )[0]
            print "seqfeature_id", data
            if not data:
                    valid_id = False
            if valid_id:
                orthogroup = data[0]

                select = 'and (taxon_id = %s' % genomes[0]
                if len(genomes) >1:
                    for i in range(0, len(genomes)-1):
                        select+= ' or taxon_id = %s' % genomes[i]
                    select+= ' or taxon_id = %s)' % genomes[-1]
                else:
                    select+= ')'
                sql3 = 'select locus_tag from orthology_detail_%s where orthogroup = "%s" %s' % (biodb, orthogroup, select)
                print sql3
                locus_tag_target_genomes = [i[0] for i in server.adaptor.execute_and_fetchall(sql3, )]

                if plot_region:
                    print "plotting!!!!!!!!!!!!!!"
                    print "locus_tag_list", locus_tag_target_genomes
                    home_dir = os.path.dirname(__file__)
                    print "home_dir", home_dir
                    temp_location = os.path.join(home_dir, "../assets")
                    print "temp loc", temp_location
                    temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
                    print "temp file", temp_file.name
                    name = os.path.basename(temp_file.name)
                    print name.split('.')
                    name_png = name.split('.')[0] + '.png'

                    locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, biodb, locus_tag_target_genomes,
                                                                                      temp_file.name, int(region_size),
                                                                                      cache)

                    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                              'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'


                    sql_locus = 'locus_tag="%s"' % locus_tags[0]
                    for locus in range(1, len(locus_tags)):
                        sql_locus += ' or locus_tag="%s"' % locus_tags[locus]

                    sql = 'select %s from orthology_detail_%s where %s' % (columns, biodb, sql_locus)
                    print sql

                    raw_data = server.adaptor.execute_and_fetchall(sql,)

                    n = 1
                    search_result = []
                    for one_hit in raw_data:
                        search_result.append((n,) + one_hit)
                        n+=1
                        print n
                    print search_result


            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = plot_region_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/plot_region.html', locals())

'''
@login_required
def plot_region(request, biodb):
    plot_form_class = make_plot_form(biodb)
    cache = get_cache('default')
    print "cache", cache
    server = manipulate_biosqldb.load_db()

    if request.method == 'POST':  # S'il s'agit d'une requête POST


        form = plot_form_class(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides






            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = plot_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/plot_region.html', locals())
'''


@login_required
def comparative_extract():
    '''
    get features present in genomes X Y Z and not in A B C

    Se baser sur la table d'orthologie pour identifier les genes

    SQL
    select orthogroup from orthology_chlam where `taxon X` > 1 and `taxon B` = 0

    :return:
    '''


@login_required
def orthogroups(request):
    if request.method == 'POST':
        form = BiodatabaseForm(request.POST)
        if form.is_valid():
            biodb = form.cleaned_data['biodatabase']
            groups = "orthogroup_size_distrib_%s.svg" % biodb
            envoi = True
    else:
        form = BiodatabaseForm()

    return render(request, 'chlamdb/orthogroups.html', locals())






def circos_main(request, biodb):

    import gbk2circos
    import circos

    cache = get_cache('default')
    server, db = manipulate_biosqldb.load_db(biodb)

    reference_taxon = int(request.GET.getlist('ref')[0])
    target_taxons = [int(i) for i in request.GET.getlist('t')]
    highlight = request.GET.getlist('h')

    remove = [
    "group_2998"
,"group_4199"
,"group_4358"
,"group_4361"
,"group_4405"
,"group_4411"
,"group_4426"
,"group_4429"
,"group_4432"
,"group_4434"
,"group_4453"
,"group_4455"
,"group_4457"
,"group_4459"
,"group_4469"
,"group_4471"
,"group_4475"
,"group_4476"
,"group_4483"
,"group_4485"
,"group_4489"
,"group_4490"
,"group_4493"
,"group_4494"
,"group_4495"
,"group_4496"
,"group_4500"
,"group_4503"
,"group_4511"
,"group_4512"
,"group_4518"
,"group_4519"
,"group_4520"
,"group_4522"
,"group_4523"
,"group_4524"
,"group_4525"
,"group_4527"
,"group_4537"
,"group_4540"
,"group_4544"
,"group_4545"
,"group_4547"
,"group_4575"
,"group_4576"
,"group_4584"
,"group_4595"
,"group_4598"
,"group_4601"
,"group_4624"
,"group_4632"
,"group_4633"
,"group_4635"
,"group_4641"
,"group_4642"
,"group_4646"
,"group_4655"
,"group_4657"
,"group_4667"
,"group_4671"
,"group_4678"
,"group_4684"
,"group_4704"
,"group_4713"
,"group_4740"
,"group_4750"
,"group_4758"
,"group_4760"
,"group_4775"
,"group_4782"
,"group_4797"
,"group_3579"
,"group_2518"
,"group_2460"
,"group_3289"
,"group_4663"
,"group_4697"
,"group_4507"
,"group_4515"
,"group_4347"
,"group_4546"
,"group_3735"
,"group_4458"
,"group_4468"
,"group_4543"
,"group_4499"
,"group_4600"
,"group_4466"
,"group_4654"
,"group_4536"
,"group_4366"
,"group_4706"
,"group_4449"
,"group_4541"
,"group_4637"
,"group_4682"
,"group_4699"
,"group_4670"
,"group_4651"
,"group_4718"
,"group_4787"
,"group_4474"
,"group_4532"
,"group_4542"
,"group_4588"
,"group_4650"
,"group_4660"
,"group_4666"
,"group_4738"
,"group_4755"
,"group_4526"
,"group_4534"
,"group_4581"
,"group_83"
,"group_4533"
,"group_4668"
,"group_4516"
,"group_4574"
,"group_4653"
,"group_4470"
,"group_4510"
,"group_4492"
,"group_4529"
,"group_4652"
,"group_4473"
,"group_4482"
,"group_4488"
,"group_4538"
,"group_4548"
,"group_4730"
,"group_4253"
,"group_4481"
,"group_4631"
,"group_4674"
,"group_4702"
,"group_4774"
,"group_4487"
,"group_4514"
,"group_4528"
,"group_4530"
,"group_4531"
,"group_4580"
,"group_4604"
,"group_4639"
,"group_4644"
,"group_4680"
,"group_4707"
,"group_4648"
,"group_4701"
,"group_4737"
,"group_4416"
,"group_4480"
,"group_4696"
,"group_4521"
,"group_4672"
,"group_4638"
,"group_4716"
,"group_4636"
,"group_4593"
,"group_4677"
,"group_4664"
,"group_4669"
,"group_2714"
,"group_4732"
,"group_4517"
,"group_4766"
,"group_4695"
,"group_4739"
,"group_4478"
,"group_4502"
,"group_4513"
,"group_4662"
,"group_4472"
,"group_4498"
,"group_4535"
,"group_4417"
,"group_4717"
,"group_4332"
,"group_4460"
,"group_4578"
,"group_4757"
,"group_4506"
,"group_4754"
,"group_4594"
,"group_4708"
,"group_4656"
,"group_4683"
,"group_4659"
,"group_4756"
,"group_4477"
,"group_4640"
,"group_4721"
,"group_4508"
,"group_4681"
,"group_4589"
,"group_4705"
,"group_4221"
,"group_4675"
,"group_4679"
,"group_4735"
,"group_4665"
,"group_4431"
,"group_4794"
,"group_4645"
,"group_4596"]



    highlight_def = []
    for i in highlight:
        if i in remove:
            continue
        else:
            highlight_def.append(i)

    print '"' + '","'.join(highlight_def) + '"'
    import time
    time.sleep(10)

    #sql = 'select locus_tag,traduction from orthology_detail_k_cosson_05_16 where orthogroup in (%s) and accession="NC_016845"' % ('"'+'","'.join(highlight)+'"')
    #print sql
    #import time
    #time.sleep(20)
    description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

    reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb) # ["NC_009648"] NZ_CP009208 NC_016845

    print "reference_accessions", reference_accessions
    record_list = []
    for accession in reference_accessions:
        print "reference accession", accession
        biorecord = cache.get(biodb + "_" + accession)

        if not biorecord:
            print biodb + "_" + accession, "NOT in memory"
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
            print biodb + "_" + accession, "In memory"
            record_list.append(biorecord)

    ref_name = ('').join(reference_accessions)

    circos_file = "circos/%s.svg" % ref_name

    querries = manipulate_biosqldb.get_genome_accessions(server, biodb)

    target_accessions = [manipulate_biosqldb.taxon_id2accessions(server,int(i),biodb)[0] for i in target_taxons]

    target_accessions += reference_accessions



    draft_data = []
    for biorecord in record_list:
        draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

    home_dir = os.path.dirname(__file__)

    temp_location = os.path.join(home_dir, "../assets/circos/")

    #sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    sql_order = 'select A.taxon_1 from (select taxon_1,median_identity from comparative_tables.shared_orthogroups_average_identity_%s where taxon_2=%s ' \
                ' union select taxon_2,median_identity from comparative_tables.shared_orthogroups_average_identity_%s ' \
                ' where taxon_1=%s order by median_identity DESC) A;' % (biodb, reference_taxon, biodb, reference_taxon)

    sql_order = 'select taxon_2 from comparative_tables.core_orthogroups_identity_msa_%s where taxon_1=%s order by identity desc;' % (biodb, reference_taxon)

    print sql_order
    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]

    '''
    print tree
    t1 = ete2.Tree(tree)

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
    print 'ordered_taxons', ordered_taxons

    myplot = circos.CircosAccession2multiplot(server,
                              db,
                              biodb,
                              record_list,
                              target_accessions,
                              locus_highlight=highlight_def,
                              out_directory=temp_location,
                              draft_fasta=draft_data,
                              href="/chlamdb/locusx/%s/" % biodb,
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

    envoi_circos = True


    envoi_region = True

    return render(request, 'chlamdb/circos.html', locals())




@login_required
def circos(request, biodb):

    import gbk2circos
    circos_form_class = make_circos_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)

    cache = get_cache('default')

    if request.method == 'POST':

        form = circos_form_class(request.POST)

        if form.is_valid():
            reference_taxon = form.cleaned_data['circos_reference']

            description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)

            print "reference_accessions", reference_accessions
            record_list = []
            for accession in reference_accessions:

                print "reference accession", accession
                biorecord = cache.get(biodb + "_" + accession)

                if not biorecord:
                    print biodb + "_" + accession, "NOT in memory"
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
                    print biodb + "_" + accession, "IN memory"
                    record_list.append(biorecord)


            if 'submit_circos' in request.POST:

                ref_name = ''
                for i in reference_accessions:
                    ref_name += i
                circos_file = "circos/%s.svg" % ref_name
                import circos
                import shell_command
                import ete2


                querries = manipulate_biosqldb.get_genome_accessions(server, biodb)
                target_taxons = form.cleaned_data['targets']



                target_accessions = [manipulate_biosqldb.taxon_id2accessions(server,int(i),biodb)[0] for i in target_taxons]

                print 'targets!', target_accessions

                target_accessions += reference_accessions
                print target_accessions

                draft_data = []
                for biorecord in record_list:
                    draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

                home_dir = os.path.dirname(__file__)
                print "home_dir", home_dir
                temp_location = os.path.join(home_dir, "../assets/circos/")

                #sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

                sql_order = 'select A.taxon_1 from (select taxon_1,median_identity from comparative_tables.shared_orthogroups_average_identity_%s where taxon_2=%s ' \
                            ' union select taxon_2,median_identity from comparative_tables.shared_orthogroups_average_identity_%s ' \
                            ' where taxon_1=%s order by median_identity DESC) A;' % (biodb, reference_taxon, biodb, reference_taxon)

                sql_order = 'select taxon_2 from comparative_tables.core_orthogroups_identity_msa_%s where taxon_1=%s order by identity desc;' % (biodb, reference_taxon)
                print sql_order
                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                '''
                print tree
                t1 = ete2.Tree(tree)

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
                print 'ordered_taxons', ordered_taxons

                myplot = circos.CircosAccession2multiplot(server,
                                          db,
                                          biodb,
                                          record_list,
                                          target_accessions,
                                          locus_highlight=[],
                                          out_directory=temp_location,
                                          draft_fasta=draft_data,
                                          href="/chlamdb/locusx/%s/" % biodb,
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

                envoi_circos = True


            envoi_region = True
    else:
        form = circos_form_class()
    return render(request, 'chlamdb/circos.html', locals())

@login_required
def alignment(request, input_fasta):
    print align
    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        pass
    return render(request, 'chlamdb/alignment.html', locals())






def format_seqfeature_values(server, biodb, seqfeature_id):


    seqfeature_data = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodb)
    if not 'translation' in seqfeature_data.keys():
        # TODO add handeling of other kind of features than CDS
        return None
    try:
        sql_family_size = 'select count(*) as `n_rows` from ' \
           ' (select * from orthology_detail_chlamydia_02_15 ' \
           ' where orthogroup = "%s" group by taxon_id) a' % (seqfeature_data["orthogroup"])
        sql_n_homologues = 'select count(*) as `n_rows` from ' \
           ' (select * from orthology_detail_chlamydia_02_15 ' \
           ' where orthogroup = "%s") a' % (seqfeature_data["orthogroup"])
        n_family = int(server.adaptor.execute_and_fetchall(sql_family_size,)[0][0])
        n_homologues = int(server.adaptor.execute_and_fetchall(sql_n_homologues,)[0][0])
    except KeyError:
        n_family = '-'
        n_homologues = '-'



    # 0 seqfeature_values["description"],
    # 1 seqfeature_values["locus_tag"],
    # 2 seqfeature_values["protein_id"],
    # 3 seqfeature_values["product"],
    # 4 seqfeature_values["orthogroup"],
    # 5 seqfeature_values["gene"],
    # 6 seqfeature_values["translation"],
    # 7 n_family,
    # 8 n_homologues
    template = ["-"] * 9
    try:
        template[0] = seqfeature_data["description"]
    except KeyError:
        pass
    try:
        template[1] = seqfeature_data["locus_tag"]
    except KeyError:
        pass
    try:
        template[2] = seqfeature_data["protein_id"]
    except KeyError:
        pass
    try:
        template[3] = seqfeature_data["product"]
    except KeyError:
        pass
    try:
        template[4] = seqfeature_data["orthogroup"]
    except KeyError:
        pass
    try:
        template[5] = seqfeature_data["gene"]
    except KeyError:
        pass
    try:
        template[6] = seqfeature_data["translation"]
    except KeyError:
        pass

    template[7] = n_family
    template[8] = n_homologues

    return template
    '''
        try:
            return []
        except KeyError:
            try:
                return [seqfeature_values["description"], "-", seqfeature_values["protein_id"], seqfeature_values["product"], seqfeature_values["orthogroup"], seqfeature_values["gene"], seqfeature_values["translation"], n_family, n_homologues]
            except KeyError:
                try:
                    return [seqfeature_values["description"], seqfeature_values["locus_tag"] , "-", seqfeature_values["product"], seqfeature_values["orthogroup"], seqfeature_values["gene"], seqfeature_values["translation"], n_family, n_homologues]
                except KeyError:
                    try:
                        return [seqfeature_values["description"], "-" , "-", seqfeature_values["product"], seqfeature_values["orthogroup"], seqfeature_values["gene"], seqfeature_values["translation"], n_family, n_homologues]
                    except KeyError:
                        try:
                            return [seqfeature_values["description"], seqfeature_values["locus_tag"] , seqfeature_values["protein_id"], seqfeature_values["product"], seqfeature_values["orthogroup"], "-", seqfeature_values["translation"], n_family, n_homologues]
                        except KeyError:
                            return [seqfeature_values["description"], seqfeature_values["locus_tag"] , seqfeature_values["protein_id"], "-", seqfeature_values["orthogroup"], "-", seqfeature_values["translation"], n_family, n_homologues]
    '''




def format_search(count, seqfeature_data):
    # [y, i["description"], i["locus_tag"] + " / " + i["protein_id"], i["product"], i["orthogroup"],  i["gene"], i["translation"]]
    pass


@login_required
def search_taxonomy(request, biodb):
    from collections import Counter
    server = manipulate_biosqldb.load_db()

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        print request.POST

        if len(request.POST['Phylum']) == 0:
            genome_accession = request.POST['Genome']
            superkingdom = request.POST['Superkingdom'].split('_')[-1]
            print "genome_accession", genome_accession
            print 'sdfsdf', superkingdom[0]
            if 'unclassified' in superkingdom:
                superkingdom = '-'

            if len(superkingdom) == 0:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s"' % (biodb,
                                                         genome_accession,
                                                         biodb,
                                                         genome_accession,
                                                         biodb,
                                                         genome_accession,
                                                         1)

            else:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s"' % (biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         1,
                                                                         superkingdom)

        else:
            genome_accession = request.POST['Genome']
            superkingdom = request.POST['Superkingdom'].split('_')[-1]

            if request.POST['Phylum'] == 'Unclassified':
                phylum = '-'
            else:
                phylum = request.POST['Phylum'].split('_')[-1]

            if len(phylum) ==0:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id ' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s"' % (biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         1,
                                                                         superkingdom)
            else:
                sql= 'select t1.locus_tag, ' \
                          ' t2.subject_taxon_id,' \
                          ' t4.nr_hit_id,' \
                          ' t1.subject_accession, ' \
                          ' t3.kingdom, ' \
                          ' t3.phylum, ' \
                          ' t3.order, ' \
                          ' t3.family, ' \
                          ' t3.genus, ' \
                          ' t3.species, ' \
                          ' t4.evalue, ' \
                          ' t4.percent_identity, ' \
                          ' t4.query_start, ' \
                          ' t4.query_end, ' \
                          ' t1.subject_title, ' \
                          ' t3.taxon_id,' \
                          ' t3.superkingdom' \
                          ' from blastnr.blastnr_hits_%s_%s as t1 ' \
                          ' inner join blastnr.blastnr_hits_taxonomy_filtered_%s_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr.blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr.blastnr_hsps_%s_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s" and t3.phylum="%s";' % (biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         biodb,
                                                                         genome_accession,
                                                                         1,
                                                                         superkingdom,
                                                                         phylum)
        print sql
        raw_data = server.adaptor.execute_and_fetchall(sql, )

        n = 1
        data = []
        families = []
        phylum = []
        superkingdom = []
        order = []
        for i, one_hit in enumerate(raw_data):
            if n == 1:
                data.append(one_hit + (n,))
                families.append(one_hit[7])
                phylum.append(one_hit[5])
                superkingdom.append(one_hit[-1])
                n+=1
            else:
                if raw_data[i][0] == raw_data[i-1][0]:
                    data.append(one_hit + (n,))
                else:
                    n+=1
                    data.append(one_hit + (n,))
                    families.append(one_hit[7])
                    phylum.append(one_hit[5])
                    superkingdom.append(one_hit[-1])
        if not 'Phylum' in request.POST and not 'Superkingdom' in request.POST:
            classif_table = dict(Counter(superkingdom))

        elif not 'Phylum' in request.POST and 'Superkingdom' in request.POST:
            classif_table = dict(Counter(phylum))

        else:
            classif_table = dict(Counter(families))

        print classif_table

        '''
        import pandas
        frame = pandas.DataFrame(data, columns= ['n',
                                                 'locus_tag',
                                                 'subject_taxon_id',
                                                 'nr_hit_id',
                                                 'subject_accession',
                                                 'kingdom',
                                                 'phylum',
                                                 'order',
                                                 'family',
                                                 'genus',
                                                 'species',
                                                 'evalue',
                                                 'percent_identity','query_start','query_end','subject_title','taxon_id'])

        print frame
        '''
        envoi = True


    return render(request, 'chlamdb/search_taxonomy.html', locals())


@login_required
def interpro(request, biodb):
    server = manipulate_biosqldb.load_db()

    interproform = make_interpro_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = interproform(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            invalid_id = False
            # Ici nous pouvons traiter les données du formulaire
            search_type = form.cleaned_data['search_type']
            search_term = form.cleaned_data['search_term']
            taxon_ids = form.cleaned_data['targets']
            #biodb = form.cleaned_data['biodatabase']
            server, db = manipulate_biosqldb.load_db(biodb)

            columns = 'accession,' \
                      'locus_tag,' \
                      'organism, ' \
                      'analysis, ' \
                      'signature_accession, ' \
                      'signature_description, ' \
                      'interpro_accession, ' \
                      'interpro_description,' \
                      'start, ' \
                      'stop, ' \
                      'score, ' \
                      'GO_terms'

            if len(taxon_ids) == 0:
                invalid_id = True
            else:

                taxon_limit = '(taxon_id=%s' % taxon_ids[0]
                if len(taxon_ids) > 1:
                    for i in range(1, len(taxon_ids)-1):
                        taxon_limit+= ' or taxon_id=%s' % taxon_ids[i]
                    taxon_limit+=' or taxon_id=%s)' % taxon_ids[-1]
                else:
                    taxon_limit += ')'

                if search_type == "description":

                    sql = 'select %s from interpro_%s where %s and (interpro_description REGEXP "%s" or signature_description REGEXP "%s")' % (columns, biodb, taxon_limit, search_term, search_term)

                if search_type == "GO":
                    sql = 'select %s from interpro_%s where %s and (GO_terms REGEXP "%s")' % (columns, biodb, taxon_limit, search_term)


                if search_type == "EC":
                    sql = 'select %s from interpro_%s where %s and (pathways REGEXP "%s")' % (columns, biodb, taxon_limit, search_term)

                if search_type == "interpro_accession":
                    sql = 'select %s from interpro_%s where %s and (interpro_accession REGEXP "%s")' % (columns, biodb, taxon_limit, search_term)


                try:
                    raw_data = server.adaptor.execute_and_fetchall(sql, )
                except IndexError:
                    invalid_id = True


            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = interproform()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/interpro.html', locals())

@login_required
def search(request, biodb):
    server = manipulate_biosqldb.load_db()
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = SearchForm(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            invalid_id = False
            # Ici nous pouvons traiter les données du formulaire
            search_type = form.cleaned_data['search_type']
            search_term = form.cleaned_data['search_term']
            #biodb = form.cleaned_data['biodatabase']
            server, db = manipulate_biosqldb.load_db(biodb)
            print "biodb", biodb

            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'

            if search_type == "gene":

                sql = 'select %s from orthology_detail_%s where gene REGEXP "%s"' % (columns, biodb, search_term)
                raw_data = server.adaptor.execute_and_fetchall(sql,)

            if search_type == "product":
                sql = 'select %s from orthology_detail_%s where product REGEXP "%s"' % (columns, biodb, search_term)
                raw_data = server.adaptor.execute_and_fetchall(sql,)

            if search_type == "locus_tag":
                sql = 'select %s from orthology_detail_%s where locus_tag REGEXP "%s"' % (columns, biodb, search_term)
                raw_data = server.adaptor.execute_and_fetchall(sql,)


            n = 1
            search_result = []
            for one_hit in raw_data:
                if one_hit[2] != '-':
                    interpro_id = one_hit[2]
                else:
                    interpro_id = one_hit[1]
                search_result.append((n,) + one_hit + (interpro_id,))
                n+=1
                print n
            print search_result


            '''
            seqfeature_data = []
            for one_id in seqfeature_ids:
                print one_id
                seqfeature_values = manipulate_biosqldb.seqfeature_id2seqfeature_qualifier_values(server, one_id, biodb)
                seqfeature_data.append(seqfeature_values)
            search_result = []
            y = 1
            for i in seqfeature_data:
                if not 'translation' in i.keys():
                    continue
                print i.keys()
                print i
                #  ['locus_tag', 'orthogroup', 'transl_table', 'product', 'translation', 'gene']
                search_result.append(format_search(y, i))


                y+=1
            '''
            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = SearchForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/search.html', locals())



@login_required
def primer_search(request, biodb):
    server = manipulate_biosqldb.load_db()
    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = PCRForm(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            from Bio.Blast.Applications import NcbiblastpCommandline
            #from StringIO import StringIO
            from tempfile import NamedTemporaryFile


            from Bio.Alphabet import IUPAC
            import os
            import shell_command
            import re
            def ExtractAlphanumeric(InputString):
                from string import ascii_letters, digits
                return "".join([ch for ch in InputString if ch in (ascii_letters + digits)])

            input_sequence = form.cleaned_data['blast_input']
            input_sequence = ExtractAlphanumeric(input_sequence)
            print input_sequence

            #biodb = form.cleaned_data['biodatabase']
            input_sequence = input_sequence.rstrip(os.linesep)
            print input_sequence
            my_record = SeqRecord(Seq(input_sequence, IUPAC.protein), id="INPUT", description="INPUT")
            print my_record
            query_file = NamedTemporaryFile()
            SeqIO.write(my_record, query_file, "fasta")
            query_file.flush()

            blastdb = settings.BASE_DIR + '/assets/blast_db/%s.faa' % biodb


            blastp_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=0.001, outfmt=0)
            stdout, stderr = blastp_cline()
            print "blast!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print stdout
            blast_file = NamedTemporaryFile()
            blast_file.write(stdout)
            mview_cmd = 'mview -in blast -ruler on -html data -css on -coloring identity %s' % blast_file.name
            stdout, stderr, code = shell_command.shell_command(mview_cmd)
            blast_result = stdout

            #blast_result = NCBIXML.parse(StringIO(stdout))
            #print blast_result
            #blast_record = next(blast_result)
            #print blast_record



            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = PCRForm()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/pcr.html', locals())


@login_required
def motif_search(request, biodb):
    import shell_command
    import os
    server = manipulate_biosqldb.load_db()

    motif_form_class = make_motif_form(biodb)

    if request.method == 'POST':

        form = motif_form_class(request.POST)

        if form.is_valid():

            from Bio.Emboss.Applications import FuzznucCommandline
            from tempfile import NamedTemporaryFile

            input_pattern = form.cleaned_data['motif_input']
            n_missmatch = form.cleaned_data['n_missmatch']
            target_taxon_id = form.cleaned_data['search_in']
            print "target", target_taxon_id

            if target_taxon_id != "all":
                accessions = manipulate_biosqldb.taxon_id2accessions(server, target_taxon_id, biodb)


            #input_pattern = ExtractAlphanumeric(input_pattern)
            print input_pattern
            '''
            fuzzpro -sequence CHUV_chr_and_plasmid.faa
            -rformat gff
            -auto
            -stdout
            -pattern "K-x(2)-[LIVF]-x(4)-[LIVF]-D-x(3)-R-x(2)-L-x(5)-[LIV]-Y"
            '''

            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))



            if target_taxon_id != "all":
                db_path = os.path.join(PROJECT_ROOT,"../assets/%s/faa/" % biodb +accessions[0] + ".faa")
            else:
                db_path = os.path.join(PROJECT_ROOT,"../assets/%s/faa/all.faa" % biodb)
            print "path", db_path

            cmd = 'fuzzpro -sequence %s -pmismatch %s -rformat seqtable -auto -stdout -pattern "%s"' % (db_path ,n_missmatch ,input_pattern)
            print cmd
            std_out, std_err, code = shell_command.shell_command(cmd)
            print std_out
            print std_err
            #fuzznuc_cline = FuzznucCommandline(sequence=genome_db, mismatch=n_missmatch, pattern=input_pattern, stdout=True)#, rformat="srspair")
            #stdout,stderr = fuzznuc_cline()


            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = motif_form_class()  # empty form

    return render(request, 'chlamdb/motifs.html', locals())



@login_required
def blast(request, biodb):
    server = manipulate_biosqldb.load_db()

    blast_form_class = make_blast_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST

        form = blast_form_class(request.POST)  # Nous reprenons les données

        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            from Bio.Blast.Applications import NcbiblastpCommandline
            from Bio.Blast.Applications import NcbiblastnCommandline
            from Bio.Blast.Applications import NcbitblastnCommandline
            from Bio.Blast.Applications import NcbiblastxCommandline
            from tempfile import NamedTemporaryFile

            from Bio.Alphabet import IUPAC
            import os
            import shell_command
            import re


            input_sequence = form.cleaned_data['blast_input']

            target_accession = form.cleaned_data['target']

            blast_type = form.cleaned_data['blast']


            input_sequence = extract_alphanumeric(input_sequence)
            print "one", input_sequence

            input_sequence = input_sequence.rstrip(os.linesep)
            print "two", input_sequence

            my_record = SeqRecord(Seq(input_sequence, IUPAC.protein), id="INPUT", description="INPUT")

            query_file = NamedTemporaryFile()
            SeqIO.write(my_record, query_file, "fasta")
            query_file.flush()


            if blast_type=='blastn_ffn':
                blastdb = settings.BASE_DIR + "/assets/%s/ffn/%s.ffn" % (biodb, target_accession)
                blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
            if blast_type=='blastn_fna':
                blastdb = settings.BASE_DIR + "/assets/%s/fna/%s.fna" % (biodb, target_accession)
                blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
            if blast_type=='blastp':
                blastdb = settings.BASE_DIR + "/assets/%s/faa/%s.faa" % (biodb, target_accession)
                blast_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
            if blast_type=='tblastn':
                blastdb = settings.BASE_DIR + "/assets/%s/fna/%s.fna" % (biodb, target_accession)
                blast_cline = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                blast_cline2 = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=5)
            if blast_type=='blastx':
                blastdb = settings.BASE_DIR + "/assets/%s/faa/%s.faa" % (biodb, target_accession)
                blast_cline = NcbiblastxCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)


            print blast_cline
            blast_stdout, blast_stderr = blast_cline()

            if blast_type=='tblastn':
                from Bio.SeqUtils import six_frame_translations
                from StringIO import StringIO
                blast_stdout2, blast_stderr2 = blast_cline2()
                from Bio.Blast import NCBIXML
                blast_records = NCBIXML.parse(StringIO(blast_stdout2))
                all_data = []
                for record in blast_records:
                    for alignment in record.alignments:
                        accession = alignment.title.split(' ')[1]
                        sql = 'select description from bioentry where accession="%s" ' % accession
                        description = server.adaptor.execute_and_fetchall(sql,)[0][0]
                        for hsp in alignment.hsps:
                            start = hsp.sbjct_start
                            end = hsp.sbjct_end
                            length = end-start
                            #print 'seq for acc', accession, start, end,
                            leng = end-start

                            print 'end', 'start', end, start, end-start

                            seq = manipulate_biosqldb.location2sequence(server, accession, biodb, start, leng)
                            print seq
                            from Bio.Seq import reverse_complement, translate
                            anti = reverse_complement(seq)
                            comp = anti[::-1]
                            length = len(seq)
                            frames = {}
                            for i in range(0, 3):
                                fragment_length = 3 * ((length-i) // 3)
                                tem1 = translate(seq[i:i+fragment_length], 1)
                                frames[i+1] = '<span style="color: #181407;">%s</span><span style="color: #bb60d5;">%s</span><span style="color: #181407;">%s</span>' % (tem1[0:100], tem1[100:len(tem1)-99], tem1[len(tem1)-99:])
                                tmp2 = translate(anti[i:i+fragment_length], 1)[::-1]
                                frames[-(i+1)] = tmp2
                            all_data.append([accession, start, end, length, frames[1], frames[2], frames[3], frames[-1], frames[-2], frames[-3], description, seq])


            no_match = re.compile('.* No hits found .*', re.DOTALL)

            if no_match.match(blast_stdout):
                print "no blast hit"
                blast_no_hits = blast_stdout
            elif len(blast_stderr) != 0:
                print "blast error"
                blast_err = blast_stderr
            else:
                print "running mview"
                blast_file = NamedTemporaryFile()

                f = open('/tmp/blast.temp', 'w')
                f.write(blast_stdout)
                f.close()



                #blast_file.write(blast_stdout)
                out, err, code = shell_command.shell_command('cat /temp/blast.temp | wc -l')
                print 'n lines', out
                if blast_type=='blastp' or blast_type=='blastn_ffn':
                    mview_cmd = 'mview -in blast -srs on -ruler on -html data -css on -coloring identity /tmp/blast.temp' #% blast_file.name
                else:
                    mview_cmd = 'mview -in blast -ruler on -html data -css on -coloring identity /tmp/blast.temp'
                stdout, stderr, code = shell_command.shell_command(mview_cmd)

                if len(stdout) == 0:
                    blast_no_hits = blast_stdout
                    blast_result = None
                else:
                    blast_result = stdout
            #blast_result = NCBIXML.parse(StringIO(stdout))
            #print blast_result
            #blast_record = next(blast_result)
            #print blast_record

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = blast_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/blast.html', locals())


def get_record_from_memory(biodb, cache_obj, record_key, accession):

        biorecord = cache_obj.get(record_key)
        if not biorecord:
            print record_key, "NOT in memory"
            new_record = biodb.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
            cache_obj.set(record_key, biorecord)
        return biorecord





@login_required
def mummer(request, biodb):

    server = manipulate_biosqldb.load_db()
    mummer_form_class = make_mummer_form(biodb)

    cache = get_cache('default')

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        #make_circos2genomes_form
        plot = True
        form = mummer_form_class(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            server, db = manipulate_biosqldb.load_db(biodb)
            reference_taxon = form.cleaned_data['reference_genome']
            query_taxon = form.cleaned_data['query_genome']

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)
            query_accessions = manipulate_biosqldb.taxon_id2accessions(server, query_taxon, biodb)


            ref_accession = manipulate_biosqldb.taxon_id2chromosome_accession(server, biodb, reference_taxon)
            query_accession = manipulate_biosqldb.taxon_id2chromosome_accession(server, biodb, query_taxon)
            #result = result[0]
            print "mummer acc", ref_accession, query_accession




            print settings.STATIC_ROOT, type(settings.STATIC_ROOT)

            reference_path = settings.BASE_DIR + '/assets/%s/fna/%s.fna' % (biodb, ref_accession)
            query_path = settings.BASE_DIR + '/assets/%s/fna/%s.fna' % (biodb, query_accession)

            rand = id_generator(5)

            out_delta = settings.BASE_DIR + '/assets/temp/promer_%s' % rand
            out_plot = settings.BASE_DIR + '/assets/temp/promer_%s' % rand




            cmd1 = 'promer -l 2 -p %s %s %s' % (out_delta, reference_path, query_path)
            cmd2 = 'mummerplot -layout -small -png -p %s %s.delta' % (out_plot, out_delta)

            print cmd1

            from shell_command import shell_command

            out, err, log = shell_command(cmd1)
            out, err, log = shell_command(cmd2)

            plot_path = 'temp/promer_%s.png' % rand

            if not os.path.exists(settings.BASE_DIR + '/assets/' + plot_path):
                plot = False

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = mummer_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/mummer.html', locals())





@login_required
def circos2genomes(request, biodb):
    import circos
    import shell_command
    server = manipulate_biosqldb.load_db()
    circos2genomes_form_class = make_circos2genomes_form(biodb)

    cache = get_cache('default')

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        #make_circos2genomes_form

        form = circos2genomes_form_class(request.POST)
        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            valid_id = True
            server, db = manipulate_biosqldb.load_db(biodb)
            reference_taxon = form.cleaned_data['reference_genome']
            query_taxon = form.cleaned_data['query_genome']
            import re
            protein_locus_list = re.sub(" ", "", form.cleaned_data['locus_list'])
            protein_locus_list = list(protein_locus_list.split(","))

            #reference_taxon = manipulate_biosqldb.description2taxon_id(server, reference_genome, biodb)
            #query_taxon = manipulate_biosqldb.description2taxon_id(server, query_genome, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)
            query_accessions = manipulate_biosqldb.taxon_id2accessions(server, query_taxon, biodb)

            reference_records = []
            for accession in reference_accessions:
                reference_records.append(get_record_from_memory(db, cache, biodb + "_" + accession, accession))

            query_records = []
            for accession in query_accessions:

                query_records.append(get_record_from_memory(db, cache, biodb + "_" + accession, accession))

            print "genomes", reference_records, query_records

            orthogroup_list = []
            print "protein_locus_list", protein_locus_list[0], len(protein_locus_list[0])
            if len(protein_locus_list[0]) > 0:
                print 'okkkkkkkk'
                for protein in protein_locus_list:
                    print "protein", protein
                    sql = 'select orthogroup from orthology_detail_%s where protein_id="%s" or locus_tag="%s"' % (biodb, protein, protein)
                    print sql
                    try:
                        protein_group = server.adaptor.execute_and_fetchall(sql,)[0][0]
                        orthogroup_list.append(protein_group)
                    except IndexError:
                        valid_id = False

                    #protein_group = manipulate_biosqldb.locus_tag2orthogroup_id(server, protein, biodb)
            if valid_id:
                print "LOCUS:", orthogroup_list

                #accession2description = manipulate_biosqldb.accession2description_dict(server, biodb)

                #print "reference_genome", reference_genome
                taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
                reference_name = taxon_id2description[reference_taxon]
                query_name = taxon_id2description[query_taxon]

                accession2taxon_id = manipulate_biosqldb.accession2taxon_id(server, biodb)
                taxon_id_reference = accession2taxon_id[reference_accessions[0]]
                taxon_id_query = accession2taxon_id[query_accessions[0]]

                reference_n_orthogroups = manipulate_biosqldb.get_genome_number_of_orthogroups(server, biodb, taxon_id_reference)
                reference_n_proteins = manipulate_biosqldb.get_genome_number_of_proteins(server, biodb, taxon_id_reference)

                query_n_orthogroups = manipulate_biosqldb.get_genome_number_of_orthogroups(server, biodb, taxon_id_query)
                query_n_proteins = manipulate_biosqldb.get_genome_number_of_proteins(server, biodb, taxon_id_query)

                n_shared_orthogroups = manipulate_biosqldb.get_number_of_shared_orthogroups(server, biodb, taxon_id_reference, taxon_id_query)

                print "n_orthogroups", reference_n_orthogroups, query_n_orthogroups
                print "n_proteins", reference_n_proteins, query_n_proteins

                print "n_shared", n_shared_orthogroups

                print reference_name, query_name


                import circos

                path = settings.BASE_DIR + "/assets/circos"



                biplot = circos.CircosAccession2biplot(server, db, biodb, reference_records, query_records,
                                                       orthogroup_list, path)

                reference_file = "circos/%s" % biplot.reference_circos
                query_file = "circos/%s" % biplot.query_circos

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = circos2genomes_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/circos2genomes.html', locals())

@login_required
def circos2genomes_main(request):
    return render(request, 'chlamdb/circos2genomes_main.html', locals())


def update_db(server):
    biodb_list = manipulate_biosqldb.get_biodatabase_list(server)
    for biodb in biodb_list:
        update_genomes_db(server, biodb)

    print

def update_genomes_db(server, biodatabase_name):
    from models import Genome
    from models import Database
    from models import GenDB

    print "updating genome Database"

    database = Database.objects.get_or_create(db_name=biodatabase_name)[0]
    genome_list = manipulate_biosqldb.get_genome_description_list(server, biodatabase_name)
    for genome_description in genome_list:
        print "updating", genome_description
        genome = Genome.objects.get_or_create(genome_name=genome_description, database=database)[0]
        print "genome", genome
        GenDB.objects.get_or_create(database=database, ref_genome=genome, query_genome=genome, genome_name=genome.genome_name, database_name=database.db_name)
        print "OKKKKKKKKKKKK"




#def get_genomes(server, biodb_name):













@login_required
def crossplot(request):

    cache = get_cache('default')
    print "cache", cache
    cache.clear()
    bioentry_in_memory = cache.get('biodb')
    print "bioentry_in_memory", bioentry_in_memory
    if not bioentry_in_memory:
        print "creating cache entry"
        cache.set("biodb", {})
    bioentry_in_memory = cache.get("biodb")
    server = manipulate_biosqldb.load_db()

    update_db(server)

    crossplot_form_class = make_crossplot_form("Chlamydia_11_14")
    if request.method == 'POST':  # S'il s'agit d'une requête POST
        #make_circos2genomes_form

        form = DBForm(request.POST) #crossplot_form_class(request.POST)


        if form.is_valid():  # Nous vérifions que les données envoyées sont valides
            server, db = manipulate_biosqldb.load_db("Chlamydia_11_14")

            update_genomes_db(server, "Chlamydia_11_14")

            print "form", form
            print dir(form)
            print form.fields
            print "data1", form.cleaned_data


            reference_genome = str(form.cleaned_data['ref_genome'])
            query_genome = str(form.cleaned_data['query_genome'])
            protein_locus = form.cleaned_data['accession']
            region_size = form.cleaned_data['region_size']

            print "asdfasfffffffffffff", reference_genome


            #description2accession = manipulate_biosqldb.description2accession(server, "saureus1")

            #reference_accession = description2accession[reference_genome]
            #query_accession = description2accession[query_genome]

            orthogroup = manipulate_biosqldb.locus_tag2orthogroup_id(server, protein_locus, "Chlamydia_11_14")
            print "orthogroup", orthogroup
            ortho_detail = list(manipulate_biosqldb.orthogroup_id2locus_tag_list(server, orthogroup, "Chlamydia_11_14"))
            locus_tag_list = []
            for i in range(0, len(ortho_detail)):
                print "ortho detail", ortho_detail[i]
                if reference_genome in ortho_detail[i] or query_genome in ortho_detail[i]:
                    locus_tag_list.append(ortho_detail[i][2])
            print "locus_lis", locus_tag_list


            if "Chlamydia_11_14" not in bioentry_in_memory.keys():
                bioentry_in_memory["Chlamydia_11_14"] = {}

            home_dir = os.path.dirname(__file__)
            print "home_dir", home_dir
            temp_location = os.path.join(home_dir, "../assets")
            print "temp loc", temp_location
            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
            print "temp file", temp_file.name
            name = os.path.basename(temp_file.name)
            bioentry_dict = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, "Chlamydia_11_14", locus_tag_list,
                                                                                  temp_file.name, int(region_size),
                                                                                  bioentry_in_memory["Chlamydia_11_14"])




            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = DBForm() #crossplot_form_class()  # Nous créons un formulaire vide
        form2 = DBForm()
    return render(request, 'chlamdb/crossplot.html', locals())




def string_page(request, biodb, cog_id, genome_accession):

    import manipulate_biosqldb
    import urllib2
    server, db = manipulate_biosqldb.load_db(biodb)

    import pandas as pd
    try:
        connect = True
        response = urllib2.urlopen("http://string-db.org/api/psi-mi-tab/interactions?identifier=%s&required_score=600&targetmode=cogs" % cog_id)
        all_cogs = []
        string_interactions = []
        for line in response:
            data = line.rstrip().split('\t')
            string_interactions.append([data[2], data[3], data[-1].split('|')[0]])
            if data[2] not in all_cogs:
                all_cogs.append(data[2])
            if data[3] not in all_cogs:
                all_cogs.append(data[3])
        print all_cogs

        cogs_in_chlamdb = []
        cogs_in_reference = []
        cog2description = {}
        for cog in all_cogs:

            sql1 = 'select functon, name from COG.cog_names_2014 where COG_id="%s"' % cog
            try:
                data = list(server.adaptor.execute_and_fetchall(sql1,)[0])
                cog2description[cog] = "%s (%s)" % (data[1], data[0])
                print data
            except:
                print server.adaptor.execute_and_fetchall(sql1,), cog
                cog2description[cog] = "-"


            try:
                sql = 'select * from COG.locus_tag2gi_hit_chlamydia_12_15 where COG_id="%s" limit 1;' % cog
                sql2 = 'select * from COG.locus_tag2gi_hit_chlamydia_12_15 where COG_id="%s" and accession="%s" limit 1;' % (cog,
                                                                                                                            genome_accession)
                print sql2
                print "############################"
                data = server.adaptor.execute_and_fetchall(sql)
                data2 = server.adaptor.execute_and_fetchall(sql2)
                print "############data", data
                if len(data)>0:
                    cogs_in_chlamdb.append(cog)
                if len(data2)>0:
                    cogs_in_reference.append(cog)
            except:
                print '%s not present in %s' % (cog, biodb)
        for i, data in enumerate(string_interactions):
            string_interactions[i] = data + [cog2description[data[0]], cog2description[data[1]]]

        print string_interactions
        cog_url = '?'
        for i in cogs_in_chlamdb:
            cog_url+= 'cog_list=%s&' % i
        cog_url = cog_url[0:-1]
    except urllib2.URLError:
        connect = False

    print "cogs_in_chlamdb", cogs_in_chlamdb
    return render(request, 'chlamdb/string.html', locals())


def multiple_COGs_heatmap(request, biodb):

    from ete2 import Tree, TextFace
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    cog_list = request.GET.getlist('cog_list')

    cog_filter = '"' + '","'.join(cog_list) + '"'

    cog_annotation_sql = 'select * from COG.cog_names_2014 where COG_id in (%s)' % cog_filter

    cog_annotation = list(server.adaptor.execute_and_fetchall(cog_annotation_sql,))

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
    print tree
    t1 = Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


    sql = 'show columns from comparative_tables.COG_%s' % biodb
    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    #print 'taxons!', ordered_taxons

    ortho_sql = '"' + '","'.join(cog_list) + '"'

    sql = 'select * from comparative_tables.COG_%s where id in (%s)' % (biodb, ortho_sql)

    profile_tuples = list(server.adaptor.execute_and_fetchall(sql,))

    taxon2group2n_homologs = {}

    for i, tuple in enumerate(profile_tuples):
        # get position of the group based on score
        # get colum of taxon i
        taxon2group2n_homologs[tuple[0]] = {}
        for i, taxon in enumerate(ordered_taxons):
            taxon2group2n_homologs[tuple[0]][taxon] = tuple[i+1]
    head = True
    for lf in t1.iter_leaves():
        lf.branch_vertical_margin = 0
        for col, value in enumerate(cog_list):
            if head:

                    'first row, print gene names'
                    #print 'ok!'
                    n = TextFace(' %s ' % str(value))
                    n.rotation= 270
                    n.margin_top = 4
                    n.margin_right = 4
                    n.margin_left = 4
                    n.margin_bottom = 4

                    n.inner_background.color = "white"
                    n.opacity = 1.
                    lf.add_face(n, col, position="aligned")

            n = TextFace(' %s ' % str(taxon2group2n_homologs[value][lf.name]))
            n.margin_top = 4
            n.margin_right = 4
            n.margin_left = 4
            n.margin_bottom = 4
            if taxon2group2n_homologs[value][lf.name] >0:
                n.inner_background.color = "#58ACFA"
            else:
                n.inner_background.color = 'white'

            lf.add_face(n, col, position="aligned")
        lf.name = taxon_id2organism_name[lf.name]
        head=False

    if len(cog_list) > 30:
        big = True
        path = settings.BASE_DIR + '/assets/temp/cog_tree.png'
        asset_path = '/assets/temp/cog_tree.png'
        t1.render(path, dpi=1200, h=600)
    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/cog_tree.svg'
        asset_path = '/assets/temp/cog_tree.svg'

        t1.render(path, dpi=800, h=600)

    return render(request, 'chlamdb/cog_tree.html', locals())

def pfam_tree(request, biodb, orthogroup):
    import ete_motifs
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_locus2protein_id = 'select locus_tag, protein_id from orthology_detail_%s where orthogroup="%s"' % (biodb, orthogroup)

    locus2protein_id= manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_locus2protein_id,))

    locus2pfam_data = ete_motifs.get_pfam_data(orthogroup, biodb)

    motif_count = {}
    for data in locus2pfam_data.values():
        for motif in data:
            print data
            try:
                if motif[4] not in motif_count:
                    motif_count[motif[4]] = [1, motif[5]]
                else:
                    motif_count[motif[4]][0]+=1
            except:
                print "motif", motif

    print "motif_count", motif_count

    sql_tree = 'select phylogeny from biosqldb_phylogenies.%s where orthogroup="%s"' % (biodb, orthogroup)
    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        no_tree = True
        return render(request, 'chlamdb/pfam_tree.html', locals())
    import manipulate_biosqldb

    sql = 'select taxon_id, family from genomes_classification;'

    taxon_id2family = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    print 'tree', tree
    print
    t, ts, leaf_number = ete_motifs.draw_pfam_tree(tree, locus2pfam_data, locus2protein_id, taxon_id2family=False)
    path = settings.BASE_DIR + '/assets/temp/pfam_tree.svg'
    asset_path = '/assets/temp/pfam_tree.svg'
    #print "path", path
    t.render(path, h=leaf_number*12, dpi=800, tree_style=ts)

    return render(request, 'chlamdb/pfam_tree.html', locals())

def TM_tree(request, biodb, orthogroup):
    print 'bonjour', request.method
    import manipulate_biosqldb
    import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)

    locus2TM_data = ete_motifs.get_TM_data(orthogroup, biodb)

    sql_tree = 'select phylogeny from biosqldb_phylogenies.%s where orthogroup="%s"' % (biodb, orthogroup)
    tree = server.adaptor.execute_and_fetchall(sql_tree,)[0]
    print 'tree', tree
    t, ts, leaf_number = ete_motifs.draw_TM_tree(tree, locus2TM_data)
    path = settings.BASE_DIR + '/assets/temp/tm_tree.svg'
    #print "path", path
    t.render(path, h=leaf_number*12, dpi=800, tree_style=ts)

    return render(request, 'chlamdb/pfam_tree.html', locals())


def multiple_orthogroup_heatmap(request, biodb, reference_orthogroup, max_distance=2.2):

    import manipulate_biosqldb
    import biosql_own_sql_tables
    import pandas
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl
    from ete2 import Tree, TextFace

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
    print tree
    t1 = Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


    sql = 'select * from comparative_tables.phylogenetic_profiles_euclidian_distance_%s' \
          ' where (group_1="%s" or group_2="%s") and euclidian_dist <=%s limit 40;' % (biodb,
                                                                          reference_orthogroup,
                                                                          reference_orthogroup,
                                                                            max_distance)
    data = list(server.adaptor.execute_and_fetchall(sql,))

    data_frame = pandas.DataFrame(data)
    sorted_data_frame = data_frame.sort(2)
    print "sorted data", sorted_data_frame

    ordered_orthogroups = []
    for i in sorted_data_frame.itertuples(index=False):
        if i[0] != reference_orthogroup:
            ordered_orthogroups.append(i[0])
        else:
            ordered_orthogroups.append(i[1])

    l = sorted(set(sorted_data_frame[sorted_data_frame.columns[2]]))
    colmap = dict(zip(l,range(len(l))[::-1]))
    norm = mpl.colors.Normalize(vmin=-2, vmax=len(l))
    cmap = cm.OrRd
    cmap_blue = cm.Blues
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)

    orthogroup2distance = {}
    distances = []

    for one_pair in sorted_data_frame.itertuples(index=False):
        distances.append(one_pair[2])
        if one_pair[0] == one_pair[1]:
            orthogroup2distance[one_pair[0]] = one_pair[2]
        elif one_pair[0] == reference_orthogroup:
            orthogroup2distance[one_pair[1]] = one_pair[2]
        elif one_pair[1] == reference_orthogroup:
            orthogroup2distance[one_pair[0]] = one_pair[2]
        else:
            raise 'Error: unexpected combination of groups'
    ordered_distances = sorted(distances)

    sql = 'show columns from orthology_%s' % biodb
    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    ortho_sql = '"' + '","'.join(orthogroup2distance.keys()) + '"' + ',"%s"' % reference_orthogroup

    sql = 'select * from orthology_%s where orthogroup in (%s)' % (biodb, ortho_sql)

    profile_tuples = list(server.adaptor.execute_and_fetchall(sql,))

    taxon2group2n_homologs = {}

    for i, tuple in enumerate(profile_tuples):
        # get position of the group based on score
        # get colum of taxon i
        taxon2group2n_homologs[tuple[0]] = {}
        for i, taxon in enumerate(ordered_taxons):
            taxon2group2n_homologs[tuple[0]][taxon] = tuple[i+1]

    head = True
    for lf in t1.iter_leaves():

        lf.branch_vertical_margin = 0

        for col, value in enumerate(ordered_orthogroups):
            #print 'value', value
            if head:

                    'first row, print gene names'
                    #print 'ok!'
                    n = TextFace(' %s ' % str(value))
                    n.rotation= 270
                    n.margin_top = 4
                    n.margin_right = 4
                    n.margin_left = 4
                    n.margin_bottom = 4
                    if value == reference_orthogroup:
                        n.inner_background.color = "red"
                    else:
                        n.inner_background.color = "white"
                    n.opacity = 1.
                    lf.add_face(n, col, position="aligned")

            n = TextFace(' %s ' % str(taxon2group2n_homologs[value][lf.name]))
            n.margin_top = 4
            n.margin_right = 4
            n.margin_left = 4
            n.margin_bottom = 4
            if taxon2group2n_homologs[value][lf.name] >0:
                if value == reference_orthogroup:
                    n.inner_background.color = "red"
                else:
                    n.inner_background.color = rgb2hex(m2.to_rgba(float(colmap[orthogroup2distance[value]])))
            else:
                n.inner_background.color = 'white'

            lf.add_face(n, col, position="aligned")
        lf.name = taxon_id2organism_name[lf.name]
        head=False


    if len(ordered_orthogroups) > 30:
        big = True
        path = settings.BASE_DIR + '/assets/temp/profile_tree_%s.png' % reference_orthogroup
        asset_path = '/assets/temp/profile_tree_%s.png' % reference_orthogroup
        t1.render(path, dpi=1200, h=800)
    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/profile_tree_%s.svg' % reference_orthogroup
        asset_path = '/assets/temp/profile_tree_%s.svg' % reference_orthogroup
        t1.render(path, dpi=800, h=600)

    # get data about orthogroups

    match_groups_data, raw_data = biosql_own_sql_tables.orthogroup_list2detailed_annotation(ordered_orthogroups, biodb)


    return render(request, 'chlamdb/profile_tree.html', locals())




def orthogroup_conservation_tree(request, biodb, orthogroup):

    import manipulate_biosqldb
    import shell_command

    server, db = manipulate_biosqldb.load_db(biodb)

    asset_path = '/assets/temp/phylo.svg'
    path = settings.BASE_DIR + asset_path
    a,b,c = shell_command.shell_command("rm %s" % path)
    print a, b, c
    import ete_heatmap_conservation
    sql_grp = 'select taxon_id,count(*) from  orthology_detail_%s where orthogroup="%s" group by organism;' % (biodb, orthogroup)
    print sql_grp
    taxid2n = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_grp,))
    tree_sql = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where t2.name="%s"' % biodb
    tree = server.adaptor.execute_and_fetchall(tree_sql,)[0][0]

    taxon_profile = []
    for i in taxid2n:
        if taxid2n[i]>0:
            taxon_profile.append(i)

    url_pattern = ''
    for i in taxon_profile:
        url_pattern += 'taxons_profile=%s&' % i
    url_pattern=url_pattern[0:-1]
    print 'pattern ############################', url_pattern
    print settings.BASE_DIR
    print path
    print taxid2n
    t1, leaf_number = ete_heatmap_conservation.plot_heat_tree(biodb, taxid2n, tree)

    shell_command.shell_command('rm %s' % path)
    print path
    t1.render(path, dpi=800, h=leaf_number*12)


    return render(request, 'chlamdb/orthogroup_conservation.html', locals())



@login_required
def priam_kegg(request, biodb):

    priam_form_class = make_priam_form(biodb)

    print 'request', request.method

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        print 'request', request.method
        form = priam_form_class(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        print 'aaa'
        if form.is_valid():
            genome = form.cleaned_data['genome']



            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = priam_form_class()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/priam_kegg.html', locals())



@login_required
def kegg_module(request, biodb):
    import ete_motifs
    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)
    module_overview_form = make_module_overview_form(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        form = module_overview_form(request.POST)
        if form.is_valid():
            #if request.method == 'POST':  # S'il s'agit d'une requête POST

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            category = form.cleaned_data['category']

            sql_pathway_count = 'select BB.module_name,count_all,count_db,count_db/count_all from (select module_id, count(*) ' \
                                ' as count_db from (select distinct ko_id from enzyme.locus2ko_%s) as t1' \
                                ' inner join enzyme.module2ko as t2 on t1.ko_id=t2.ko_id group by module_id) AA ' \
                                ' right join (select t1.module_id,module_name, count_all from (select module_id, count(*) as count_all ' \
                                'from enzyme.module2ko group by module_id) t1 inner join enzyme.kegg_module as t2 ' \
                                'on t1.module_id=t2.module_id where module_sub_cat="%s")BB on AA.module_id=BB.module_id;' % (biodb, category) # where pathway_category!="1.0 Global and overview maps"

            map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))
            print 'map2count', map2count

            # C.pathway_category,taxon_id, A.pathway_name,A.n_enzymes, C.description
            sql = 'select B.module_sub_cat,A.taxon_id,B.module_name,A.n,B.description from ' \
                                ' (select taxon_id, module_id, count(*) as n from ' \
                                ' (select distinct taxon_id,ko_id from enzyme.locus2ko_%s) t1 ' \
                                ' inner join enzyme.module2ko as t2 on t1.ko_id=t2.ko_id group by taxon_id, module_id) A ' \
                                ' inner join enzyme.kegg_module as B on A.module_id=B.module_id where module_sub_cat="%s";' % (biodb, category)

            print sql
            pathway_data = server.adaptor.execute_and_fetchall(sql,)
            all_maps = []
            category2maps = {}
            # pathway cat 2 taxon_id 2 pathway_map 2 [count, pathway description]
            pathway_category2taxon2map = {}
            for one_row in pathway_data:
                # first pathway category
                if one_row[0] not in pathway_category2taxon2map:
                    category2maps[one_row[0]] = [[one_row[2],one_row[4]]]
                    all_maps.append(one_row[2])
                    pathway_category2taxon2map[one_row[0]] = {}
                    pathway_category2taxon2map[one_row[0]][one_row[1]] = {}
                    pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
                else:

                    if one_row[2] not in all_maps:
                        category2maps[one_row[0]].append([one_row[2],one_row[4]])
                        all_maps.append(one_row[2])
                    # if noew taxon
                    if one_row[1] not in pathway_category2taxon2map[one_row[0]]:
                        pathway_category2taxon2map[one_row[0]][one_row[1]] = {}

                        pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
                    # if new map for existing taxon
                    else:

                        pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]

            tree = ete_motifs.pathways_heatmap(biodb,
                                              category2maps,
                                              pathway_category2taxon2map,map2count)


            #except:
            #    tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, merged_dico)

            if len(all_maps) > 1000:
                big = True
                path = settings.BASE_DIR + '/assets/temp/metabo_tree.png'
                asset_path = '/assets/temp/metabo_tree.png'
                tree.render(path, dpi=1200, h=600)
            else:
                big = False
                path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
                asset_path = '/assets/temp/metabo_tree.svg'

                tree.render(path, dpi=800, h=600)

            envoi = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = module_overview_form()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/module_overview.html', locals())


@login_required
def module_comparison(request, biodb):

    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        form = comp_metabo_form(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            print 'db id', database_id

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


            sql_pathway_count = 'select BB.module_name,count_all,count_db,count_db/count_all from (select module_id, count(*) ' \
                        ' as count_db from (select distinct ko_id from enzyme.locus2ko_%s) as t1' \
                        ' inner join enzyme.module2ko as t2 on t1.ko_id=t2.ko_id group by module_id) AA ' \
                        ' right join (select t1.module_id,module_name, count_all from (select module_id, count(*) as count_all ' \
                        'from enzyme.module2ko group by module_id) t1 inner join enzyme.kegg_module as t2 ' \
                        'on t1.module_id=t2.module_id)BB on AA.module_id=BB.module_id;' % (biodb)

            map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))
            print 'map2count', map2count
            category2maps = {}

            sql_category2maps = 'select module_sub_cat,module_name,description from enzyme.kegg_module;'

            data = server.adaptor.execute_and_fetchall(sql_category2maps,)

            for one_map in data:
                if one_map[0] not in category2maps:
                    category2maps[one_map[0]] = [[one_map[1], one_map[2]]]
                else:
                    category2maps[one_map[0]].append([one_map[1], one_map[2]])

            print "category2maps", category2maps

            sql = 'select distinct module_name,module_sub_sub_cat from enzyme.kegg_module;'
            module2sub_category = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            taxon_maps = []
            for taxon in taxon_list:
                database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

                sql = 'select module_name, n from (select B.module_id,count(*) as n from ' \
                      ' (select * from enzyme.locus2ko_%s where taxon_id=%s) A ' \
                      ' left join enzyme.module2ko as B on A.ko_id=B.ko_id group by module_id) AA ' \
                      ' right join enzyme.kegg_module as BB on AA.module_id=BB.module_id;' % (biodb, taxon)

                print sql
                map2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                taxon_maps.append(map2count_taxon)


            envoi_comp = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = comp_metabo_form()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/module_comp.html', locals())



@login_required
def metabo_overview(request, biodb):
    import ete_motifs
    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)

    #if request.method == 'POST':  # S'il s'agit d'une requête POST

    sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

    print 'db id', database_id

    sql_pathway_count = 'select PATH2.pathway_name,PATH2.n,PATH1.n,PATH1.n/PATH2.n from (select pathway_name,count(*) as n ' \
                        ' from (select distinct pathway_name,t4.ec_id from enzyme.locus2ec_%s as t3 ' \
                        ' inner join enzyme.kegg2ec as t4 on t3.ec_id=t4.ec_id inner join enzyme.kegg_pathway as t5 ' \
                        ' on t4.pathway_id=t5.pathway_id where pathway_category!="1.0 Global and overview maps") A ' \
                        ' group by pathway_name) ' \
                        ' PATH1 right join ' \
                        ' (select pathway_name,count(*) as n from enzyme.kegg2ec as t1 ' \
                        ' inner join enzyme.kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                        '  group by pathway_name) ' \
                        ' PATH2 on PATH2.pathway_name=PATH1.pathway_name;' % biodb # where pathway_category!="1.0 Global and overview maps"

    map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))
    print 'map2count', map2count

    sql = 'select C.pathway_category,taxon_id, A.pathway_name,A.n_enzymes, C.description from ' \
          '( select distinct taxon_id,pathway_name, count(*) as n_enzymes from (' \
          'select distinct taxon_id, ec_id  from enzyme.locus2ec_%s as b1 ' \
          'inner join biosqldb.bioentry as b2 on b1.accession=b2.accession where biodatabase_id=%s) ' \
          't1 inner join enzyme.kegg2ec as t2  on t1.ec_id=t2.ec_id ' \
          'inner join enzyme.kegg_pathway as t3 on t2.pathway_id=t3.pathway_id ' \
          'group by taxon_id,pathway_name) A ' \
          'inner join enzyme.kegg_pathway as C  on A.pathway_name=C.pathway_name;' % (biodb, database_id)

    print sql
    pathway_data = server.adaptor.execute_and_fetchall(sql,)
    all_maps = []
    category2maps = {}
    # pathway cat 2 taxon_id 2 pathway_map 2 [count, pathway description]
    pathway_category2taxon2map = {}
    for one_row in pathway_data:
        # first pathway category
        if one_row[0] not in pathway_category2taxon2map:
            category2maps[one_row[0]] = [[one_row[2],one_row[4]]]
            all_maps.append(one_row[2])
            pathway_category2taxon2map[one_row[0]] = {}
            pathway_category2taxon2map[one_row[0]][one_row[1]] = {}
            pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
        else:

            if one_row[2] not in all_maps:
                category2maps[one_row[0]].append([one_row[2],one_row[4]])
                all_maps.append(one_row[2])
            # if noew taxon
            if one_row[1] not in pathway_category2taxon2map[one_row[0]]:
                pathway_category2taxon2map[one_row[0]][one_row[1]] = {}

                pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
            # if new map for existing taxon
            else:

                pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]

    tree = ete_motifs.pathways_heatmap(biodb,
                                      category2maps,
                                      pathway_category2taxon2map,map2count)


    #except:
    #    tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, merged_dico)

    if len(all_maps) > 1000:
        big = True
        path = settings.BASE_DIR + '/assets/temp/metabo_tree.png'
        asset_path = '/assets/temp/metabo_tree.png'
        tree.render(path, dpi=1200, h=600)
    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
        asset_path = '/assets/temp/metabo_tree.svg'

        tree.render(path, dpi=800, h=600)

    envoi = True

    #else:  # Si ce n'est pas du POST, c'est probablement une requête GET
    #    pass  # Nous créons un formulaire vide

    return render(request, 'chlamdb/metabo_overview.html', locals())


@login_required
def metabo_comparison(request, biodb):

    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        form = comp_metabo_form(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            print 'db id', database_id

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


            sql_pathway_count = 'select PATH2.pathway_name,PATH2.n,PATH1.n,PATH1.n/PATH2.n from (select pathway_name,count(*) as n ' \
                                ' from (select distinct pathway_name,t4.ec_id from enzyme.locus2ec_%s as t3 ' \
                                ' inner join enzyme.kegg2ec as t4 on t3.ec_id=t4.ec_id inner join enzyme.kegg_pathway as t5 ' \
                                ' on t4.pathway_id=t5.pathway_id where pathway_category!="1.0 Global and overview maps") A ' \
                                ' group by pathway_name) ' \
                                ' PATH1 right join ' \
                                ' (select pathway_name,count(*) as n from enzyme.kegg2ec as t1 ' \
                                ' inner join enzyme.kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                                ' where pathway_category!="1.0 Global and overview maps" group by pathway_name) ' \
                                ' PATH2 on PATH2.pathway_name=PATH1.pathway_name;' % biodb

            map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))
            print 'map2count', map2count
            category2maps = {}

            sql_category2maps = 'select pathway_category,pathway_name,description from  enzyme.kegg2ec as t1 ' \
                                ' inner join enzyme.kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                                ' where pathway_category!="1.0 Global and overview maps" group by pathway_name;'

            data = server.adaptor.execute_and_fetchall(sql_category2maps,)

            for one_map in data:
                if one_map[0] not in category2maps:
                    category2maps[one_map[0]] = [[one_map[1], one_map[2]]]
                else:
                    category2maps[one_map[0]].append([one_map[1], one_map[2]])

            print "category2maps", category2maps

            taxon_maps = []
            for taxon in taxon_list:

                sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

                biodatabase_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

                sql = 'select PATH2.pathway_name,PATH1.n from (select pathway_name,count(*) as n from ' \
                      ' (select distinct pathway_name,t4.ec_id from enzyme.locus2ec_%s as t3 ' \
                      ' inner join enzyme.kegg2ec as t4 on t3.ec_id=t4.ec_id inner join enzyme.kegg_pathway as t5 ' \
                      ' on t4.pathway_id=t5.pathway_id inner join biosqldb.bioentry as t6 on t3.accession=t6.accession ' \
                      ' where biodatabase_id=%s and pathway_category!="1.0 Global and overview maps" and t6.taxon_id=%s) A ' \
                      ' group by pathway_name) PATH1 right join (select pathway_name,count(*) as n ' \
                      ' from enzyme.kegg2ec as t1 inner join enzyme.kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                      ' where pathway_category!="1.0 Global and overview maps" ' \
                      ' group by pathway_name) PATH2 on PATH2.pathway_name=PATH1.pathway_name order by n;' % (biodb, biodatabase_id, taxon)
                map2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                taxon_maps.append(map2count_taxon)


            envoi_comp = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = comp_metabo_form()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/metabo_comp.html', locals())

@login_required
def pfam_comparison(request, biodb):

    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        form = comp_metabo_form(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            print 'db id', database_id

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            columns = '`' + '`,`'.join(taxon_list) + '`'
            filter = '(`' + '`>0 or`'.join(taxon_list) + '`>0)'


            sql = 'select * from (select id from comparative_tables.Pfam_%s where %s group by id) A' \
                  ' inner join (select distinct signature_accession,signature_description,count(*) as n ' \
                  ' from interpro_%s where analysis="Pfam" group by signature_accession) B on A.id = B.signature_accession' % (biodb,filter, biodb)

            sql_pathway_count = 'select distinct signature_accession,signature_description,count(*) as n ' \
                                ' from interpro_%s where analysis="Pfam" group by signature_accession;' % biodb

            pfam_data_raw = server.adaptor.execute_and_fetchall(sql,)
            pfam2data = {}
            for one_pfam_entry in pfam_data_raw:
                pfam2data[one_pfam_entry[0]] = one_pfam_entry[1:]




            taxon_dicos = []
            for taxon in taxon_list:

                sql = 'select distinct signature_accession,count(*) as n ' \
                      ' from interpro_%s where analysis="Pfam" and taxon_id=%s ' \
                      ' group by signature_accession;' % (biodb,taxon)

                accession2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


                for accession in pfam2data:
                    if accession not in accession2count_taxon:
                        accession2count_taxon[accession] = 0
                taxon_dicos.append(accession2count_taxon)

            envoi_comp = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = comp_metabo_form()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/pfam_comp.html', locals())


@login_required
def orthogroup_comparison(request, biodb):

    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        form = comp_metabo_form(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            print 'db id', database_id

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            columns = '`' + '`,`'.join(taxon_list) + '`'
            filter = '(`' + '`>1 or`'.join(taxon_list) + '`>1)'


            sql = 'select orthogroup,count(*) from orthology_detail_%s group by orthogroup' % (biodb)
            print sql

            orthogroups2total_count= manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select orthogroup,%s from comparative_tables.orthology_%s where %s' % (columns, biodb, filter)
            print sql

            orthogroups2counts = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select orthogroup,product, count(*) from orthology_detail_%s group by orthogroup,product;' % biodb

            group2annot = {}
            for i in server.adaptor.execute_and_fetchall(sql,):
                #print i
                if i[0] not in group2annot:
                    group2annot[i[0]] = [i[1:]]
                else:
                    if len(group2annot[i[0]])<5:
                        group2annot[i[0]].append(i[1:])
                    else:
                        if ['...', '%s homologs' % orthogroups2total_count[i[0]]] not in group2annot[i[0]]:
                            group2annot[i[0]].append(['...', '%s homologs' % orthogroups2total_count[i[0]]])
                        else:
                            continue
            envoi_comp = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = comp_metabo_form()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/ortho_comp.html', locals())

@login_required
def ko_comparison(request, biodb):

    print 'request', request.method
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        form = comp_metabo_form(request.POST)  # Nous reprenons les données
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            print 'db id', database_id

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            columns = '`' + '`,`'.join(taxon_list) + '`'
            filter = '(`' + '`>0 or`'.join(taxon_list) + '`>0)'


            sql = 'select ko_id, count(*) from enzyme.locus2ko_%s group by ko_id;' % (biodb)
            print sql

            ko2total_count= manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select id,%s from comparative_tables.ko_%s where %s' % (columns, biodb, filter)
            print sql

            ko2counts = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select ko_id,definition from enzyme.ko_annotation'

            ko2annot = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            envoi_comp = True

    else:  # Si ce n'est pas du POST, c'est probablement une requête GET
        form = comp_metabo_form()  # Nous créons un formulaire vide

    return render(request, 'chlamdb/ko_comp.html', locals())