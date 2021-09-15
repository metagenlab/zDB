# -*- coding: utf-8 -*-

# todo circos gc file curently written in home directory, move it to other place
# todo save temp files in temp folder




#from django.shortcuts import render
#from datetime import datetime
#from bin.metagenlab_libs.metagenlab_libs.db_utils import DB
from django.shortcuts import render
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import chlamdb.plots
import seaborn as sns
import numpy 
import re
# from django.core.cache import cache
#import pylibmc
#from django.core.cache import cache
import os
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from chlamdb.forms import make_contact_form
from chlamdb.forms import make_plot_form
from chlamdb.forms import SearchForm
from chlamdb.forms import BiodatabaseForm
from chlamdb.forms import make_circos_form
from chlamdb.forms import make_circos2genomes_form
from chlamdb.forms import make_mummer_form
from chlamdb.forms import make_blast_form
from chlamdb.forms import make_crossplot_form
from chlamdb.forms import ConnexionForm
from chlamdb.forms import DBForm
from chlamdb.forms import make_motif_form
from chlamdb.forms import PCRForm
from chlamdb.forms import make_extract_form
from chlamdb.forms import make_circos_orthology_form
from chlamdb.forms import make_interpro_from
from chlamdb.forms import make_metabo_from
from chlamdb.forms import make_module_overview_form
from chlamdb.forms import make_extract_region_form
from chlamdb.forms import make_venn_from
from chlamdb.forms import make_priam_form
from chlamdb.forms import AnnotForm
from chlamdb.forms import hmm_sets_form
from chlamdb.forms import hmm_sets_form_circos
from chlamdb.forms import make_blastnr_form
from chlamdb.forms import make_genome_selection_form
from chlamdb.forms import make_comment_from
from chlamdb.forms import locus_int_form
from chlamdb.forms import LocusInt
from chlamdb.forms import make_species_curation_form
from chlamdb.forms import get_LocusAnnotForm
from chlamdb.forms import make_pathway_overview_form
from chlamdb.forms import make_interpro_taxonomy
from chlamdb.forms import BlastProfileForm
from chlamdb.forms import make_pairwiseid_form
from chlamdb.forms import make_locus2network_form
from chlamdb.forms import heatmap_form
from chlamdb.forms import blast_sets_form
from chlamdb.forms import make_kegg_form
from chlamdb.forms import transporters_superfam_form
from chlamdb.forms import make_pairwiseCDS_length_form
from django.contrib.auth import logout
from django.conf import settings
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.forms.utils import flatatt
from chlamdb.forms import make_blastnr_best_non_top_phylum_form



from chlamdb.biosqldb import manipulate_biosqldb
from chlamdb.biosqldb import mysqldb_plot_genomic_feature
#from django.core.cache import caches
from django.core.cache import cache
from tempfile import NamedTemporaryFile
from Bio import SeqIO
from chlamdb.biosqldb.gbk2table import Record
import chlamdb.models
import simplejson
import string
import random
import json
from django.shortcuts import render
from celery.result import AsyncResult
from django.http import HttpResponse
from django.http import StreamingHttpResponse
from chlamdb.forms import GenerateRandomUserForm
from chlamdb.tasks import run_circos
from chlamdb.tasks import run_circos_main
from chlamdb.tasks import extract_interpro_task
from chlamdb.tasks import plot_neighborhood_task
from chlamdb.tasks import TM_tree_task
from chlamdb.tasks import pfam_tree_task
from chlamdb.tasks import phylogeny_task
from chlamdb.tasks import plot_heatmap_task
from chlamdb.tasks import KEGG_map_ko_task
from chlamdb.tasks import KEGG_map_ko_organism_task
from chlamdb.tasks import basic_tree_task
from chlamdb.celeryapp import app as celery_app

from metagenlab_libs import db_utils
from metagenlab_libs.ete_phylo import EteTree, SimpleColorColumn, ModuleCompletenessColumn
from metagenlab_libs.ete_phylo import KOAndCompleteness
from metagenlab_libs.ete_phylo import Column
from metagenlab_libs.KO_module import ModuleParser
from metagenlab_libs.chlamdb import search_bar as sb

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Seq import reverse_complement, translate
from tempfile import NamedTemporaryFile
from io import StringIO
from Bio.Blast import NCBIXML
          
            
import os
from chlamdb.biosqldb import shell_command
import re

from ete3 import Tree
from reportlab.lib import colors
from ete3 import TextFace, StackedBarFace
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation


with db_utils.DB.load_db_from_name(settings.BIODB_DB_PATH) as db:
    hsh_config = db.get_config_table(ret_mandatory=True)
    optional2status = { name: value for name, (mandatory, value) in hsh_config.items() if not mandatory}
    missing_mandatory = [name for name, (mandatory, value) in hsh_config.items()
            if mandatory and not value]
    

def my_locals(local_dico):
    local_dico["optional2status"] = optional2status
    local_dico["missing_mandatory"] = missing_mandatory
    return local_dico

@celery_app.task(bind=True)
def debug_task(self):
    print('Request: {0!r}'.format(self.request))


def get_task_info(request):
    
    task_id = request.GET.get('task_id', None)
    
    if task_id is not None:
        task = AsyncResult(task_id)
        data = {
            'state': task.state,
            'result': task.result,
        }
        
        return HttpResponse(json.dumps(data), content_type='application/json')
    else:
        return HttpResponse('No job id given.')


if settings.DEBUG == True:
    debug_mode = 1
else:
    debug_mode = 0





def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))


def extract_alphanumeric(input_string):
    from string import ascii_letters, digits
    import string
    return "".join([ch for ch in input_string if ch in (ascii_letters + digits + '_-.')])


def choose_db(request):
    
    server, db = manipulate_biosqldb.load_db(biodb)
    if request.method == 'POST': 

        form = BiodatabaseForm(request.POST)

        if form.is_valid():  

            biodb = form.cleaned_data['biodatabase']
            envoi = True
    else:  
        form = BiodatabaseForm()

    return render(request, 'chlamdb/choose_db.html', my_locals(locals()))

def help(request):

    return render(request, 'chlamdb/help.html', my_locals(locals()))


def about(request):
    import re
    import bibtexparser
    path = settings.BASE_DIR + '/assets/bibliography/references.bib'
    with open(path) as bibtex_file:
        bib_database = bibtexparser.load(bibtex_file)

    entry_list = []

    for entry in bib_database.entries:

        """
        Sophie S Abby, Jean Cury, Julien Guglielmini, Bertrand Néron, Marie Touchon, and Eduardo PC
        Rocha. Identification of protein secretion systems in bacterial genomes. Scientific reports, 6, 2016.
        53, 57, 61, 164

        Peter JA Cock, Tiago  Antao, Jeffrey T Chang, Brad
        Biopython: freely    available    python    tools    for computational molecular biology and bioinformatics.Bioinformatics, 25(11):1422–1423, 2009.

        {'publisher': 'Oxford University Press', 'year': '2009', 'pages': '1422--1423', 'number': '11', 'volume': '25',
         'journal': 'Bioinformatics',
         'author': 'Cock, Peter JA and Antao, Tiago and Chang, Jeffrey T and Chapman, Brad A and Cox, Cymon J and Dalke, Andrew and Friedberg, Iddo and Hamelryck, Thomas and Kauff, Frank and Wilczynski, Bartek and others',
         'title': 'Biopython: freely available Python tools for computational molecular biology and bioinformatics',
         'ENTRYTYPE': 'article', 'ID': 'cock2009biopython'}

        """
        string = ("<b>%s</b></br> %s, %s, %s(%s):%s, %s" % (re.sub('[{}]','', entry["title"]),
                                                         entry["author"],
                                                            entry["journal"],
                                                            entry["volume"],
                                                            entry["number"],
                                                            entry["pages"],
                                                            entry["year"],
                           ))
        url = entry["url"]
        entry_list.append([string, url])

    return render(request, 'chlamdb/credits.html', my_locals(locals()))

def create_user(username, mail, password, first_name, last_name, staff=True):
    from django.contrib.auth.models import User
    user = User.objects.create_user(username, mail, password)
    user.first_name, user.last_name = first_name, last_name
    user.is_staff = staff
    user.save()



def chlamdb_login(request):
    error = False

    if request.method == "POST":
        form = ConnexionForm(request.POST)
        if form.is_valid():
            username = form.cleaned_data["username"]
            password = form.cleaned_data["password"]
            user = authenticate(username=username, password=password)  # Nous vérifions si les données sont correctes
            if user:  # Si l'objet renvoyé n'est pas None
                login(request, user)  # nous connectons l'utilisateur
                #return HttpResponseRedirect("/chlamdb/home")
            else: # sinon une erreur sera affichée
                error = True

        if form.is_valid():  

            biodb = form.cleaned_data['biodatabase']
            envoi = True
    else:
        form = ConnexionForm()

    return render(request, 'chlamdb/login.html', my_locals(locals()))

def logout_view(request):
    logout(request)
    return render(request, 'chlamdb/logout.html', my_locals(locals()))
    return render(request, 'chlamdb/logout.html', my_locals(locals()))


class StackedBarColumn(Column):
    def __init__(self, values, header, colours=None, relative=False,
            face_params=None, header_params=None):
        super().__init__(header, face_params, header_params)
        self.values=values
        self.colours = colours
        self.relative = relative
        if relative:
            self.max = max(values.values())
            self.min = min(values.values())

    def get_face(self, index):
        val = self.values[int(index)]

        if self.relative and self.max!=self.min:
            val = 100*float(val-self.min)/(self.max-self.min)
        elif self.relative:
            val = 100

        face = StackedBarFace([val, 100-val], width=50, height=9, colors=self.colours)
        self.set_default_params(face)
        face.inner_border.color = "black"
        face.inner_border.width = 0
        return face


def home(request):

    INTRO=settings.INTRO
    TITLE=settings.TITLE
    SUBTITLE=settings.SUBTITLE

    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100*x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x/pow(10,6), 2))

    data_table_header = ["Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot", "n_contigs", "length", "coding_density"]].values.tolist()

    # plot phylo only of not already in assets
    tree = db.get_reference_phylogeny()
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    if not R is None:
        t1.set_outgroup(R)
    t1.ladderize()

    e_tree = EteTree(t1)
    header_params = {"rotation": -30}
    stacked_face_params = {"margin_right": 8, "margin_left": 5}

    tree_params = [
            # serie_name, header, color and is_relative
            ["length", "Size (Mbp)", "#91bfdb", True],
            ["gc", "GC %", "#fc8d59", False],
            ["coding_density", "Coding density %", "#99d594", False],
            ["completeness", "Completeness", "#d7191c", False],
            ["contamination", "Contamination", "black", False]]
        
    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours = [col, "white"],
                relative=is_relative, header=header,
                header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)


    number_of_files=db.count_files()
    print("number_of_files", number_of_files)
     
    orthogroups_freq=db.get_all_orthogroups( min_size=None)
    df_ort=pd.DataFrame(orthogroups_freq, columns=["Orthogroup", "freq"])
    number_ort= df_ort.shape[0]
    print("number_of_orthogroups", number_ort)

    description_db = db.get_genomes_description()
    print("description_db",description_db)
   
    taxids = list(description_db.index)

    df_hits = db.get_og_count(taxids, search_on="taxid")
    missing_entries = df_hits[df_hits == 0].count(axis=1)
    core = len(missing_entries[missing_entries == 0])
    print("core",core)

    return render(request, 'chlamdb/home.html', my_locals(locals()))


def curated_taxonomy(request):
    
    from chlamdb.phylo_tree_display import phylo_tree_bar
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)  
      
    sql = 'select distinct t5.AssemblyAccession,t1.accession,t1.taxon_id as assembly_id,t1.description,t3.* from bioentry t1' \
            ' inner join taxid2species t2 on t1.taxon_id=t2.taxon_id ' \
            ' inner join species_curated_taxonomy t3 on t2.species_id=t3.species_id ' \
            ' left join bioentry2assembly t4 on t1.bioentry_id=t4.bioentry_id ' \
            ' left join assembly_metadata t5 on t4.assembly_id=t5.assembly_id;'
            
    data = server.adaptor.execute_and_fetchall(sql,)

    header2taxon2text = {}
    for n, row in enumerate(data):
        taxon_id = row[2]
        if n == 0:
            header2taxon2text["species_id"] = {}
            header2taxon2text["phylum"] = {}
            header2taxon2text["order"] = {}
            header2taxon2text["family"] = {}
            header2taxon2text["genus"] = {}
            header2taxon2text["species"] = {}

        header2taxon2text["species_id"][taxon_id] = row[4]
        header2taxon2text["phylum"][taxon_id] = row[5]
        header2taxon2text["order"][taxon_id] = row[6]
        header2taxon2text["family"][taxon_id] = row[7]
        header2taxon2text["genus"][taxon_id] = row[8]
        header2taxon2text["species"][taxon_id] = row[9]

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 ' \
               ' on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    tree, style = phylo_tree_bar.plot_tree_text_metadata(tree,
                                                         header2taxon2text,
                                                         ["species_id","phylum", "order", "family", "genus", "species"],
                                                         biodb)

    path1 = settings.BASE_DIR + '/assets/temp/interpro_tree2.svg'
    asset_path1 = '/temp/interpro_tree2.svg'
    tree.render(path1, dpi=550, tree_style=style)

    return render(request, 'chlamdb/curated_taxonomy.html', my_locals(locals()))


def edit_species_taxonomy(request, species_id):
    
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    curation_form_class = make_species_curation_form(biodb, species_id)

    if request.method == 'POST': 

        form = curation_form_class(request.POST)
        if form.is_valid():
            phylum = form.cleaned_data['phylum']
            order = form.cleaned_data['order']
            family = form.cleaned_data['family']
            genus = form.cleaned_data['genus']
            species = form.cleaned_data['species']
            
            sql1 = 'update species_curated_taxonomy set phylum="%s" where species_id=%s' % (phylum, species_id)
            sql2 = 'update species_curated_taxonomy set `order`="%s" where species_id=%s' % (order, species_id)
            sql3 = 'update species_curated_taxonomy set family="%s" where species_id=%s' % (family, species_id)
            sql4 = 'update species_curated_taxonomy set genus="%s" where species_id=%s' % (genus, species_id)
            sql5 = 'update species_curated_taxonomy set species="%s" where species_id=%s' % (species, species_id)
            
            server.adaptor.execute(sql1,)
            server.adaptor.execute(sql2,)
            server.adaptor.execute(sql3,)
            server.adaptor.execute(sql4,)
            server.adaptor.execute(sql5,)
            server.commit()
            
            saved = True

            return curated_taxonomy(request)


    else: 
        form = curation_form_class()  

    return render(request, 'chlamdb/edit_species_taxonomy.html', my_locals(locals()))


def circos_homology(request):

    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    circos_orthology_form_class = make_circos_orthology_form(biodb)

    if request.method == 'POST': 

        form = circos_orthology_form_class(request.POST)

        if form.is_valid():

            accession = form.cleaned_data['accession']

            sql = 'select accession from bioentry' \
                  ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
                  ' and biodatabase.name = "%s"' % biodb



            all_accession = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
            columns = 'orthogroup, accession, start, stop'
            sql_ref = 'select %s from orthology_detail where locus_tag = "%s" or protein_id = "%s" or orthogroup = "%s"' % (columns,
                                                                                                                            accession,
                                                                                                                            accession, 
                                                                                                                            accession)




            ref_record = server.adaptor.execute_and_fetchall(sql_ref,)[0]

            orthogroup = ref_record[0]

            columns = 'accession, start, stop'
            sql_targets = 'select %s from orthology_detail where orthogroup ="%s"' % (columns,
                                                                                      orthogroup)

            target_records = server.adaptor.execute_and_fetchall(sql_targets,)

            record_list = []
            for accession in all_accession:
                if accession == "CP001848" or accession == "BX119912":
                    continue
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


            path = settings.BASE_DIR + "/assets/circos"

            chlamdb.plots.circos_orthology.circos_orthology(record_list, ref_record[1:], target_records, location=path)
            circos_file = "circos/circos_ortho.svg"
            envoi_circos = True

    else:  
        form = circos_orthology_form_class()
    return render(request, 'chlamdb/circos_homology.html', my_locals(locals()))


def format_lst_to_html(lst_elem, add_count=True, format_func=None):
    if format_func==None:
        format_func = lambda x:x

    dict_elem = {}
    for elem in lst_elem:
        if pd.isna(elem):
            elem = "-"
        cnt = dict_elem.get(elem, 0)
        dict_elem[elem] = cnt+1

    elems = []
    for k, v in dict_elem.items():
        if k != "-":
            token = format_func(k)
        else:
            token = k
        if add_count and k!="-":
            elems.append(f"{token} ({v})")
        else:
            elems.append(f"{token}")
    return "<br/>".join(elems)


def get_optional_annotations(db, seqids):
    header = []
    config_table = db.get_config_table()
    annotations = []
    if config_table.get("KEGG", False):
        header.append("KO")
        ko_hits = db.get_ko_hits(seqids, search_on="seqid", indexing="seqid")
        annotations.append(ko_hits)
    if config_table.get("COG", False):
        header.append("COG")
        cog_hits = db.get_cog_hits(seqids, indexing="seqid", search_on="seqid")
        annotations.append(cog_hits)

    if len(annotations)==2:
        return header, annotations[0].join(annotations[1], how="outer")
    elif len(annotations)==1:
        return header, annotations[0]
    return header, pd.DataFrame()


def get_table_details(db, annotations):
    header = ["Orthogroup", "Organism", "Locus", "Gene", "Product"]
    hsh_organisms = db.get_organism(annotations.index.tolist())
    infos = []
    for seqid, data in annotations.iterrows():
        organism = hsh_organisms[seqid]
        og = format_orthogroup(data.orthogroup, to_url=True)
        gene = data.gene
        if pd.isna(data.gene):
            gene = "-"
        locus = format_locus(data.locus_tag, to_url=True)
        entry = [og, organism, locus, gene, data["product"]]
        infos.append(entry)
    return header, infos


def extract_orthogroup(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    extract_form_class = make_extract_form(db, plasmid=True)
    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        form = extract_form_class()
        # add error message in web page
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))
    
    include_taxids, include_plasmids = form.get_include_choices()
    exclude_taxids, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()
    single_copy = "checkbox_single_copy" in request.POST

    sum_include_lengths = len(include_taxids)
    if not include_plasmids is None:
        sum_include_lengths += len(include_plasmids)

    if n_missing>=sum_include_lengths:
        wrong_n_missing = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    og_counts_in = db.get_og_count(include_taxids, plasmids=include_plasmids)
    if not single_copy:
        og_counts_in["presence"] = og_counts_in[og_counts_in > 0].count(axis=1)
    else:
        og_counts_in["presence"] = og_counts_in[og_counts_in == 1].count(axis=1)

    og_counts_in["selection"] = og_counts_in.presence >= (sum_include_lengths-n_missing)

    sum_exclude_lengths = len(exclude_taxids)
    if not exclude_plasmids is None:
        sum_exclude_lengths += len(exclude_plasmids)
    if sum_exclude_lengths>0:
        mat_exclude = db.get_og_count(exclude_taxids, plasmids=exclude_plasmids)
        mat_exclude["presence"] = mat_exclude[mat_exclude > 0].count(axis=1)
        mat_exclude["exclude"] = mat_exclude.presence > 0
        neg_index = mat_exclude[mat_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    pos_index = og_counts_in[og_counts_in.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        no_match = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    count_all_genomes = db.get_og_count(selection, search_on="orthogroup")

    if not single_copy:
        orthogroup2count_all = count_all_genomes[count_all_genomes > 0].count(axis=1)
    else:
        orthogroup2count_all = count_all_genomes[count_all_genomes == 1].count(axis=1)
    max_n = orthogroup2count_all.max()
    match_groups_data = []

    all_taxids = include_taxids
    if not include_plasmids is None:
        all_taxids += include_plasmids
    annotations = db.get_genes_from_og(orthogroups=selection, taxon_ids=all_taxids,
        terms=["gene", "product", "locus_tag"])
    if annotations.empty:
        no_match = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))
    

    opt_header, optional_annotations = get_optional_annotations(db, seqids=annotations.index.tolist())
    details_header, details_data = get_table_details(db, annotations)
    annotations = annotations.join(optional_annotations)
    grouped = annotations.groupby("orthogroup")
    genes = grouped["gene"].apply(list)
    products = grouped["product"].apply(list)

    if "COG" in opt_header:
        cogs = grouped["cog"].apply(list)

    if "KO" in opt_header:
        kos = grouped["ko"].apply(list)

    table_headers = ["Orthogroup", "Genes", "Products"]
    table_headers.extend(opt_header)
    table_headers.extend([f"Present in {sum_include_lengths}", f"Freq complete database (/{max_n})"])

    for row, count in orthogroup2count_all.iteritems():
        cnt_in = og_counts_in.presence.loc[row]
        g = genes.loc[row]
        gene_data = format_lst_to_html("-" if pd.isna(entry) else entry for entry in g)
        prod_data = format_lst_to_html(products.loc[row])
        column_header = format_orthogroup(row)
        optional = []
        if "KO" in opt_header and row in kos:
            optional.append(format_lst_to_html(kos.loc[row], add_count=True, format_func=format_ko_url))
        if "COG" in opt_header and row in cogs:
            optional.append(format_lst_to_html(cogs.loc[row], add_count=True, format_func=format_cog_url))
        entry = [column_header, gene_data, prod_data, *optional, cnt_in, count]
        match_groups_data.append(entry)

    envoi_extract = True
    return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))


def locus_list2orthogroups(request):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST': 

        form = AnnotForm(request.POST)

        if form.is_valid():  
            from chlamdb.biosqldb import biosql_own_sql_tables
            from chlamdb.phylo_tree_display import ete_motifs
            server, db = manipulate_biosqldb.load_db(biodb)

            match_locus = [i.rstrip() for i in form.cleaned_data['orthogroups'].rstrip().split('\n')]

            sql = 'select locus_tag, orthogroup from orthology_detail'

            locus2group = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            group2count = {}

            for locus in match_locus:
                if locus2group[locus] not in group2count:
                    group2count[locus2group[locus]] = 1
                else:
                    group2count[locus2group[locus]] += 1

            total = len(group2count)

            envoi_annot = True
    else:  

        form = AnnotForm()

    return render(request, 'chlamdb/locus_list2orthogroups.html', my_locals(locals()))



def orthogroup_annotation(request, display_form):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST': 

        form = AnnotForm(request.POST)

        if form.is_valid():  
            from chlamdb.biosqldb import biosql_own_sql_tables
            from chlamdb.phylo_tree_display import ete_motifs
            server, db = manipulate_biosqldb.load_db(biodb)

            match_groups = [i.rstrip() for i in form.cleaned_data['orthogroups'].rstrip().split('\n')]

            match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(match_groups, biodb)
            taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, match_groups, type="orthogroup")

            labels = match_groups
            tree, style = ete_motifs.multiple_profiles_heatmap(biodb, match_groups,taxon2orthogroup2count)


            big = False
            path = settings.BASE_DIR + '/assets/temp/tree.svg'
            asset_path = '/temp/tree.svg'

            tree.render(path, dpi=800, tree_style=style)

            envoi_annot = True
            envoi_annot = True
    else:  
        if display_form == "True":
            form = AnnotForm()
        else:

            from chlamdb.phylo_tree_display import ete_motifs
            from chlamdb.biosqldb import biosql_own_sql_tables


            server, db = manipulate_biosqldb.load_db(biodb)

            match_groups = target_taxons = [i for i in request.GET.getlist('g')]

            match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(match_groups, biodb)
            taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, match_groups, type="orthogroup")

            labels = match_groups
            tree, style = ete_motifs.multiple_profiles_heatmap(biodb, match_groups,taxon2orthogroup2count)


            big = False
            path = settings.BASE_DIR + '/assets/temp/tree.svg'
            asset_path = '/temp/tree.svg'

            tree.render(path, dpi=800, tree_style=style)

            envoi_annot = True

    return render(request, 'chlamdb/orthogroup_annotation.html', my_locals(locals()))

def test():
    from chlamdb.plots import plot_genomic_feature

def locus_annotation(request, display_form):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    form_class = get_LocusAnnotForm(biodb)

    if request.method == 'POST': 

        form = form_class(request.POST)
      

        #form2 = ContactForm(request.POST)
        if form.is_valid():  
            from chlamdb.biosqldb import biosql_own_sql_tables
            from chlamdb.phylo_tree_display import ete_motifs
            from ete3 import Tree
            import copy

            server, db = manipulate_biosqldb.load_db(biodb)

            match_locus = [i.rstrip() for i in form.cleaned_data['locus_list'].rstrip().split('\n')]
            circos_target = form.cleaned_data['circos_target']

            print("match_locus", match_locus)

            filter = '"'+'","'.join(match_locus)+'"'
            # inner join orthology_detail_
            # left join COG on seqfeature_id
            if db_driver == 'mysql':
                sql = 'select t1.locus_tag, t2.accession, t2.start, t2.stop, t2.gene, t2.product, t2.n_genomes, t2.orthogroup, ' \
                    ' CHAR_LENGTH(t2.translation), t4.COG_name,t6.code,t4.description, t2.taxon_id ' \
                    ' from custom_tables_locus2seqfeature_id t1 ' \
                    ' inner join orthology_detail t2 on t1.seqfeature_id=t2.seqfeature_id' \
                    ' left join COG_seqfeature_id2best_COG_hit t3 on t1.seqfeature_id=t3.seqfeature_id' \
                    ' left join COG_cog_names_2014 t4 on t3.hit_COG_id=t4.COG_id' \
                    ' left join COG_cog_id2cog_category t5 on t4.COG_id=t5.COG_id' \
                    ' left join COG_code2category t6 on t5.category_id=t6.category_id' \
                    ' where t1.locus_tag in (%s)' % (filter)
            if db_driver == 'sqlite':
                sql = 'select t1.locus_tag, t2.accession, t2.start, t2.stop, t2.gene, t2.product, t2.n_genomes, t2.orthogroup, ' \
                    ' LENGTH(t2.translation), t4.COG_name,t6.code,t4.description, t2.taxon_id ' \
                    ' from custom_tables_locus2seqfeature_id t1 ' \
                    ' inner join orthology_detail t2 on t1.seqfeature_id=t2.seqfeature_id' \
                    ' left join COG_seqfeature_id2best_COG_hit t3 on t1.seqfeature_id=t3.seqfeature_id' \
                    ' left join COG_cog_names_2014 t4 on t3.hit_COG_id=t4.COG_id' \
                    ' left join COG_cog_id2cog_category t5 on t4.COG_id=t5.COG_id' \
                    ' left join COG_code2category t6 on t5.category_id=t6.category_id' \
                    ' where t1.locus_tag in (%s)' % (filter)                
            locus_annot = [list(i) for i in server.adaptor.execute_and_fetchall(sql,)]
            locus2taxon = {}
            for row in locus_annot:
                locus2taxon[row[0]] = row[-1]

            '''
            sql = 'select seqfeature_id,t1.COG_id,function,name from COG_locus_tag2gi_hit t1 ' \
                  ' inner join COG_cog_names_2014 t2 on t1.COG_id=t2.COG_id;'


            locus2COG_data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for one_locus_annot in locus_annot:
                try:
                    one_locus_annot.append(locus2COG_data[one_locus_annot[0]][0])
                    one_locus_annot.append(locus2COG_data[one_locus_annot[0]][1])
                    one_locus_annot.append(locus2COG_data[one_locus_annot[0]][2])
                except:
                    one_locus_annot.append('-')
                    one_locus_annot.append('-')
                    one_locus_annot.append('-')

            sql = 'select locus_tag, taxon_id from orthology_detail where locus_tag in (%s)' % (biodb, filter)

            locus2taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            '''

            taxon2locus2identity_closest = ete_motifs.get_locus2taxon2identity(biodb, match_locus)

            #taxon2orthogroup2count_all = ete_motifs.get_taxon2orthogroup2count(biodb, orthogroups)
            locus2annot, \
            locus_tag2cog_catego, \
            locus_tag2cog_name, \
            locus_tag2ko, \
            pathway2category, \
            module2category, \
            ko2ko_pathways, \
            ko2ko_modules,\
            locus2interpro = get_locus_annotations(biodb, match_locus)

            labels = match_locus

            sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                       ' where t2.name="%s";' % biodb
                       
            server, db = manipulate_biosqldb.load_db(biodb)
            tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
            
            # ssubtree:
            
            acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']

            filter = '"' + '","'.join(acc_list) + '"'

            sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter

            taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            t1 = Tree(tree)
            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            
            try:
                t2 = t1.prune(taxon_list)
            except:
                pass
            t1.ladderize()


            fasta_url = '?l=' + '&l='.join(match_locus)

            tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                        labels,
                                                        taxon2locus2identity_closest,
                                                        identity_scale=True,
                                                        show_labels=False,
                                                        reference_taxon=locus2taxon,
                                                        tree=t1, 
                                                        rotate=False)


            big = False
            path = settings.BASE_DIR + '/assets/temp/tree.svg'
            asset_path = '/temp/tree.svg'
            #style2.rotation = 90
            tree2.render(path, dpi=800, tree_style=style2)


            sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                       ' where t2.name="%s";' % biodb
                       
            server, db = manipulate_biosqldb.load_db(biodb)
            
            tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
            
            acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
            filter = '"' + '","'.join(acc_list) + '"'
            
            sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
            
            taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            t1 = Tree(tree)
            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            try:
                t2 = t1.prune(taxon_list)
            except:
                pass
            t1.ladderize()

            taxon2locus2n_paralogs = ete_motifs.get_locus2taxon2n_paralogs(biodb, 
                                                                           match_locus)

            tree3, style3 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                 labels,
                                                                 taxon2locus2n_paralogs,
                                                                 reference_taxon=locus2taxon,
                                                                 tree=t1,
                                                                 identity_scale=False,
                                                                 show_labels=True,
                                                                 column_scale=True,
                                                                 as_float=False,
                                                                 rotate=False)

            show_on_circos_url = "/%s?" % circos_target +  '&l='.join(locus2taxon.keys())

            path2 = settings.BASE_DIR + '/assets/temp/tree2.svg'
            asset_path2 = '/temp/tree2.svg'
            #style3.rotation = 90
            tree3.render(path2, dpi=800, tree_style=style3)

            envoi_annot = True
            
            # plot species tree
            # retrieve species tree
            from chlamdb.phylo_tree_display import species_tree
            
            tree_complete, tree_species = species_tree.get_species_tree(biodb)
            
            sql = 'select A.locus_tag,t6.species,AVG(A.identity) from (select t1.*,t2.locus_tag ' \
                  ' from comparative_tables.identity_closest_homolog2_2019_06_PVC t1' \
                  ' inner join annotation.seqfeature_id2locus_2019_06_PVC t2 on t1.locus_1=t2.seqfeature_id ' \
                  ' where t2.locus_tag in ("%s")) A' \
                  ' inner join biosqldb.taxid2species_2019_06_PVC t4 on A.taxon_2=t4.taxon_id ' \
                  ' inner join biosqldb.species_curated_taxonomy_2019_06_PVC t6 on t4.species_id=t6.species_id ' \
                  ' group by A.locus_tag,t6.species; ' % '","'.join(match_locus)
            print(sql)
            data = server.adaptor.execute_and_fetchall(sql,)
            
            locus2taxon2identity = {}
            for row in data:
                locus_tag, species, average_identity = row
                if locus_tag not in locus2taxon2identity:
                    locus2taxon2identity[locus_tag] = {}
                locus2taxon2identity[locus_tag][species] = round(average_identity, 2)
            # group2taxon2count        
            tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                 match_locus,
                                                                 locus2taxon2identity,
                                                                 identity_scale=True,
                                                                 show_labels=True,
                                                                 reference_taxon=False,
                                                                 tree=tree_species, 
                                                                 rotate=False)


            big = False
            path = settings.BASE_DIR + '/assets/temp/tree.svg'
            asset_path = '/temp/tree.svg'
            #style2.rotation = 90
            tree2.render(path, dpi=800, tree_style=style2)
            
            tree_complete, tree_species = species_tree.get_species_tree(biodb)
            
            sql = 'select orthogroup_id, locus_tag from annotation.seqfeature_id2locus_2019_06_PVC t1 ' \
                  ' inner join orthology.seqfeature_id2orthogroup_2019_06_PVC t2 on t1.seqfeature_id=t2.seqfeature_id' \
                  ' where locus_tag in ("%s");' % '","'.join(match_locus)
            print(sql)
                  
            group_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            
            sql = 'select A.orthogroup_id,A.species, max(n_paralogs) as max_paralogs from' \
                  ' (select t2.orthogroup_id,t4.species,t1.taxon_id,count(*) as n_paralogs ' \
                  ' from annotation.seqfeature_id2locus_2019_06_PVC t1  ' \
                  ' inner join orthology.seqfeature_id2orthogroup_2019_06_PVC t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                  ' inner join biosqldb.taxid2species_2019_06_PVC t3 on t1.taxon_id=t3.taxon_id ' \
                  ' inner join biosqldb.species_curated_taxonomy_2019_06_PVC t4 on t3.species_id=t4.species_id ' \
                  ' where t2.orthogroup_id in (%s) group by t2.orthogroup_id,t4.species,t1.taxon_id) A' \
                  ' group by A.orthogroup_id,A.species;' % (','.join(group_id2locus_tag.keys()))
            print(sql)

            locus2taxon2n_paralogs = {}
            
            
            label_list = []
            for row in server.adaptor.execute_and_fetchall(sql,):
                orthogroup_id, species, max_paralog = row 
                locus_tag = "%s (%s)" % (orthogroup_id, group_id2locus_tag[str(orthogroup_id)])
                if locus_tag not in label_list:
                    label_list.append(locus_tag)
                if locus_tag not in locus2taxon2n_paralogs:
                    locus2taxon2n_paralogs[locus_tag] = {}
                locus2taxon2n_paralogs[locus_tag][species] = max_paralog
                

            tree3, style3 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                 label_list,
                                                                 locus2taxon2n_paralogs,
                                                                 reference_taxon=False,
                                                                 tree=tree_species,
                                                                 identity_scale=False,
                                                                 show_labels=True,
                                                                 column_scale=True,
                                                                 as_float=False,
                                                                 rotate=False)

            show_on_circos_url = "/%s?" % circos_target +  '&l='.join(locus2taxon.keys())

            path2 = settings.BASE_DIR + '/assets/temp/tree2.svg'
            asset_path2 = '/temp/tree2.svg'
            #style3.rotation = 90
            tree3.render(path2, dpi=800, tree_style=style3)

    else:  
        if display_form == "True":
            form = form_class()
        else:

            from chlamdb.phylo_tree_display import ete_motifs
            from chlamdb.biosqldb import biosql_own_sql_tables


            server, db = manipulate_biosqldb.load_db(biodb)

            match_groups = target_taxons = [i for i in request.GET.getlist('g')]

            match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(match_groups, biodb)
            
            taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, match_groups, type="orthogroup")

            labels = match_groups
            tree, style = ete_motifs.multiple_profiles_heatmap(biodb, 
                                                               match_groups,
                                                               taxon2orthogroup2count)


            big = False
            path = settings.BASE_DIR + '/assets/temp/tree.svg'
            asset_path = '/temp/tree.svg'

            tree.render(path, dpi=800, tree_style=style)

            envoi_annot = True

    return render(request, 'chlamdb/locus_annotation.html', my_locals(locals()))


def venn_orthogroup(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    venn_form_class = make_venn_from(db, limit=6)
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_orthogroup.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():
        return render(request, 'chlamdb/venn_orthogroup.html', my_locals(locals()))

    targets = form_venn.get_taxids()
    genomes = db.get_genomes_description()
    og_count = db.get_og_count(targets)
    fmt_data = []
    for taxon in og_count:
        ogs = og_count[taxon]
        ogs_str = ",".join(f"{to_s(format_orthogroup(og))}" for og, cnt in ogs.iteritems() if cnt > 0)
        genome = genomes.loc[int(taxon)].description
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{ogs_str}] }}")
    series = "[" + ",".join(fmt_data) + "]"

    og_list = og_count.index.tolist()
    annotations = db.get_genes_from_og(orthogroups=og_list, taxon_ids=genomes.index.tolist())
    grouped = annotations.groupby("orthogroup")
    genes = grouped["gene"].apply(list)
    products = grouped["product"].apply(list)

    orthogroup2description = []
    for og in og_list:
        forbidden = "\""
        gene_data = "-"
        if og in genes.index:
            g = genes.loc[og]
            gene_data = format_lst_to_html(g, add_count=False)
        prod_data = "-"
        if og in products.index:
            p = products.loc[og]
            prod_data = format_lst_to_html(p, add_count=False)
        og_info = "[\"" + gene_data + "\",\"" + prod_data + "\"]"
        og_item = f"h[{to_s(format_orthogroup(og))}] = {og_info};"
        orthogroup2description.append(og_item)
    orthogroup2description = "\n".join(orthogroup2description)
    envoi_venn = True
    return render(request, 'chlamdb/venn_orthogroup.html', my_locals(locals()))


def format_pfam(pfam_id):
    return f"PF{pfam_id:04d}"


def extract_pfam(request, classification="taxon_id"):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, plasmid=True, label="Pfam domains")

    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if not include_plasmids is None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if not exclude_plasmids is None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing>=sum_include_length:
        ctx = {"wrong_n_missing" : True}
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))

    pfam_include = db.get_pfam_hits(include, plasmids=include_plasmids, 
            search_on="taxid", indexing="taxid")
    if sum_exclude_length > 0:
        pfam_exclude = db.get_pfam_hits(exclude, plasmids=exclude_plasmids,
                search_on="taxid", indexing="taxid")
        pfam_exclude["sum_pos"] = pfam_exclude[cog_exclude > 0].count(axis=1)
        pfam_exclude["exclude"] = pfam_exclude.sum_pos > 0
        neg_index = pfam_exclude[pfam_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    pfam_include["sum_pos"] = pfam_include[pfam_include > 0].count(axis=1)
    pfam_include["selection"] = pfam_include.sum_pos >= len(include)-n_missing
    pos_index = pfam_include[pfam_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    pfam_defs = db.get_pfam_def(selection)

    if len(selection) == 0:
        ctx = {no_match: True}
        return render(request, 'chlamdb/extract_ko.html', my_locals(ctx))

    all_database = db.get_pfam_hits(pfam_include.index.tolist(), search_on="pfam", indexing="taxid")
    sums = all_database.sum(axis=1)
    sum_group = len(selection)

    match_groups_data = []
    for no, pfam in enumerate(selection):
        count = sums.loc[pfam]
        pfam_def = pfam_defs["def"].loc[pfam]
        data = [no+1, format_pfam(pfam), pfam_def,
                pfam_include.sum_pos.loc[pfam], sums.loc[pfam]]
        match_groups_data.append(data)

    ctx = {"envoi_extract": True,
            "sum_group": sum_group,
            "n_genomes": sum_include_length,
            "max_n": sums.max(),
            "match_groups_data": match_groups_data,
            "form": form}
    return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))


def format_ko(ko_id, as_url=False):
    base = f"K{int(ko_id):05d}"
    if not as_url:
        return base
    return f"<a href=\"/fam_ko/{base}\">{base}</a>"

def format_ko_url(ko_id):
    return format_ko(ko_id, as_url=True)


def format_ko_path(hsh_pathways, ko, as_list=False):
    pathways = hsh_pathways.get(ko, [])
    if len(pathways) == 0:
        if as_list:
            return []
        return "-"
    fmt_lst = (f"<a href=\"/KEGG_mapp_ko/map{i:05d}\">{d}</a>" for i, d in pathways)
    
    if as_list:
        return list(fmt_lst)
    return "<br>".join(fmt_lst)


def format_ko_module(module_id, module_desc=None):
    if module_desc is None:
        return f"<a href=\"/KEGG_module_map/M{module_id:05d}\">M{module_id:05d}</a>"
    else:
        return f"<a href=\"/KEGG_module_map/M{module_id:05d}\">{module_desc}</a>"


def format_ko_modules(hsh_modules, ko):
    modules = hsh_modules.get(ko, [])
    if len(modules) == 0:
        return "-"
    return "<br>".join([format_ko_module(i, d) for i, d in modules])


def extract_ko(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, plasmid=True, label="Kegg Orthologs")

    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if not include_plasmids is None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if not exclude_plasmids is None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing>=sum_include_length:
        hsh_var = {"wrong_n_missing": True, "form": form}
        return render(request, 'chlamdb/extract_ko.html', my_locals(hsh_var))

    mat_include = db.get_ko_hits(include, plasmids=include_plasmids)
    if len(exclude) > 0:
        mat_exclude = db.get_ko_hits(exclude, plasmids=exclude_plasmids)
        mat_exclude["sum_pos"] = mat_exclude[mat_exclude > 0].count(axis=1)
        mat_exclude["exclude"] = mat_exclude.sum_pos > 0
        neg_index = mat_exclude[mat_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    mat_include["sum_pos"] = mat_include[mat_include > 0].count(axis=1)
    mat_include["selection"] = mat_include.sum_pos >= len(include)-n_missing
    pos_index = mat_include[mat_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        hsh_var = {"no_match": True, "form": form}
        return render(request, 'chlamdb/extract_ko.html', my_locals(hsh_var))

    all_database = db.get_ko_hits(selection, search_on="ko", indexing="taxid")
    ko_total_count = all_database.sum(axis=1)
    ko_desc = db.get_ko_desc(selection)
    ko_mod = db.get_ko_modules(selection)
    ko_path = db.get_ko_pathways(selection)
    match_groups_data = []
    for ko in selection:
        kof = format_ko(ko)
        kod = ko_desc.get(ko, "-")
        kop = format_ko_path(ko_path, ko)
        kom = format_ko_modules(ko_mod, ko)
        kot = ko_total_count.loc[ko]
        data = [kof, kod, kop, kom, mat_include.sum_pos.loc[ko], kot]
        match_groups_data.append(data)

    # should contain the loci of the KO in the reference genome
    max_n = ko_total_count.max()
    locus_list = [] # [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
    target_circos_taxons = include + exclude

    # url to get the barchart of selected KO
    taxons_in_url = "?i="+("&i=").join(map(str, include)) + '&m=%s' % str(n_missing)
    taxon_out_url = "&o="+("&o=").join(map(str, exclude))
    envoi_extract = True
    mm = 'module'
    pp = 'pathway'

    locus2annot = {}
    locus_tag2cog_catego = {}
    locus_tag2cog_name = {}
    locus_tag2ko = {}
    pathway2category = {}
    module2category = {}
    ko2ko_pathways = {}
    ko2ko_modules = {}
    locus2interpro = {}
    return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))


def extract_EC(request):
    biodb = settings.BIODB

    '''

    :param request:
    :param biodb:
    :param classification: either taxon_id or accession (merging plasmids or not)
    :return:
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    extract_form_class = make_extract_form(biodb, label="EC numbers")

    if request.method == 'POST': 

        form = extract_form_class(request.POST)

        if form.is_valid():  
            from chlamdb.biosqldb import biosql_own_sql_tables

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
                                                                              freq_missing,
                                                                              cache=cache)

                match_groups = mat.index.tolist()

                if len(match_groups) == 0:
                    no_match = True
                else:

                    # get count in subgroup
                    ec2count = dict((mat > 0).sum(axis=1))
                    # get count in complete database
                    ec2count_all = dict((mat_all > 0).sum(axis=1))

                    max_n = max(ec2count_all.values())

                    # GET max frequency for template
                    sum_group = len(match_groups)


            sql = 'select ec, value, pathway_name, pathway_category, description from ' \
            ' (select enzyme_id, ec,value from enzyme_enzymes as t1 inner join enzyme_enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id ' \
            ' where line="description") A left join enzyme_kegg2ec as B on A.enzyme_id=B.ec_id ' \
            ' left join enzyme_kegg_pathway on B.pathway_id=kegg_pathway.pathway_id;'
            
            ec2description_raw = server.adaptor.execute_and_fetchall(sql,)

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

            locus_list_sql = 'select locus_tag from (select taxon_id,locus_tag,ec_id from enzyme_locus2ec as t1  ' \
                             ' inner join bioentry as t2 on t1.accession=t2.accession ' \
                             ' where biodatabase_id=%s) A inner join enzyme_enzymes as B on A.ec_id=B.enzyme_id' \
                             ' where A.taxon_id=%s and B.ec in (%s);' % (biodb_id,
                                                                         reference_taxon,
                                                                         EC_list)

            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]

            circos_url = '?ref=%s&' % reference_taxon
            circos_url+= "t="+('&t=').join((include + exclude)) + '&h=' + ('&h=').join(locus_list)

            target_circos_taxons = include + exclude


            # get phylogenetic profile of match if not too big
            if len(match_groups) < 50:
                from chlamdb.phylo_tree_display import ete_motifs
                sql = 'select distinct ec,orthogroup from enzyme_locus2ec as t1 ' \
                      ' inner join enzyme_enzymes as t2 on t1.ec_id=t2.enzyme_id where ec in (%s);' % ('"' + '","'.join(match_groups) + '"')
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
                tree2, style2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                                     labels,
                                                                     taxon2orthogroup2count,
                                                                     taxon2enzyme2count,
                                                                     ec2orthogroups)



                if len(labels) > 40:
                    big = True
                    path = settings.BASE_DIR + '/assets/temp/profil_tree.png'
                    asset_path = '/temp/profil_tree.png'
                    tree2.render(path, dpi=1200)



                else:
                    big = False

                    path2 = settings.BASE_DIR + '/assets/temp/profil_tree.svg'
                    asset_path = '/temp/profil_tree.svg'

                    tree2.render(path2, dpi=800, tree_style=style2)





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






    else:  
        form = extract_form_class()

    return render(request, 'chlamdb/extract_EC.html', my_locals(locals()))


def venn_pfam(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    venn_form_class = make_venn_from(db, limit=6)
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_Pfam.html', my_locals({"form_venn": form_venn}))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():  
        # add error message
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_Pfam.html', my_locals(locals()))

    targets = form_venn.get_taxids()
    pfam_hits = db.get_pfam_hits(targets, search_on="taxid", indexing="taxid")
    data = db.get_pfam_def(pfam_hits.index.tolist())
    genomes_desc = db.get_genomes_description().description.to_dict()

    series_tab = []
    for target in targets:
        pfams = pfam_hits[target]
        non_zero = pfams[pfams > 0]
        str_fmt = ",".join(f"\"{format_pfam(pfam)}\"" for pfam, _ in non_zero.iteritems())
        series_tab.append(f"{{name: \"{genomes_desc[target]}\", data: [{str_fmt}]}}")
    series = "[" + ",".join(series_tab) + "]"

    descriptions = []
    for pfam, pfam_info in data.iterrows():
        pfam_def = pfam_info["def"]
        descriptions.append(f"h[\"{format_pfam(pfam)}\"] = \"{pfam_def}\"")

    ctx = {"envoi_venn": True,
            "series": series,
            "pfam2description": ";".join(descriptions)}
    return render(request, 'chlamdb/venn_Pfam.html', my_locals(ctx))


def venn_EC(request):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    venn_form_class = make_venn_from(biodb, limit=6)
    if request.method == 'POST': 

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
                sql ='select id from comparative_tables_EC where `%s` > 0' % (target)
                ec_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_ec_list += ec_list
                data = '"' + '","'.join(ec_list) + '"'
                series+=template_serie % (taxon_id2genome[target], ec_list) + ','
            series = series[0:-1] + ']'


            ec2description = ''
            sql = 'select ec, value, pathway_name, pathway_category, description from ' \
                  ' (select enzyme_id, ec,value from enzyme_enzymes as t1 inner join enzyme_enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id ' \
                  ' where line="description") A left join enzyme_kegg2ec as B on A.enzyme_id=B.ec_id ' \
                  ' left join enzyme_kegg_pathway on B.pathway_id=kegg_pathway.pathway_id;'
                  
            ec2description_raw = server.adaptor.execute_and_fetchall(sql,)
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


            envoi_venn = True
    else:  
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_EC.html', my_locals(locals()))



def extract_interpro(request, classification="taxon_id"):
    biodb = settings.BIODB
    '''

    :param request:
    :param biodb:
    :param classification: either taxon_id or accession (merging plasmids or not)
    :return:
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    extract_form_class = make_extract_form(biodb, label="Interpro entries", plasmid=True)

    if request.method == 'POST': 

        form = extract_form_class(request.POST)

        #form2 = ContactForm(request.POST)
        if form.is_valid():  
            print("valid form!")
            from chlamdb.biosqldb import biosql_own_sql_tables

            include = form.cleaned_data['orthologs_in']
            exclude = form.cleaned_data['no_orthologs_in']
            n_missing = form.cleaned_data['frequency']
            reference_taxon = form.cleaned_data['reference']
            if reference_taxon == "None":
                reference_taxon = include[0]

            freq_missing = (len(include)-float(n_missing))/len(include)

            try:
                accessions = request.POST['checkbox_accessions']
                accessions = True
                fasta_url='?a=T'
            except:
                print("accessions = False")
                accessions = False
                fasta_url='?a=F'
                accession2taxon = manipulate_biosqldb.accession2taxon_id(server, biodb)
                
                include = [str(accession2taxon[i]) for i in include]
                exclude = [str(accession2taxon[i]) for i in exclude]
                
                reference_taxon = accession2taxon[reference_taxon]

            
            print("running task!")
            task = extract_interpro_task.delay(biodb, 
                                                include,
                                                exclude,
                                                freq_missing,
                                                reference_taxon,
                                                accessions,
                                                n_missing)
            task_id = task.id

            print("task ok!")

            return HttpResponse(json.dumps({'task_id': task.id}), content_type='application/json')
        else:
            return HttpResponse(json.dumps({'task_id': None}), content_type='application/json')

    else:  
        form = extract_form_class()

    return render(request, 'chlamdb/extract_interpro.html', my_locals(locals()))


def venn_interpro(request):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    venn_form_class = make_venn_from(biodb, limit=6)
    if request.method == 'POST': 

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
                sql ='select id from comparative_tables_interpro where `%s` > 0' % (target)
                cogs = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                all_pfam_list += cogs
                data = '"' + '","'.join(cogs) + '"'
                series+=template_serie % (taxon_id2genome[target], cogs) + ','
            series = series[0:-1] + ']'


            interpro2description = ''
            sql = 'select interpro_accession,interpro_description, count(*) from interpro group by interpro_accession,interpro_description;'
            data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            for i in data:
                if i in all_pfam_list:
                    interpro2description+='h["%s"] = "%s </td><td>%s";' % (i, data[i][0], data[i][1])
                else:
                    print ('pas ok')

            envoi_venn = True
    else:  
        form_venn = venn_form_class()
    return render(request, 'chlamdb/venn_interpro.html', my_locals(locals()))


def extract_cog(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, plasmid=True, label="COG")

    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if not include_plasmids is None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if not exclude_plasmids is None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing>=sum_include_length:
        wrong_n_missing = True
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    cog_include = db.get_cog_hits(include, plasmids=include_plasmids, 
            search_on="taxid", indexing="taxid")
    if sum_exclude_length > 0:
        cog_exclude = db.get_cog_hits(exclude, plasmids=exclude_plasmids,
                search_on="taxid", indexing="taxid")
        cog_exclude["sum_pos"] = cog_exclude[cog_exclude > 0].count(axis=1)
        cog_exclude["exclude"] = cog_exclude.sum_pos > 0
        neg_index = cog_exclude[cog_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    cog_include["sum_pos"] = cog_include[cog_include > 0].count(axis=1)
    cog_include["selection"] = cog_include.sum_pos >= len(include)-n_missing
    pos_index = cog_include[cog_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        no_match = True
        return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))

    all_database = db.get_cog_hits(cog_include.index.tolist(), search_on="cog", indexing="taxid")
    sums = all_database.sum(axis=1)

    cat_count = {}
    cogs_summaries = db.get_cog_summaries(sums.index.tolist())
    cogs_funct = db.get_cog_code_description()
    cog_data = []
    for cog_id in selection:
        count = sums.loc[cog_id]

        # some cogs do not have a description, skip those
        if cog_id not in cogs_summaries:
            continue

        data = [format_cog(cog_id)]
        func_acc = []
        for func, func_descr, cog_descr in cogs_summaries[cog_id]:
            func_acc.append((func, func_descr))
            inc, not_incl = cat_count.get(func, (0, 0))
            cat_count[func] = (inc+cog_include.sum_pos.loc[cog_id], not_incl)
        funcs = "<br>".join(f"{func} ({func_desc})" for func, func_desc in func_acc)
        data = (format_cog(cog_id), funcs, cog_descr,
                cog_include.sum_pos.loc[cog_id], str(count))
        cog_data.append(data)

    # get the categories for all cogs
    for cog_id, details_lst in cogs_summaries.items():
        for func, func_descr, cog_descr in details_lst:
            inc, not_incl = cat_count.get(func, (0, 0))
            cat_count[func] = (inc, not_incl+sums.loc[cog_id])

    max_n = sums.max(axis=0)
    sum_group = len(selection)

    # Code to generate the barchart diagrams
    cat_map_str = ",".join([f"\"{func}\": \"{descr}\"" for func, descr in cogs_funct.items()])
    category_map = f"var category_description = {{ {cat_map_str} }};"
    ttl_sel = sum([c1 for func, (c1, c2) in cat_count.items()])
    ttl_all = sum([c2 for func, (c1, c2) in cat_count.items()])

    cat_count_comp = ",".join([f"\"{func}\" : [\"{c2}\", \"{c1}\"]" for func, (c1, c2) in cat_count.items()])
    category_count_complete = f"var category_count_complete = {{ {cat_count_comp} }};"

    serie_selection_val = [str(round(float(c1)/ttl_sel, 2)) for func, (c1, c2) in cat_count.items()]
    serie_all_val = [str(round(float(c2)/ttl_all, 2)) for func, (c1, c2) in cat_count.items()]
    serie_selection = f"{{ labels: \"selection\", values: [{','.join(serie_selection_val)}] }}"
    serie_all = f"{{ labels: \"complete genomes\", values: [{','.join(serie_all_val)}] }}"
    series_str = ",".join([serie_all, serie_selection])
    series = f"[ {series_str} ]"
    labels_str = ",".join([f"\"{funct}\"" for funct in cat_count.keys()])
    labels = f"[ {labels_str} ]"

    locus2annot = {}
    locus_tag2cog_catego = {}
    locus_tag2cog_name = {}
    locus_tag2ko = {}
    pathway2category = {}
    module2category = {}
    ko2ko_pathways = {}
    ko2ko_modules = {}
    locus2interpro = {} #  get_locus_annotations(biodb, locus_list)
    envoi_extract = True
    return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))


def venn_ko(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    venn_form_class = make_venn_from(db, limit=6)
    display_form = True
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():  
        # add error message
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    taxids    = form_venn.get_taxids()
    genomes   = db.get_genomes_description().description.to_dict()
    ko_counts = db.get_ko_count(taxids)

    fmt_data = []
    ko_list = ko_counts.index.get_level_values("KO").unique().to_list()
    for taxid in taxids:
        kos = ko_counts.loc[taxid].index.values
        kos_str = ",".join(f"{to_s(format_ko(ko))}" for ko in kos)
        genome = genomes[taxid]
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{kos_str}] }}")
    series = "[" + ",".join(fmt_data) + "]"

    ko_descriptions = db.get_ko_desc(ko_list)
    ko2description = []
    for ko, ko_desc in ko_descriptions.items():
        forbidden = "\""
        ko_item = f"h[{to_s(format_ko(ko))}] = [{forbidden}{ko_desc}{forbidden}];"
        ko2description.append(ko_item)
    ko2description = "\n".join(ko2description)
    envoi_venn = True
    return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))


def format_cog(cog_id, as_url=False):
    base = f"COG{int(cog_id):04d}"
    if as_url==False:
        return base
    return f"<a href=\"/fam_cog/{base}\">{base}</a>"

def format_cog_url(cog_id):
    return format_cog(cog_id, as_url=True)

def venn_cog(request, sep_plasmids=False):
    """
    Will need to modify the signature of the method to remove the sep_plasmid 
    parameter as it is not taken into account. Or put back the differentiate
    plasmid parameter in the web page.
    """

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    display_form = True
    venn_form_class = make_venn_from(db, limit=6, label="COG")
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():
        # TODO: add error message
        return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))

    targets = form_venn.get_taxids()
    cog_hits = db.get_cog_hits(targets, indexing="taxid", search_on="taxid")
    data = db.get_cog_summaries(cog_hits.index.tolist(), only_cog_desc=True, as_df=True)
    genome_desc = db.get_genomes_description().description.to_dict()

    # necessary as some COG do not have a description
    # --> filter them out
    cog_hits = cog_hits.reindex(data.index)

    series_tab = []
    for target in targets:
        cogs = cog_hits[target]
        non_zero_cogs = cogs[cogs > 0]
        str_fmt = ",".join(f"\"{format_cog(cog)}\"" for cog, count in non_zero_cogs.iteritems())
        series_tab.append( f"{{name: \"{genome_desc[target]}\", data: [{str_fmt}]}}" )
    series = "[" + ",".join(series_tab) + "]"

    cog2description_l = []
    cog_codes = db.get_cog_code_description()
    for cog, data in data.iterrows():
        name = data.description
        func = data.function
        functions = ",".join(f"\"{abbr}\"" for abbr in func)
        cog2description_l.append(f"h[\"{format_cog(cog)}\"] = [[{functions}], \"{name}\"]")

    cog_func_dict = (f"\"{func}\": \"{descr}\"" for func, descr in cog_codes.items())
    cog_func_dict = "{"+",".join(cog_func_dict)+"}"
    cog2description = ";".join(cog2description_l)
    envoi_venn = True
    return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))


def pmid_associations(request, bioentry_id, pmid, data_type):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if data_type == "PaperBlast":
        sql = 'select t1.description from biosqldb.bioentry t1 inner join biosqldb.biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id where t1.bioentry_id=%s and t2.name="%s"' % (bioentry_id, biodb)
        genome = server.adaptor.execute_and_fetchall(sql,)[0][0]

        sql = 'select title,journal,year from string.pmid2data_paperblast where pmid=%s;' % (pmid)
        paper_data = server.adaptor.execute_and_fetchall(sql,)[0]

        sql = 'select distinct locus_tag from string.seqfeature_id2paperblast t1 ' \
            ' inner join string.paperblast2pmid t2 on t1.paperblast_id=t2.paperblast_id ' \
            ' inner join annotation.seqfeature_id2locus_%s t3 on t1.seqfeature_id=t3.seqfeature_id ' \
            ' where PMID=%s and bioentry_id=%s;' % (biodb, pmid, bioentry_id)

    if data_type == "STRING":
        sql = 'select t1.description from biosqldb.bioentry t1 ' \
              ' inner join biosqldb.biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id where t1.bioentry_id=%s and t2.name="%s"' % (bioentry_id, biodb)
        genome = server.adaptor.execute_and_fetchall(sql,)[0][0]

        sql = 'select title,journal,year from string.pmid2data_stringdb where pmid=%s;' % (pmid)
        paper_data = server.adaptor.execute_and_fetchall(sql,)[0]

        sql = 'select distinct locus_tag from string.seqfeature_id2string_protein_mapping t1 ' \
            ' inner join string.string_protein2pmid t2 on t1.string_protein_id=t2.string_protein_id ' \
            ' inner join annotation.seqfeature_id2locus_%s t3 on t1.seqfeature_id=t3.seqfeature_id ' \
            ' where PMID=%s and bioentry_id=%s;' % (biodb, pmid, bioentry_id)

    locus_list =[i[0] for i in server.adaptor.execute_and_fetchall(sql,)] 

    locus2annot, \
    locus_tag2cog_catego, \
    locus_tag2cog_name, \
    locus_tag2ko, \
    pathway2category, \
    module2category, \
    ko2ko_pathways, \
    ko2ko_modules, \
    locus2interpro = get_locus_annotations(biodb, locus_list)


    return render(request, 'chlamdb/pmid_associations.html', my_locals(locals()))


def pmid_associations_orthogroups(request, pmid, data_type):
    from chlamdb.biosqldb import biosql_own_sql_tables

    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if data_type == "PaperBlast":
        sql = 'select t1.description from biosqldb.bioentry t1 inner join biosqldb.biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id where t1.bioentry_id=%s and t2.name="%s"' % (bioentry_id, biodb)
        genome = server.adaptor.execute_and_fetchall(sql,)[0][0]

        sql = 'select title,journal,year from string.pmid2data_paperblast where pmid=%s;' % (pmid)
        paper_data = server.adaptor.execute_and_fetchall(sql,)[0]

        sql = 'select distinct locus_tag from string.seqfeature_id2paperblast t1 ' \
            ' inner join string.paperblast2pmid t2 on t1.paperblast_id=t2.paperblast_id ' \
            ' inner join annotation.seqfeature_id2locus_%s t3 on t1.seqfeature_id=t3.seqfeature_id ' \
            ' where PMID=%s and bioentry_id=%s;' % (biodb, pmid, bioentry_id)

    if data_type == "STRING":
        sql = 'select t1.pmid,orthogroup_name from string.pmid2data_stringdb t1 ' \
                ' inner join string.string_protein2pmid t2 on t1.pmid=t2.pmid ' \
                ' inner join string.seqfeature_id2string_protein_mapping t3 on t2.string_protein_id=t3.string_protein_id ' \
                ' inner join orthology.seqfeature_id2orthogroup_%s t4 on t3.seqfeature_id=t4.seqfeature_id ' \
                ' inner join orthology.orthogroup_%s t5 on t4.orthogroup_id=t5.orthogroup_id' \
                ' where t1.pmid in (%s) group by t1.pmid,orthogroup_name' % (biodb, biodb, pmid)

        orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql2 = 'select title,journal,year from string.pmid2data_stringdb where pmid=%s;' % (pmid)
        paper_data = server.adaptor.execute_and_fetchall(sql2,)[0]

    orthogroup_list =[i[1] for i in server.adaptor.execute_and_fetchall(sql,)] 

    match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(orthogroup_list,
                                                                                                  biodb,
                                                                                                  taxon_filter=False,
                                                                                                  accessions=False)


    return render(request, 'chlamdb/pmid_associations_orthogroups.html', my_locals(locals()))




def pmid(request, seqfeature_id):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = f'select distinct accession,db_name,organism,description,query_cov,hit_cov,identity,evalue,score,t4.* from string.seqfeature_id2paperblast t1' \
            f' inner join string.paperblast_entry t2 on t1.paperblast_id=t2.id inner join string.paperblast2pmid t3 on t1.paperblast_id=t3.paperblast_id' \
            f' inner join string.pmid2data_paperblast t4 on t3.pmid=t4.pmid where t1.seqfeature_id={seqfeature_id};'
    
    print(sql1)

    pubmed_count = 0
    try:
        paperblast_data = server.adaptor.execute_and_fetchall(sql1,)
    except:
        paperblast_data = False

    sql1 = f'select bioentry_id,locus_tag from annotation_seqfeature_id2locus where seqfeature_id={seqfeature_id}' 

    data = server.adaptor.execute_and_fetchall(sql1,)[0]
    bioentry_id = data[0]
    locus_tag = data[1]

    if paperblast_data:
        # retrieve bioentry

        pmid_list =[str(i[9]) for i in paperblast_data] 
        pmid_filter = ','.join(pmid_list)

        sql_2 = 'select pmid,count(*) from (select distinct pmid,t3.seqfeature_id from string_seqfeature_id2paperblast t1 ' \
                ' inner join string_paperblast2pmid t2 on t1.paperblast_id=t2.paperblast_id ' \
                ' inner join annotation_seqfeature_id2locus_%s t3 on t1.seqfeature_id=t3.seqfeature_id ' \
                ' where bioentry_id=%s) A where A.pmid in (%s) group by pmid;'  % (bioentry_id, 
                                                                                   pmid_filter)
                
        pmid2n_associated_proteins = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_2,))
        pmid2n_associated_proteins = {int(k):int(v) for k,v in pmid2n_associated_proteins.items()}

    sql2 = f'select distinct accession,organism,description,query_cov,hit_cov,identity,evalue,score,t4.* ' \
           f' from string_seqfeature_id2string_protein_mapping t1 ' \
           f' inner join string_string_protein_entry t2 on t1.string_protein_id=t2.id ' \
           f' inner join string_string_protein2pmid t3 on t1.string_protein_id=t3.string_protein_id ' \
           f' inner join string_pmid2data_stringdb t4 on t3.pmid=t4.pmid where t1.seqfeature_id={seqfeature_id};'

    print(sql2)    
    try:
        string_data = server.adaptor.execute_and_fetchall(sql2,)
    except:
        string_data = False


    if string_data:
        
        pmid_list =[str(i[8]) for i in string_data] 
        pmid_filter = ','.join(pmid_list)

        sql_2 = 'select pmid,count(*) from (select distinct pmid,t3.seqfeature_id from string.seqfeature_id2string_protein_mapping t1 ' \
                ' inner join string.string_protein2pmid t2 on t1.string_protein_id=t2.string_protein_id ' \
                ' inner join annotation.seqfeature_id2locus_%s t3 on t1.seqfeature_id=t3.seqfeature_id ' \
                ' where bioentry_id=%s) A where A.pmid in (%s) group by pmid;'  % (biodb, bioentry_id, pmid_filter)
        
        pmid2n_associated_proteins_string = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_2,))
        pmid2n_associated_proteins_string = {int(k):int(v) for k,v in pmid2n_associated_proteins_string.items()}
        print(pmid2n_associated_proteins_string)
    return render(request, 'chlamdb/pmid.html', my_locals(locals()))



def extract_region(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100*x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x/pow(10,6), 2))

    filenames_tax_id= db.get_filenames_to_taxon_id()
    filenames_tax_id_db= pd.DataFrame.from_dict(list(filenames_tax_id.items()))
    filenames_tax_id_db.columns = ['filename','taxon_id']
    filenames_tax_id_db.index= list(filenames_tax_id_db['taxon_id'])
    filenames_list= list(filenames_tax_id_db["filename"])
    path_pre="temp/"
    path_suf_faa=".faa"
    path_faa=[path_pre + filename + path_suf_faa for filename in filenames_list]

    path_suf_fna=".fna"
    path_fna=[path_pre + filename + path_suf_fna for filename in filenames_list]

    path_suf_ffn=".ffn"
    path_ffn=[path_pre + filename + path_suf_ffn for filename in filenames_list]

    path_suf_gbk=".gbk"
    path_gbk=[path_pre + filename + path_suf_gbk for filename in filenames_list]

    filenames_tax_id_db['path_to_faa']= path_faa
    filenames_tax_id_db['path_to_fna']= path_fna
    filenames_tax_id_db['path_to_ffn']= path_ffn
    filenames_tax_id_db['path_to_gbk']= path_gbk
    print("path_gbk",path_gbk)
    filenames_tax_id_db = filenames_tax_id_db[["path_to_faa", "path_to_fna", "path_to_ffn", "path_to_gbk" ]]
    genomes_data=genomes_data.join(filenames_tax_id_db, on= "taxon_id")
    print("genomes_data", genomes_data)


    data_table_header = [ "Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding", "N plasmid contigs", "faa seq", "fna seq", "ffn seq", "gbk file"]
    data_table = genomes_data[[ "id", "description", "gc", "n_prot", "n_contigs", "length", "coding_density", "has_plasmid", "path_to_faa", "path_to_fna", "path_to_ffn", "path_to_gbk" ]].values.tolist()

    ext_list= [".faa", ".fna", ".ffn"]
    for i in list(filenames_tax_id.keys()): 
        fasta_src_faa = settings.FOLDER_PATH + "blast_DB/faa/faa_SEQ/" + i + ".faa"
        fasta_dst_faa = settings.BASE_DIR + "/assets/temp/" + i + ".faa"
      
        fasta_src_fna = settings.FOLDER_PATH + "blast_DB/fna/fna_SEQ/" + i + ".fna"
        fasta_dst_fna = settings.BASE_DIR + "/assets/temp/" + i + ".fna"

        fasta_src_ffn = settings.FOLDER_PATH + "blast_DB/ffn/ffn_seq/" + i + ".ffn"
        fasta_dst_ffn = settings.BASE_DIR + "/assets/temp/" + i + ".ffn"

        src_gbk = settings.FOLDER_PATH + "data/prokka_output_filtered/"  + i + ".gbk"
        dst_gbk = settings.BASE_DIR + "/assets/temp/" + i + ".gbk"
        try:
            os.symlink(fasta_src_faa, fasta_dst_faa)
            os.symlink(fasta_src_fna, fasta_dst_fna)
            os.symlink(fasta_src_ffn, fasta_dst_ffn)
            os.symlink(src_gbk, dst_gbk)
            break
        except FileExistsError:
        
                os.remove(fasta_dst_faa)
                os.remove(fasta_dst_fna)
                os.remove(fasta_dst_ffn)
                os.remove(dst_gbk)
                os.symlink(fasta_src_faa, fasta_dst_faa)
                os.symlink(fasta_src_fna, fasta_dst_fna)
                os.symlink(fasta_src_ffn, fasta_dst_ffn)
                os.symlink(src_gbk, dst_gbk)

       
    return render(request, 'chlamdb/extract_region.html', my_locals(locals()))

def extract_contigs(request, genome):
    
    taxid = int(genome)
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    foo = db.get_proteins_info([taxid], search_on="taxid", as_df=True)
    seqids = foo.index.tolist()

    ogs = db.get_og_count(seqids, search_on="seqid")
    taxid_seqid = db.get_taxid_from_seqid(seqids)
    taxid_seqid_db = pd.DataFrame.from_dict(list(taxid_seqid.items()))
    taxid_seqid_db.columns = ['seqid','taxon_id']

    loc = db.get_gene_loc(seqids, as_hash=True) 
    loc_db= pd.DataFrame.from_dict(list(loc.items()))
    loc_db.columns = ['seqid','all']
    loc_info= loc_db["all"]
    loc_info.columns = ['all']
    loc_info = pd.DataFrame(loc_info.to_list(), columns=['strand','start', 'stop'])
    loc_info.strand = [format(num).rstrip('0').rstrip('.')for num in loc_info.strand]
    loc_info.start = [format(num).rstrip('0').rstrip('.')for num in loc_info.start] 
    loc_info.stop = [format(num).rstrip('0').rstrip('.')for num in loc_info.stop] 

    contigs=db.get_contigs_to_seqid(taxid)
    bar = foo.join(ogs).join(taxid_seqid_db).join(loc_info).join(contigs)

    data_table_header = [ "Name", "Gene",   "Product",  "Locus_tag", "Orthogroup", "Contig" ,"Strand", "Start", "Stop", ]
    data_table = bar[["description", "gene",  "product", "locus_tag", "orthogroup", "contig" ,"strand", "start", "stop" ]].values.tolist()
    
    return render(request, 'chlamdb/extract_contigs.html', my_locals(locals()))

def format_lst(lst):
    hsh_values = {}
    for item in lst:
        val = hsh_values.get(item, 0)
        hsh_values[item] = val+1
    return hsh_values


def tab_homologs(db, infos, hsh_organism, ref_seqid=None, og=None):
    n_genomes = len(set(hsh_organism.values()))
    if n_genomes==1:
        n_genomes = "1 genome"
    else:
        n_genomes = f"{n_genomes} genomes"

    headers = ["", "Locus tag", "Source", "Gene", "Product"]
    identities = None
    if ref_seqid!=None:
        identities = db.get_og_identity(og, ref_seqid)
        headers.insert(2 , "Identity")

    homologues = []
    index = 0
    for seqid, data in infos.iterrows():
        organism = hsh_organism[seqid]
        locus_fmt = format_locus(data.locus_tag, to_url=True)
        entry = [index+1, locus_fmt, organism, format_gene(data.gene), data["product"]]
        if ref_seqid!=None:
            if seqid==ref_seqid:
                ident = "-"
            else:
                ident = round(identities.loc[seqid].identity, 1)
            entry.insert(2, ident)

        homologues.append(entry)
        index += 1

    return {"orthogroup": orthogroup,
            "n_genomes": n_genomes,
            "headers": headers,
            "homologues": homologues }


def tab_lengths(n_homologues, annotations):
    import plotly.figure_factory as ff

    length_distrib = n_homologues > 1
    if not length_distrib:
        return {"length_distrib": False}

    lengths = annotations["length"]
    max_protein_length = lengths.max()
    std_protein_length = f"{lengths.std():.1f}"
    min_protein_length = lengths.min()
    mean_protein_length = f"{lengths.mean():.1f}"
    median_protein_length = f"{lengths.median():.1f}"
    if len(lengths.unique()) > 1:
        fig1 = ff.create_distplot([lengths.tolist()], ["Sequence length"], bin_size=20)
        fig1.update_xaxes(range=[0, max_protein_length])
        fig1.layout.margin.update({"l": 80, "r": 20, "b": 40, "t": 20, "pad": 10, })
        html_plot_prot_length = manipulate_biosqldb.make_div(fig1, div_id="distplot")
    else:
        return { "length_distrib": True, "single_length": True, "prot_length": lengths.iloc[0]}

    return {"length_distrib": True,
            "max_protein_length": max_protein_length,
            "std_protein_length": std_protein_length,
            "min_protein_length": min_protein_length,
            "mean_protein_length": mean_protein_length,
            "median_protein_length": median_protein_length,
            "html_plot_prot_length": html_plot_prot_length }


class SimpleTextColumn(Column):
    def __init__(self, header=None):
        super().__init__(header)

    def get_face(self, index):
        return TextFace(index, fsize=7)


def tab_og_phylogeny(db, og_id, annot):
    og_phylogeny = db.get_og_phylogeny(og_id)

    tree = Tree(og_phylogeny)

    locuses = [branch.name for branch in tree.iter_leaves()]
    locus_to_genome = db.get_locus_to_genomes(locuses)
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    e_tree.add_column(SimpleTextColumn("Locus tag"))
    e_tree.rename_leaves(locus_to_genome, leaf_name_type=str)

    asset_path = f"/temp/og_phylogeny{og_id}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    e_tree.render(path, dpi=1200)
    return {"og_phylogeny": asset_path}


def tab_og_conservation_tree(db, group, compare_to=None):
    from metagenlab_libs.ete_phylo import EteTree, SimpleColorColumn, ModuleCompletenessColumn
    from ete3 import Tree
    
    ref_phylogeny = db.get_reference_phylogeny()
    leaf_to_name = db.get_genomes_description().description.to_dict()
    
    # Note: as the orthogroup may either be in a plasmid or in the chromosome
    # of the bacteria, we need to index by taxon to group them (index on taxon)
    count = db.get_og_count([group], search_on="orthogroup")

    tree = Tree(ref_phylogeny)
    R = tree.get_midpoint_outgroup()
    if not R is None:
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    e_tree.add_column(SimpleColorColumn.fromSeries(count.loc[group], header=format_orthogroup(group)))
    if not compare_to is None:
        identity_matrix = db.get_og_identity(group, compare_to)
        seqids = identity_matrix.index.tolist()

        # to get the taxid of the reference seqid, so as to exclude it from 
        # the phylogenetic tree
        seqids.append(compare_to)
        seqid_to_taxon  = db.get_taxid_from_seqid(seqids)
        identity_matrix["taxid"] = identity_matrix.index.map(seqid_to_taxon)
        max_identity = identity_matrix.groupby("taxid").max()
        col = LocusHeatmapColumn.fromSeries(max_identity["identity"],
                header="Identity", cls=LocusHeatmapColumn)
        col.ref_taxon = seqid_to_taxon[compare_to]
        e_tree.add_column(col)

    e_tree.rename_leaves(leaf_to_name)

    dpi = 1200
    asset_path = f"/temp/og_conservation{group}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    e_tree.render(path, dpi=dpi)
    return {"asset_path": asset_path}


def og_tab_get_kegg_annot(db, seqids):
    ko_hits = db.get_ko_hits(seqids, search_on="seqid", indexing="seqid")
    if ko_hits.empty:
        return {}

    n_occurences = ko_hits["ko"].value_counts()
    ko_ids = n_occurences.index.tolist()
    ko_descr = db.get_ko_desc(ko_ids)
    ko_pathways = db.get_ko_pathways(ko_ids)
    ko_modules = db.get_ko_modules(ko_ids)
        
    ko_entries = []
    for ko_id, count in n_occurences.iteritems():
        entry = [format_ko(ko_id, as_url=True), count]
        descr = ko_descr.get(ko_id, "-")
        entry.append(descr)
        entry.append(format_ko_path(ko_pathways, ko_id))
        entry.append(format_ko_modules(ko_modules, ko_id))
        ko_entries.append(entry)

    ko_header = ["KO", "Occurences", "Description", "Modules", "Pathways"]
    return {
        "ko_entries": ko_entries,
        "ko_header": ko_header
    }


def og_tab_get_cog_annot(db, seqids):
    cog_hits = db.get_cog_hits(seqids, indexing="seqid", search_on="seqid")

    if cog_hits.empty:
        return {}

    n_entries = cog_hits["cog"].value_counts()
    cog_summ = db.get_cog_summaries(n_entries.index.tolist())
    cog_entries = []
    for cog_id, count in n_entries.iteritems():
        entry = [format_cog(cog_id, as_url=True), count]
        funcs = []
        func_descrs = []
        cog_descrs = []
        for func, func_descr, cog_descr in cog_summ[cog_id]:
            funcs.append(func)
            func_descrs.append(func_descr)
            cog_descrs.append(cog_descr)
        entry.append(cog_descrs.pop())
        entry.append("<br>".join(funcs))
        entry.append("<br>".join(func_descrs))
        cog_entries.append(entry)

    cog_header = ["COG", "Occurences", "Description", "Category", "Category description"]
    return {
        "cog_header": cog_header,
        "cog_entries": cog_entries
    }


def orthogroup(request, og):
    tokens = og.split("_")
    try:
        og_id = int(tokens[1])
    except:
        menu = True
        invalid_id = True
        return render(request, "chlamdb/og.html", my_locals(locals()))
    else:
        valid_id = True

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    optional2status = db.get_config_table()
    og_counts = db.get_og_count([og_id], search_on="orthogroup")
    if len(og_counts.index) == 0:
        valid_id = False
        return render(request, "chlamdb/og.html", my_locals(locals()))

    annotations = db.get_genes_from_og(orthogroups=[og_id],
            terms=["locus_tag", "gene", "product", "length"])

    hsh_organisms = db.get_organism(annotations.index.tolist())
    hsh_genes = format_lst(annotations["gene"].tolist())
    hsh_products = format_lst(annotations["product"].tolist())
    n_homologues = og_counts.loc[og_id].sum()

    gene_annotations = []
    for index, values in enumerate(hsh_genes.items()):
        gene, cnt = values
        if pd.isna(gene):
            gene = "-"
        gene_annotations.append([index+1, gene, cnt])

    product_annotations = []
    for index, values in enumerate(hsh_products.items()):
        product, cnt = values
        if pd.isna(product):
            product = "-"
        product_annotations.append([index+1, product, cnt])

    cog_ctx, kegg_ctx = {}, {}
    if optional2status.get("COG", False):
        cog_ctx = og_tab_get_cog_annot(db, annotations.index.tolist())

    if optional2status.get("KEGG", False):
        kegg_ctx = og_tab_get_kegg_annot(db, annotations.index.tolist())

    og_conserv_ctx = tab_og_conservation_tree(db, og_id)
    length_tab_ctx = tab_lengths(n_homologues, annotations)
    homolog_tab_ctx = tab_homologs(db, annotations, hsh_organisms)
    context = {
        "valid_id": valid_id,
        "optional2status": optional2status,
        "n_homologues": n_homologues,
        "og": og,
        "menu": True,
        "gene_annotations": gene_annotations,
        "product_annotations": product_annotations,
        **homolog_tab_ctx,
        **length_tab_ctx,
        **og_conserv_ctx,
        **cog_ctx,
        **kegg_ctx
    }
    return render(request, "chlamdb/og.html", my_locals(context))


def tab_general(seqid, hsh_organism, gene_loc, annot):
    organism = hsh_organism[seqid]
    strand, beg, end = gene_loc[seqid]
    gene = annot.loc[seqid].gene
    if pd.isna(gene):
        gene = "-"
    length = annot.loc[seqid].length
    product = annot.loc[seqid]["product"]
    locus_tag = annot.loc[seqid].locus_tag
    return {
        "locus_tag": locus_tag,
        "organism": organism,
        "strand": strand,
        "gene": gene,
        "start": beg,
        "end": end,
        "nucl_length": end-beg,
        "length": length,
        "prot": product
    }


# to be moved somewhere else at some point
def to_color_code(c):
    red = int(256*c.red)
    green = int(256*c.green)
    blue = int(256*c.blue)
    return f"#{red:x}{green:x}{blue:x}"


class LocusHeatmapColumn(SimpleColorColumn):
    def __init__(self, values, ref_taxon=None, header=None):
        super().__init__(values, header)
        self.ref_taxon = ref_taxon
        self.min_val = min(v for k, v in values.items())
        self.max_val = max(v for k, v in values.items())


    def get_face(self, index):
        index = int(index)
        if index==self.ref_taxon:
            text_face = TextFace(" - ")
            text_face.inner_background.color = EteTree.GREEN
            return text_face

        val = self.values.get(index, None)
        if val is None:
            return TextFace(" - ")

        color = colors.linearlyInterpolatedColor(colors.gray,
                colors.firebrick, self.min_val, self.max_val, val)
        text_face = TextFace(int(val))
        text_face.inner_background.color = to_color_code(color)
        return text_face


def get_sequence(db, seqid, flanking=0):
    loc      = db.get_gene_loc([seqid])
    bioentry = db.get_bioentry(from_val=seqid)
    seq      = db.get_DNA_sequence(bioentry)
    strand, start, stop = loc[seqid]
    start -= 1

    if start < 50:
        start_w_flank = 0
        red_start = start
    else:
        start_w_flank = start-flanking
        red_start = 50

    if stop+flanking > len(seq):
        stop_w_flank = len(seq)-1
    else:
        stop_w_flank = stop+flanking
    red_stop = red_start + stop-start
    fet = SeqFeature(FeatureLocation(start_w_flank, stop_w_flank, strand=strand))
    extracted = fet.extract(seq)
    return extracted[0:red_start] + "<font color='red'>" + \
            extracted[red_start:red_stop] + "</font>" + extracted[red_stop:]


def og_tab_get_pfam_annot(db, seqid):
    pfam_hits = db.get_pfam_hits_info(seqid)
    feature_viewer_fet = []
    pfam_grouped = pfam_hits.groupby(["pfam"])
    pfam_starts = pfam_grouped["start"].apply(list)
    pfam_ends = pfam_grouped["end"].apply(list)
    pfam_defs_df = db.get_pfam_def(pfam_hits.pfam.tolist())

    pfam_defs = []
    for pfam, starts in pfam_starts.iteritems():
        ends = pfam_ends.loc[pfam]
        name = format_pfam(pfam)
        data = "["+",".join(f"{{x:{start}, y:{end}}}" for start, end in zip(starts, ends))+"]"
        feature = (
            f"{{ data: {data}, "
            f" name: \"{name}\","
            "  color: \"#0F8292\","
            "  type : \"rect\","
            "}"
        )
        pfam_def = pfam_defs_df["def"].loc[pfam]
        pfam_defs.append((name, pfam_def))
        feature_viewer_fet.append(feature)

    return {"pfam_domains": "[" + ",".join(feature_viewer_fet) + "]",
            "pfam_def": pfam_defs}


def locusx(request, locus=None, menu=True):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    
    if locus==None:
        valid_id = False
        return render(request, 'chlamdb/locus.html', my_locals(locals()))

    try:
        seqid = db.get_seqid(locus_tag=locus)
    except:
        valid = False
        return render(request, 'chlamdb/locus.html', my_locals(locals()))
    else:
        valid_id = True

    optional2status = db.get_config_table()
    gene_loc = db.get_gene_loc([seqid])
    og_inf   = db.get_og_count([seqid], search_on="seqid")
    og_id    = int(og_inf.loc[seqid].orthogroup) # need to convert from numpy64 to int
    og_annot = db.get_genes_from_og(orthogroups=[og_id],
            terms=["locus_tag", "gene", "product", "length"])
    all_og_c = db.get_og_count([og_id], search_on="orthogroup")
    all_org  = db.get_organism(og_annot.index.tolist(), as_hash=True)
    n_homologues = all_og_c.loc[og_id].sum()
    translation = db.get_translation(seqid)
    sequence = get_sequence(db, seqid, flanking=50)

    homolog_tab_ctx = tab_homologs(db, og_annot, all_org, seqid, og_id)
    general_tab     = tab_general(seqid, all_org, gene_loc, og_annot)
    og_conserv_ctx  = tab_og_conservation_tree(db, og_id, compare_to=seqid)

    try:
        og_phylogeny_ctx = tab_og_phylogeny(db, og_id, og_annot)
    except:
        og_phylogeny_ctx = {}

    kegg_ctx, cog_ctx, pfam_ctx = {}, {}, {}
    if optional2status.get("KEGG", False):
        kegg_ctx = og_tab_get_kegg_annot(db, [seqid])

    if optional2status.get("COG", False):
        cog_ctx = og_tab_get_cog_annot(db, [seqid])

    if optional2status.get("pfam", False):
        pfam_ctx = og_tab_get_pfam_annot(db, [seqid])

    context = {
        "valid_id": valid_id,
        "optional2status": optional2status,
        "menu": True,
        "data": ["foo", "bar"],
        "n_homologues": n_homologues,
        "og_id": format_orthogroup(og_id, to_url=True),
        "translation": translation,
        "seq": sequence,
        **cog_ctx,
        **kegg_ctx,
        **homolog_tab_ctx,
        **general_tab,
        **og_conserv_ctx,
        **og_phylogeny_ctx,
        **pfam_ctx
    }
    return render(request, 'chlamdb/locus.html', my_locals(context))


def str_if_none(s):
    if s is None:
        return "-"
    return s


def search_bar_helper(fam_name, fam_url, results):
    val = getattr(results, fam_name)
    if val is None:
        val = "-"
    else:
        val = f"<a href=\"/{fam_url}/{val}\">{val}</a>"
    return val


def search_bar(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    option2status = db.get_config_table()

    index = sb.ChlamdbIndex.use_index(settings.SEARCH_INDEX)
    user_query = request.GET.get("accession")

    results = list(index.search(user_query, limit=100))

    if len(results) == 0:
        ctx = {"search_failed": True, "search_term": user_query}
        return render(request, "chlamdb/search.html", my_locals(ctx))

    search_results = []
    for result in results:
        locus_tag = format_locus(result.locus_tag, to_url=True)
        gene = str_if_none(result.gene)
        product = str_if_none(result.product)
        orthogroup = format_orthogroup(result.og, to_url=True, from_str=True)
        line = [locus_tag, gene, product, orthogroup]

        if option2status.get("COG", False):
            cog = search_bar_helper("cog", "fam_cog", result)
            line.append(cog)
        if option2status.get("KEGG", False):
            ko = search_bar_helper("ko", "fam_ko", result)
            line.append(ko)
        if option2status.get("pfam", False):
            pfam = search_bar_helper("pfam", "fam_pfam", result)
            line.append(pfam)
        line.append(result.organism)
        search_results.append(line)

    header = ["accession", "gene", "product", "Orthogroup", "organism"]
    insert_index = 4
    if option2status.get("COG", False):
        header.insert(insert_index, "COG")
        insert_index += 1
    if option2status.get("KEGG", False):
        header.insert(insert_index, "KO")
        insert_index += 1
    if option2status.get("pfam", False):
        header.insert(insert_index, "PFAM")

    ctx = {"search_term": user_query,
            "header": header,
            "search_results": search_results }
    return render(request, "chlamdb/search.html", my_locals(ctx))


def locusx_legacy(request, locus=None, menu=True):
    
    biodb = settings.BIODB
    print ('-- locus or search term: %s -- biodb %s' % (locus, biodb))
    print(request.method)
    if request.method == 'GET': 

        server, db = manipulate_biosqldb.load_db(biodb)

        if locus == None:
            menu = True
            locus = request.GET.get('accession').strip()

        if locus[-1] == '+':
            locus = re.sub("\+","",locus)
        elif "%" in locus:
            locus = re.sub("%", " ", locus)
        
        # check if Chlamydia trachomatis locus
        m = re.match("CT([0-9]+)", locus)
        if m:
            locus = "CT_%s" % (m.group(1))
            print("updated_search", locus)
        
        try:
            locus = int(locus)
            sql = 'select locus_tag from custom_tables_locus2seqfeature_id where seqfeature_id=%s' % (locus)
            locus = server.adaptor.execute_and_fetchall(sql,)[0][0]
        except:
            pass
        
        sql0 = 'select locus_tag from custom_tables_seqfeature_id2old_locus_tag t1 ' \
               ' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id where old_locus_tag="%s" ' % (locus)

        try:
            data = server.adaptor.execute_and_fetchall(sql0, )[0][0]
            old_locus_tag = locus
            locus = data
            input_type = 'locus_tag'
        except IndexError:
            sql0 = 'select locus_tag from orthology_detail where (locus_tag="%s" or protein_id="%s")' % (locus,
                                                                                                         locus)
        try:
            locus = server.adaptor.execute_and_fetchall(sql0, )[0][0]
            input_type = 'locus_tag'
        except IndexError:

            sql1 = 'select orthogroup from orthology_detail where orthogroup="%s"' % (locus)
            try:
                locus = server.adaptor.execute_and_fetchall(sql1, )[0][0]
                input_type = 'orthogroup'
            except IndexError:
                print ('not a valid id, trying search')
                return search(request)
        
        valid_id = True

        # retrieve basic data about locus
        columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
        sql2 = 'select %s from orthology_detail where %s="%s"' % (columns, input_type, locus)
        data = list(server.adaptor.execute_and_fetchall(sql2, )[0])
        orthogroup = data[0] 
        sql_old = 'select old_locus_tag from custom_tables_seqfeature_id2old_locus_tag t1 ' \
                  'inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id where locus_tag="%s" ' % (data[1])
        try:
            data_old = server.adaptor.execute_and_fetchall(sql_old, )[0][0]
            old_locus_tag = data_old
        except:
            pass
        # check if orthogroup
        m = re.match("group_([0-9]+)", locus)
        if m:
            input_type = 'orthogroup'
        else:
            input_type = 'locus_tag'

        if input_type == 'locus_tag':
            from chlamdb.plots import uniprot_feature_viewer

            print("input type locus tag ----------------")

            sql4 = 'select accession from orthology_detail where locus_tag="%s" limit 1' % (locus)
            genome_accession = server.adaptor.execute_and_fetchall(sql4,)[0][0]

            sql3 = 'select COG_name,code,t2.description,t1.query_start,t1.query_end, t1.hit_start, t1.hit_end, t1.query_coverage, t1.hit_coverage, t1.identity, t1.evalue, t1.bitscore ' \
                   ' from COG_seqfeature_id2best_COG_hit t1 inner join COG_cog_names_2014 t2 on t1.hit_cog_id=t2.COG_id ' \
                   ' inner join COG_cog_id2cog_category t3 on t2.COG_id=t3.COG_id inner ' \
                   ' join COG_code2category t4 on t3.category_id=t4.category_id inner join annotation_seqfeature_id2locus t6 on t1.seqfeature_id=t6.seqfeature_id ' \
                   ' where locus_tag="%s";' % (locus)

            sql4 = 'select A.analysis_name,A.signature_accession, A.signature_description, start, stop,score, B.name,B.description from ' \
                   ' (select t1.*,t2.signature_accession,t2.signature_description,t2.interpro_id,t3.* from interpro_interpro t1 ' \
                   ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
                   ' inner join interpro_analysis t3 on t2.analysis_id=t3.analysis_id ' \
                   ' inner join annotation_seqfeature_id2locus t4 on t1.seqfeature_id=t4.seqfeature_id where locus_tag="%s") A left join ' \
                   ' interpro_entry B on A.interpro_id=B.interpro_id;' % (locus)

            sql5 = 'select t3.ko_accession, t3.name, t3.definition, t3.pathways, t3.modules, t2.thrshld, t2.score, t2.evalue from ' \
                   ' custom_tables_locus2seqfeature_id t1 ' \
                   ' inner join enzyme_seqfeature_id2ko t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                   ' inner join enzyme_ko_annotation t3 on t2.ko_id=t3.ko_id where t1.locus_tag="%s";' % (locus)

            sql6 = 'select uniprot_id from locus_tag2uniprot_hit where locus_tag="%s";' % (locus)

            sql7 = 'select A.ko_id,pathway_name,pathway_category,description from (select ko_id from enzyme_locus2ko' \
                   ' where locus_tag="%s") A inner join enzyme_pathway2ko_v1 B on A.ko_id=B.ko_id ' \
                   ' inner join enzyme_kegg_pathway C on B.pathway_id=C.pathway_id where pathway_category !="1.0 Global and overview maps";' % (locus)

            sql8 = 'select t4.module_name, t4.module_sub_sub_cat, t4.description from custom_tables_locus2seqfeature_id t1  ' \
                   ' inner join enzyme_seqfeature_id2ko t2 on t1.seqfeature_id=t2.seqfeature_id' \
                   ' inner join enzyme_module2ko t3 on t2.ko_id=t3.ko_id ' \
                   ' inner join enzyme_kegg_module t4 on t3.module_id=t4.module_id where t1.locus_tag="%s";' % (locus)

            sql9 = 'select mol_weight,isoelectric_point,aromaticity,instability_index,fraction_helix,fraction_turn,' \
                   ' fraction_sheet from custom_tables_locus2pepstats where locus_tag="%s";' % (locus)

            sql10 = 'select operon_id from custom_tables_locus2seqfeature_id t1 ' \
                    ' inner join custom_tables_DOOR2_operons t2 on t1.seqfeature_id=t2.seqfeature_id' \
                    ' where t1.locus_tag="%s"' % (locus)

            sql11 = 'select db_xref_name,db_accession from custom_tables_locus2seqfeature_id as t1 ' \
                    ' inner join custom_tables_uniprot_id2seqfeature_id as t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                    ' inner join custom_tables_uniprot_db_xref as t3 on t2.uniprot_id=t3.uniprot_id ' \
                    ' inner join custom_tables_db_xref as t4 on t3.db_xref_id=t4.db_xref_id ' \
                    ' where locus_tag="%s" and db_xref_name not in ("GO","InterPro", "Pfam");' % (locus)

            sql12 = 'select uniprot_status,annotation_score,gene,recommendedName_fullName,comment_function,ec_number,' \
                    ' comment_similarity,comment_subunit,comment_catalyticactivity,proteinExistence,uniprot_accession ' \
                    ' from custom_tables_locus2seqfeature_id as t1 inner join custom_tables_uniprot_id2seqfeature_id as t2 ' \
                    ' on t1.seqfeature_id=t2.seqfeature_id inner join custom_tables_uniprot_annotation as t3 ' \
                    ' on t2.seqfeature_id=t3.seqfeature_id where locus_tag="%s";' % (locus)


            sql13 = 'select go_term_id, term_type, go_description from (select go_term_id, go_description ' \
                    ' from custom_tables_locus2seqfeature_id as t1  ' \
                    'inner join custom_tables_uniprot_go_terms as t2 on t1.seqfeature_id=t2.seqfeature_id  ' \
                    ' where t1.locus_tag="%s") A inner join gene_ontology.term as B on A.go_term_id=B.acc;' % (locus)

            sql15 = 'select count(*) from custom_tables_locus2seqfeature_id t1 inner join blastnr_blastnr t2' \
              ' on t1.seqfeature_id=t2.seqfeature_id where locus_tag="%s";' % (locus)

            sql16 = 'select count(*) from custom_tables_locus2seqfeature_id t1 ' \
              ' inner join blastnr_blast_swissprot t2 on t1.seqfeature_id=t2.seqfeature_id' \
              ' where locus_tag="%s";' % (locus)

            sql17 = 'select phylogeny from phylogenies_BBH where orthogroup="%s"' % (orthogroup)

            sql18 = 'select signature_accession,start,stop from interpro where analysis="Phobius" and locus_tag="%s" ' \
                    ' and signature_accession in ("TRANSMEMBRANE",' \
                    ' "SIGNAL_PEPTIDE_C_REGION","SIGNAL_PEPTIDE_N_REGION", "SIGNAL_PEPTIDE", "SIGNAL_PEPTIDE_H_REGION");' % (locus)

            sql19 = 'select signature_accession, interpro_description,start, stop from interpro ' \
                    ' where analysis="SUPERFAMILY" and locus_tag="%s";' % (locus)

            sql20 = 'select t8.uniprot_accession,t2.evalue, t2.bitscore_first_hsp, t2.identity, t2.query_TMS, t2.hit_TMS, ' \
                    ' t2.query_cov, t2.hit_cov,t4.tc_name as transporter_name, t4.description as transporter_description, ' \
                    ' t5.tc_name as superfamily, t5.description as superfamily_description, ' \
                    ' t6.tc_name as family_name, t6.description as family_description, t7.tc_name as subfamily_name, ' \
                    ' t7.description as subfamily_description, t8.tcdb_description, t8.organism  from custom_tables_locus2seqfeature_id t1 ' \
                    ' inner join transporters_transporters t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                    ' inner join transporters.transporter_table t3 on t2.transporter_id=t3.transporter_id ' \
                    ' inner join transporters.tc_table t4 on t3.transporter_id=t4.tc_id ' \
                    ' inner join transporters.tc_table t5 on t3.superfamily=t5.tc_id ' \
                    ' inner join transporters.tc_table t6 on t3.family=t6.tc_id ' \
                    ' inner join transporters.tc_table t7 on t3.subfamily=t7.tc_id ' \
                    ' inner join transporters.uniprot_table t8 on t2.hit_uniprot_id=t8.uniprot_id ' \
                    ' where t1.locus_tag="%s";' % (locus)

            sql21 = 'select seqfeature_id, taxon_id from custom_tables_locus2seqfeature_id where locus_tag="%s"' % (locus)

            seqfeature_data = server.adaptor.execute_and_fetchall(sql21,)[0]
            taxon_id = seqfeature_data[1]
            seqfeature_id = seqfeature_data[0]

            sql22 = 'select temporal_class, EB_proteome, eggNOG, comment, hpi_2_1,hpi_2_2, ' \
                    ' hpi_2_3, hpi_48_1,hpi_48_2,hpi_48_3,hpi_96_1,hpi_96_2,hpi_96_3, extracellular_1, ' \
                    ' extracellular_2, extracellular_3 from rnaseq_%s where seqfeature_id=%s' % (taxon_id, seqfeature_id)

            sql23 = 'select pmid from string_seqfeature_id2pmid where seqfeature_id=%s;' % (seqfeature_id)

            sql24 = 'select hash from annotation_hash2seqfeature_id t1 ' \
                    ' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id where t2.locus_tag="%s";' % (locus)
            try:
                interpro_hash = server.adaptor.execute_and_fetchall(sql24,)[0][0]
            except:
                interpro_hash = None

            sql25 = f'select distinct t2.effectors,SVM_value from annotation_seqfeature_id2locus t1' \
                    f' left join effectors_predicted_BPBAac t2 on t1.seqfeature_id=t2.seqfeature_id where locus_tag="{locus}";'
            sql26 = f'select distinct effectors from annotation_seqfeature_id2locus t1' \
                    f' left join effectors_predicted_DeepT3 t2 on t1.seqfeature_id=t2.seqfeature_id where locus_tag="{locus}";'
            sql27 = f'select distinct effectors,probability from annotation_seqfeature_id2locus t1' \
                    f' left join effectors_predicted_T3MM t2 on t1.seqfeature_id=t2.seqfeature_id where locus_tag="{locus}";'
            sql28 = f'select distinct t2.effectors,score from annotation_seqfeature_id2locus t1' \
                    f' left join effectors_predicted_effectiveT3 t2 on t1.seqfeature_id=t2.seqfeature_id where locus_tag="{locus}";'

            sql29 = f'select final_prediction,score from annotation_seqfeature_id2locus t1' \
                    f' inner join custom_tables_seqfeature_id2psortb t2 on t1.seqfeature_id=t2.seqfeature_id where locus_tag="{locus}";'

            sql30 = f'select paperblast_id from string.seqfeature_id2paperblast where seqfeature_id={seqfeature_id};'

            sql31 = f'select * from custom_tables_seqfeature_id2pdb_BBH where seqfeature_id={seqfeature_id};'

            sql32 = f'select distinct ec,value from enzyme_seqfeature_id2ec t1 ' \
                    f' inner join enzyme_enzymes t2 on t1.ec_id=t2.enzyme_id ' \
                    f' inner join enzyme_enzymes_dat t3 on t1.ec_id=t3.enzyme_dat_id ' \
                    f' where seqfeature_id="{seqfeature_id}" and line="description";'

            sql33 = f'select EC,definition,thrshld,score,evalue,t2.ko_accession from enzyme_seqfeature_id2ko t1 ' \
                    f' inner join enzyme_ko_annotation t2 on t1.ko_id=t2.ko_id ' \
                    f' inner join enzyme_ko2ec t3 on t2.ko_accession=t3.ko_id where seqfeature_id={seqfeature_id}; '
            
            try:
                ec_priam_data = server.adaptor.execute_and_fetchall(sql32,)
            except:
                ec_priam_data = False

            try:
                ec_ko_data = server.adaptor.execute_and_fetchall(sql33,)
            except:
                ec_ko_data = False

            try:
                pdb_data = server.adaptor.execute_and_fetchall(sql31,)[0]
            except:
                pdb_data = False

            try:
                psort_data = server.adaptor.execute_and_fetchall(sql29,)[0]
            except:
                psort_data = False
            
            try:
                BPBAac_data = list(server.adaptor.execute_and_fetchall(sql25,)[0])
                DeepT3_data = list(server.adaptor.execute_and_fetchall(sql26,)[0])
                T3MM_data = list(server.adaptor.execute_and_fetchall(sql27,)[0])
                effectiveT3_data = list(server.adaptor.execute_and_fetchall(sql28,)[0])

                effector_sum = 0
                if BPBAac_data[0] == 1:
                    effectors = True
                    effector_sum+=1
                if DeepT3_data[0] == 1:
                    effectors= True
                    effector_sum+=1
                if T3MM_data[0] == 1:
                    effectors = True
                    effector_sum+=1
                if effectiveT3_data[0] == 1:
                    effectors = True 
                    effector_sum+=1            
            except:
                effector_sum = 0
                effectors = False
            pubmed_count = 0
            try:
                string_data = server.adaptor.execute_and_fetchall(sql23,)
                pubmed_count+=len(string_data)
            except:
                string_data = False
            try:
                paperblast_data = server.adaptor.execute_and_fetchall(sql30,)
                pubmed_count+=len(paperblast_data)
            except:
                paperblast_data = False

            try:
                rnaseq_data = server.adaptor.execute_and_fetchall(sql22,)[0]
            except:
                rnaseq_data = False

            try:
                transporter_data = [str(i) for i in server.adaptor.execute_and_fetchall(sql20, )[0]]
                transporter_data[16]= ' '.join(transporter_data[16].split(' ')[1:]).split("OS=")[0]
                transporter_data[17] = transporter_data[17].split("(")[0]
                # remove species name in case already in description
                transporter_data[16] = re.sub(transporter_data[17],"",transporter_data[16])

            except:
                transporter_data = False

            try:

                superfamily_data =  server.adaptor.execute_and_fetchall(sql19, )

                superfamily_features = uniprot_feature_viewer.superfamily_data2features_string(superfamily_data)


            except:

                superfamily_features = ""

            try:

                phobius_data = server.adaptor.execute_and_fetchall(sql18, )

                phobius_features = uniprot_feature_viewer.phobius_data2features_string(phobius_data)

            except:

                phobius_features = ""
            try:
                trial = server.adaptor.execute_and_fetchall(sql17, )[0][0]
                closest_phylo = True
            except:
                closest_phylo = False

            try:
                n_blastnr_hits = server.adaptor.execute_and_fetchall(sql15, )[0][0]
            except:
                n_blastnr_hits = 0
            try:
                n_swissprot_hits = server.adaptor.execute_and_fetchall(sql16, )[0][0]
            except:
                n_swissprot_hits = 0
            try:
                uniprot_go_terms = server.adaptor.execute_and_fetchall(sql13, )
            except:
                uniprot_go_terms = False
            if not uniprot_go_terms or len(uniprot_go_terms) == 0:
                try:
                    # go terms from interpro data
                    sql14 = 'select t4.acc,t4.term_type,t4.name from (select interpro_accession from interpro ' \
                          ' where locus_tag="%s" and interpro_accession != "0" group by interpro_accession) t1 ' \
                          ' inner join interpro_entry as t2 on t1.interpro_accession=t2.name ' \
                          ' inner join interpro_interpro2gene_ontology as t3 on t2.interpro_id=t3.interpro_id ' \
                          ' inner join gene_ontology.term as t4 on t3.go_id=t4.id group by acc;' % (locus)

                    interpro_go_terms = server.adaptor.execute_and_fetchall(sql14, )
                except:

                    uniprot_go_terms = False

            try:
                uniprot_annotation = list(server.adaptor.execute_and_fetchall(sql12, )[0])
                if uniprot_annotation[0] == '1':
                    uniprot_annotation[0] = 'Reviewed'
                if uniprot_annotation[0] == '0':
                    uniprot_annotation[0] = 'Unreviewed'  

            except:
                uniprot_annotation = False
            try:
                dbxref_uniprot_data = server.adaptor.execute_and_fetchall(sql11, )
                for i in dbxref_uniprot_data:
                    if i[0] == 'Swiss-Prot' and not '_' in i[1]:
                        uniprot_accession = i[1]
            except:
                dbxref_uniprot_data = False

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
                           ' inner join custom_tables_locus2seqfeature_id C on B.seqfeature_id=C.seqfeature_id' % (locus)



                    operon_ofs = server.adaptor.execute_and_fetchall(sqlo, )
                    operon_locus = [i[0] for i in operon_ofs]

                    locus2annot, \
                    locus_tag2cog_catego, \
                    locus_tag2cog_name, \
                    locus_tag2ko, \
                    pathway2category, \
                    module2category, \
                    ko2ko_pathways, \
                    ko2ko_modules,\
                    locus2interpro = get_locus_annotations(biodb, operon_locus)

                    operon = False


                except:
                    operon = False
                    operon_ofs = False

            #operon=False
            temp_location = os.path.join(settings.BASE_DIR, "assets/temp/")
            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
            name = 'temp/' + os.path.basename(temp_file.name)

            if operon or operon_ofs:
                if operon_locus[0] == '-':
                    lst = [i[3] for i in operon]
                else:
                    lst = operon_locus
            else:
                lst = []

            locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, 
                                                                                            db, 
                                                                                            biodb, 
                                                                                            [locus],
                                                                                            temp_file.name, 
                                                                                            15000,
                                                                                            cache, 
                                                                                            color_locus_list=lst)


            try:
                cog_data = server.adaptor.execute_and_fetchall(sql3, )[0]
            # TODO: distinguish keyerror from missing table
            except:
                cog_data = False

            try:
                interpro_data_detail = server.adaptor.execute_and_fetchall(sql4, )
            # TODO: distinguish keyerror from missing table
            # TODO if missing table, hide interpro tab and domain scheme
            except:
                interpro_data_detail= []

            try:
                ko_data = server.adaptor.execute_and_fetchall(sql5, )[0]
            except:
                ko_data= False
            # TODO: distinguish keyerror from missing table
            try:
                ko_pathway_data = server.adaptor.execute_and_fetchall(sql7, )
            except:
                ko_pathway_data = False
            try:
                ko_module_data = server.adaptor.execute_and_fetchall(sql8, )
            except:
                ko_module_data = False
            try:
                uniprot_id = server.adaptor.execute_and_fetchall(sql6, )[0][0]
            except:
                uniprot_id = False
            try:
                protparams_data = list(server.adaptor.execute_and_fetchall(sql9,)[0])
                protparams_data[0] = round(float(protparams_data[0])/1000, 2)
                protparams_data[1] = round(protparams_data[1], 2)
                protparams_data[2] = round(protparams_data[2], 2)
                protparams_data[3] = round(protparams_data[3], 2)
                protparams_data[4] = round(protparams_data[4]*100, 2)
                protparams_data[5] = round(protparams_data[5]*100, 2)
                protparams_data[6] = round(protparams_data[6]*100, 2)
            except:
                protparams_data = False
            try:
                sql_interpro = 'select interpro_accession, interpro_description from interpro' \
                               ' where locus_tag="%s" and interpro_accession !="0"' \
                               ' group by interpro_accession;' % (locus)
                interpro_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_interpro, )]

            except:
                interpro_data = False

            if interpro_data:
                interpro2taxononmy = {}
                for one_entry in interpro_data:
                    sql = 'select p_bacteria,p_eukaryote,p_archae,p_virus, bacteria,eukaryote,archae,virus ' \
                          ' from interpro_entry t1 inner join interpro_interpro_taxonomy_v_60 t2 on t1.interpro_id=t2.interpro_id ' \
                          ' where name="%s";' % one_entry[0]
                    try:
                        interpro2taxononmy[one_entry[0]] = server.adaptor.execute_and_fetchall(sql,)[0]
                    except:
                        pass


            try:
                sql_pfam = 'select signature_accession, signature_description,start,stop' \
                           ' from interpro where locus_tag="%s" ' \
                           ' and analysis="Pfam";' % (locus)
                pfam_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pfam, )]
                from chlamdb.plots import uniprot_feature_viewer
                features_js = uniprot_feature_viewer.intero_data2features_string(pfam_data)

            except:
                pfam_data = False

            try:
                sql_pathway = 'select locus_tag, pathways, interpro_description from interpro where ' \
                              ' locus_tag="%s" and pathways!="0" group by pathways;'  % (locus)

                pathway_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pathway, )]
                all_path = {}
                for one_path in pathway_data:
                    data_path = one_path[1].split('|')
                    for one_db in data_path:
                        one_db_info = one_db.split(':')
                        db = one_db_info[0]
                        if db not in all_path:
                            all_path[db] = [one_db_info[1][1:]]
                        else:
                            if one_db_info[1][1:] not in all_path[db]:
                                all_path[db].append(one_db_info[1][1:])

            except:
                pathways_data = False

            seq_start = int(data[3])
            seq_end = int(data[4])
            strand = int(data[5])
            leng = (seq_end-seq_start)+100

            nucl_length = seq_end-seq_start+1
            aa_length = nucl_length/3

            seq = manipulate_biosqldb.location2sequence(server, genome_accession, biodb, seq_start-50, leng)

            if strand == -1:

                from Bio.Seq import Seq
                seq_obj = Seq(seq)
                seq = str(seq_obj.reverse_complement())
                seq = seq[0:49] + '<font color="red">' + seq[49:-50] + '</font>' + seq[-50:len(seq)]
            else:
                seq = seq[0:50] + '<font color="red">' + seq[50:-50] + '</font>' + seq[-50:len(seq)]


        if input_type == 'orthogroup':
            
            
            # consensus annotation
            sql_group1 = 'select `rank`,count,description from orthology_orthogroup2gene t1 ' \
                         ' inner join orthology_orthogroup t2 on t1.group_id=t2.orthogroup_id where t2.orthogroup_name="%s";' % (locus)
                         
            sql_group2 = 'select `rank`,count,description from orthology_orthogroup2product t1 ' \
                         ' inner join orthology_orthogroup t2 on t1.group_id=t2.orthogroup_id where t2.orthogroup_name="%s";' % (locus)
                         
            sql_group3 = 'select `rank`, COG_name, t3.description, count, code, t5.description from orthology_orthogroup2cog t1 ' \
                         ' inner join orthology_orthogroup t2 on t1.group_id=t2.orthogroup_id ' \
                         ' inner join COG_cog_names_2014 t3 on t1.COG_id=t3.COG_id ' \
                         ' inner join COG_cog_id2cog_category t4 on t3.COG_id=t4.COG_id' \
                         ' inner join COG_code2category t5 on t4.category_id=t5.category_id where t2.orthogroup_name="%s";' % (locus)
                         
            sql_group4 = 'select `rank`,ko_accession,count,name,definition,EC,pathways,modules from orthology_orthogroup2ko t1 ' \
                         ' inner join orthology_orthogroup t2 on t1.group_id=t2.orthogroup_id ' \
                         ' inner join enzyme_ko_annotation t3 on t1.ko_id=t3.ko_id where t2.orthogroup_name="%s";' % (locus)
                         
            sql_group5 = f'select distinct t1.seqfeature_id,start,stop,signature_accession,signature_description from interpro_interpro t1 ' \
                         f' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
                         f' inner join interpro_analysis t3 on t2.analysis_id=t3.analysis_id ' \
                         f' inner join orthology_seqfeature_id2orthogroup t4 on t1.seqfeature_id=t4.seqfeature_id ' \
                         f' inner join orthology_orthogroup t5 on t4.orthogroup_id=t5.orthogroup_id ' \
                         f' where t5.orthogroup_name="{locus}" and analysis_name="Pfam" order by start;'
            
            if db_driver == 'mysql':
                sql_group6 = f'select char_length(translation) as len from orthology_detail where orthogroup="{locus}"'
            if db_driver == 'sqlite':
                sql_group6 = f'select length(translation) as len from orthology_detail where orthogroup="{locus}"'
        
            # protein length distribution
            sql_group7 = f'select uniprot_accession,uniprot_status,annotation_score,gene,recommendedName_fullName from orthology_orthogroup t1 ' \
                         f' inner join orthology_seqfeature_id2orthogroup t2 on t1.orthogroup_id=t2.orthogroup_id ' \
                         f' left join custom_tables_uniprot_id2seqfeature_id t3 on t2.seqfeature_id=t3.seqfeature_id ' \
                         f' left join custom_tables_uniprot_annotation t4 on t2.seqfeature_id=t4.seqfeature_id ' \
                         f' where orthogroup_name="{locus}";'
            
            
            sql_group8 = 'select phylogeny from phylogenies_BBH where orthogroup="%s"' % (locus)
            # TM domains 
            
            # Signal Peptide
            try:
                trial = server.adaptor.execute_and_fetchall(sql_group8, )[0][0]
                closest_phylo = True
            except:
                closest_phylo = False
            
            gene_annotations = server.adaptor.execute_and_fetchall(sql_group1,)
            product_annotations = server.adaptor.execute_and_fetchall(sql_group2,)
            try:
                COG_annotations = server.adaptor.execute_and_fetchall(sql_group3,)
            except:
                COG_annotations = []
            server.adaptor.execute_and_fetchall(sql_group4,)
            try:
                KO_annotations = [list(i) for i in server.adaptor.execute_and_fetchall(sql_group4,)]
            except:
                KO_annotations = []
            try:
                pfam_annotations = server.adaptor.execute_and_fetchall(sql_group5,)
            except:
                pfam_annotations = []
            protein_lengths = [int(i[0]) for i in server.adaptor.execute_and_fetchall(sql_group6,)]
            try:
                uniprot_annotations = server.adaptor.execute_and_fetchall(sql_group7,)
                unreviewed = len([i for i in uniprot_annotations if i[1] == "unreviewed"])
                reviewed = [i for i in uniprot_annotations if i[1] == "reviewed"]
                mapped = unreviewed + len(reviewed)
                unmapped = len([i for i in uniprot_annotations if i[1] == None])
            except:
                uniprot_annotations = False

            for row in KO_annotations:
                row[6] = row[6].replace("ko", "map")
                row[6] = row[6].split(",")
                row[7] = row[7].split(",")     
            
            if len(pfam_annotations) > 0:
                domain_list = []
                domain2description = {}
                seqfeature_id2pfam_domains = {}
                for row in pfam_annotations:
                    if row[3] not in domain_list:
                        domain_list.append(row[3])
                        domain2description[row[3]] = row[4]
                    if row[0] not in seqfeature_id2pfam_domains:
                        seqfeature_id2pfam_domains[row[0]] = [row[3]]
                    else:
                        seqfeature_id2pfam_domains[row[0]].append(row[3])
                domain2color = {}
                color_list = sns.color_palette("hls", len(domain_list)).as_hex()
                for n, domain in enumerate(domain_list):
                    domain2color[domain]=color_list[n]
                l = list(tuple(i) for i in seqfeature_id2pfam_domains.values())
                s = set(l)
                topo2count = dict([(x,l.count(x)) for x in set(l)])
                longest_topo = len(max(list(topo2count.keys()), key=len))

            import plotly.graph_objects as go
            from collections import Counter
            import plotly.figure_factory as ff
            
            if len(protein_lengths) >2:
                length_distrib = True
                mean_protein_length = round(numpy.mean(protein_lengths),2)
                std_protein_length = round(numpy.std(protein_lengths),2)
                min_protein_length = min(protein_lengths)
                max_protein_length = max(protein_lengths)
                median_protein_length = round(numpy.median(protein_lengths),2)
                fig1 = ff.create_distplot([protein_lengths], ["sequence length"], bin_size=20)
                fig1.update_xaxes(range=[0, max(protein_lengths)])
                fig1.layout.margin.update({"l": 80,
                "r": 20,
                "b": 40,
                "t": 20,
                "pad": 10,
                })
                html_plot_prot_length = manipulate_biosqldb.make_div(fig1, div_id="distplot")
            else:
                length_distrib = False
            
            
            # TM statistics
            # case when no interpro data loaded
            try:
                sql = 'select TM from orthology_detail where orthogroup="%s";' % (locus)
                TM_counts = Counter([int(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)])
                
                for i in range(max(TM_counts.keys())):
                    if i not in TM_counts:
                        TM_counts[i] = 0
                
                        
                plot_data = [go.Bar(
                            x= ["%s TM" % i for i in range(0, len(TM_counts.keys() ) )], #["%s TM" % i for i in TM_counts.keys()],
                            y=[TM_counts[i] for i in range(0, len(TM_counts.keys() ) )]
                    )]

                layout = go.Layout(
                    title='',
                    yaxis=go.layout.YAxis(
                        title=go.layout.yaxis.Title(
                            text='Number of occurences',
                            font=dict(
                                family='Courier New, monospace',
                                size=15,
                                color='#7f7f7f'
                            )
                        )
                    ),
                    xaxis=go.layout.XAxis(
                        title=go.layout.xaxis.Title(
                            text='Number of TM',
                            font=dict(
                                family='Courier New, monospace',
                                size=15,
                                color='#7f7f7f'
                            )
                        )
                    )
                )
                plot_width = 120 + (max(TM_counts.keys()) ) * 50
                fig = go.Figure(data=plot_data, 
                                layout=layout)

                fig.layout.margin.update({"l": 80,
                                        "r": 20,
                                        "b": 40,
                                        "t": 20,
                                        "pad": 10,
                                        })

                html_plot = manipulate_biosqldb.make_div(fig, div_id="barplot")
            except:
                html_plot = ''
                

        if data[2] == '-':
            data[2] = data[1]

        orthogroup = data[0]

        fasta = "%s_fasta/%s.txt" % (biodb, orthogroup)
        alignment = "%s_fasta/%s.html" % (biodb, orthogroup)
        alignment_fasta = "%s_fasta/%s.faa" % (biodb, orthogroup)
        alignment_fasta_nucl = "%s_fasta_nucl/%s_nucl.txt" % (biodb, orthogroup)
        tree_unrooted = "%s_fasta/%s_tree.svg" % (biodb, orthogroup)
        tree_rooted = "%s_fasta/%s_tree_reroot.svg" % (biodb, orthogroup)
        tree_file = "%s_fasta/%s.phy_phyml_tree.txt" % (biodb, orthogroup)


        if not os.path.exists(settings.BASE_DIR + '/assets/' + alignment):
            align_path = False
        else:
            align_path = True

        #columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
        #      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
        #sql3 = 'select %s from orthology_detail where orthogroup = "%s" ' % (columns, biodb, orthogroup)

        #homologues = list(server.adaptor.execute_and_fetchall(sql3, ))
        sql_groups = 'select count(*) from orthology_seqfeature_id2orthogroup t1 ' \
                     ' inner join orthology_orthogroup t2 on t1.orthogroup_id=t2.orthogroup_id where t2.orthogroup_name="%s";' % (orthogroup)

        homologues = server.adaptor.execute_and_fetchall(sql_groups, )[0][0]

        # check if one of the homolog has TM(s) domains
        sql_TM_SP = 'select count(*) from orthology_seqfeature_id2orthogroup t1 ' \
              ' inner join orthology_orthogroup t2 on t1.orthogroup_id=t2.orthogroup_id ' \
              ' inner join interpro_interpro t3 on t1.seqfeature_id=t3.seqfeature_id' \
              ' inner join interpro_signature t4 on t3.signature_id=t4.signature_id ' \
              ' where signature_accession in ("TRANSMEMBRANE", "SIGNAL_PEPTIDE_C_REGION", "SIGNAL_PEPTIDE", "SIGNAL_PEPTIDE_N_REGION") ' \
              ' and t2.orthogroup_name="%s" ; ' % (orthogroup)
        try:
            tm_count = server.adaptor.execute_and_fetchall(sql_TM_SP, )[0][0]
        except:
            tm_count = 0
        if tm_count > 0:
            show_tm_tree = True

        if int(homologues) >1:
            orthologs = True
        else:
            orthologs = False

        sql = 'select t1.locus_tag, t1.annotation from manual_annotation as t1 ' \
              ' inner join orthology_detail as t2 on t1.locus_tag=t2.locus_tag where orthogroup="%s";' % (data[0])
        
        try:
            cmt = server.adaptor.execute_and_fetchall(sql,)
        except:
            cmt = ()

        cmt_format = []
        for i in cmt:
            cmt_format.append([i[0], i[1].replace('\n', '<br />')])

        home_dir = os.path.dirname(os.path.realpath(__file__))
        local_file = "/../assets/%s/interpro/%s.html" % (biodb, data[2])
        interpro_check = home_dir + local_file
        if os.path.isfile(interpro_check):
            interpro_protein = True
        else:
            interpro_protein = False

        envoi = True

    print(my_locals(locals())["optional2status"])
    return render(request, 'chlamdb/locus.html', my_locals(locals()))


def hydropathy(request, locus):
    biodb = settings.BIODB
    print ('hydropathy -- %s --%s' % (biodb, locus))

    from chlamdb.plots import hydrophobicity_plots

    fig = hydrophobicity_plots.locus2hydrophobicity_plot(biodb, locus)

    path = settings.BASE_DIR + '/assets/temp/hydro.png'
    fig.savefig(path,dpi=500)
    asset_path = '/temp/hydro.png'

    return render(request, 'chlamdb/hydropathy.html', my_locals(locals()))



def aa_comp_locus(request, locus_tag):
    biodb = settings.BIODB
    from chlamdb.plots import pca_seq_composition
    import numpy
    server, db = manipulate_biosqldb.load_db("%s" % biodb)

    sql1 = 'select t2.taxon_id, t2.seqfeature_id from custom_tables_locus2seqfeature_id t1 ' \
           ' inner join custom_tables_aa_usage_count t2 ' \
           ' on t1.seqfeature_id=t2.seqfeature_id where t1.locus_tag="%s"' % (locus_tag)

    data = server.adaptor.execute_and_fetchall(sql1,)[0]
    taxon_id = data[0]
    seqfeature_id = data[1]

    sql2 = 'select * from custom_tables_aa_usage_count where taxon_id=%s' % (taxon_id)
    data = [i for i in server.adaptor.execute_and_fetchall(sql2,)]

    for n, row in enumerate(data):
        if row[1] == seqfeature_id:
            target_feature = n + 1
    mat = numpy.array([list(i)[3:23] for i in data])
    path = settings.BASE_DIR + '/assets/temp/hydro.png'
    asset_path = '/temp/hydro.png'
    pca_seq_composition.aa_composition_pca(mat, target_feature, path)


    return render(request, 'chlamdb/aa_pca_locus.html', my_locals(locals()))



def rnaseq_class(request, temporal_class, taxon_id):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db("%s" % biodb)

    sql2 = 'select t1.*, t2.locus_tag from rnaseq_%s t1  left join custom_tables_locus2seqfeature_id t2 ' \
           ' on t1.seqfeature_id=t2.seqfeature_id where temporal_class="%s"' % (taxon_id, temporal_class)

    data = server.adaptor.execute_and_fetchall(sql2,)

    return render(request, 'chlamdb/rnaseq_temporal_class.html', my_locals(locals()))



def gc_locus(request, locus_tag):
    biodb = settings.BIODB
    from chlamdb.plots import pairwiseid_plots
    import numpy

    print ('gc locus -- %s -- %s' % (biodb, locus_tag))

    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = 'select t2.taxon_id, t2.seqfeature_id, t2.gc_percent, t2.gc_1, t2.gc_2, t2.gc_3 from custom_tables_locus2seqfeature_id t1 ' \
           ' inner join custom_tables_gc_content t2 ' \
           ' on t1.seqfeature_id=t2.seqfeature_id where t1.locus_tag="%s"' % (locus_tag)

    data = server.adaptor.execute_and_fetchall(sql1,)[0]
    taxon_id = data[0]
    seqfeature_id = data[1]
    gc_locus = data[2]
    gc_1_locus = data[3]
    gc_2_locus = data[4]
    gc_3_locus = data[5]


    sql2 = 'select gc_percent, gc_1, gc_2, gc_3 from custom_tables_gc_content where taxon_id=%s' % (taxon_id)

    data = server.adaptor.execute_and_fetchall(sql2,)
    gc_all = []
    gc_1 = []
    gc_2 = []
    gc_3 = []
    for row in data:
        gc_all.append(row[0])
        gc_1.append(row[1])
        gc_2.append(row[2])
        gc_3.append(row[3])

    gc_all_mean = round(numpy.mean(gc_all),2)
    gc_all_median = numpy.median(gc_all)

    gc_1_mean = round(numpy.mean(gc_1),2)
    gc_1_median = numpy.median(gc_1)

    gc_2_mean = round(numpy.mean(gc_2),2)
    gc_2_median = numpy.median(gc_2)

    gc_3_mean = round(numpy.mean(gc_3),2)
    gc_3_median = numpy.median(gc_3)

    path = settings.BASE_DIR + '/assets/temp/gc.svg'
    asset_path = '/temp/gc.svg'

    pairwiseid_plots.density_plot([gc_all],
                                  ['GC content'],
                                  header="",
                                  xlab="GC (%)",
                                  ylab="density",
                                  output_path=path,
                                  abline_list=[gc_locus],
                                  show_median=False,
                                  max_value=1)
    path2 = settings.BASE_DIR + '/assets/temp/gc2.svg'
    asset_path2 = '/temp/gc2.svg'
    pairwiseid_plots.density_plot([gc_1, gc_2, gc_3],
                                  ['gc 1', 'gc 2', 'gc 3'],
                                  header="",
                                  xlab="GC (%)",
                                  ylab="density",
                                  output_path=path2,
                                  abline_list=[gc_1_locus,gc_2_locus,gc_3_locus],
                                  show_median=False,
                                  max_value=1)

    return render(request, 'chlamdb/gc_locus.html', my_locals(locals()))


def get_all_prot_infos(db, seqids, orthogroups):
    hsh_gene_locs = db.get_gene_loc(seqids, as_hash=True)
    hsh_prot_infos = db.get_proteins_info(seqids)
    hsh_organisms = db.get_organism(seqids, as_hash=True)
    group_count = []
    all_locus_data = []

    for index, seqid in enumerate(seqids):
        # NOTE: all seqids are attributed an orthogroup, the case where
        # seqid is not in orthogroups should therefore not arise.
        og = orthogroups.loc[seqid].orthogroup
        if og not in group_count:
            group_count.append(og)

        strand, start, end = hsh_gene_locs[seqid]
        organism = hsh_organisms[seqid]
        locus, prot_id, gene, product = hsh_prot_infos[seqid]
        if gene==None:
            gene = ""
        data = (index, og, locus, prot_id, start, end, strand, gene, product, organism)
        all_locus_data.append(data)
    return all_locus_data, group_count


def format_orthogroup(og, to_url=False, from_str=False):
    base_str = og
    if not from_str:
        base_str = f"group_{og}"
    if to_url:
        return f"<a href=\"/orthogroup/{base_str}\">{base_str}</a>"
    return base_str


def format_locus(locus, to_url=False):
    if to_url:
        return f"<a href=\"/locusx/{locus}\">{locus}</a>"
    return locus


class FamCogColorFunc:
    def __init__(self, og, red_color):
        self.og = og
        self.red_color = red_color

    def get_color(self, taxid):
        if (self.og, taxid) in self.red_color:
            return "#FA5858"
        else:
            return EteTree.GREEN


# TODO : add error handling
def fam_cog(request, cog_id):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)
    cog_id = int(cog_id[3:])

    if request.method != "GET":
        return render(request, 'chlamdb/fam.html', my_locals(locals()))

    df_seqid_to_cog = db.get_cog_hits([cog_id], indexing="seqid", search_on="cog", keep_taxid=True)
    if len(df_seqid_to_cog)==0:
        return render(request, 'chlamdb/fam.html', {"msg": f"No entry for {format_cog(cog_id)}"})

    seqids      = df_seqid_to_cog.index.tolist()
    
    # the (group, taxid) in this dataframe are those that should be colored in red
    # in the profile (correspondance between a cog entry and an orthogroup)
    orthogroups = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    cog_info    = db.get_cog_summaries([cog_id], only_cog_desc=True, as_df=True)
    all_locus_data, group_count = get_all_prot_infos(db, seqids, orthogroups)
    red_color = set(tuple(entry) for entry in orthogroups.to_numpy())
    ref_names = db.get_genomes_description().description.to_dict()
    ref_tree  = db.get_reference_phylogeny()
    cog_func = db.get_cog_code_description()

    df_og_count  = db.get_og_count([int(i) for i in group_count], search_on="orthogroup").T
    df_cog_count = df_seqid_to_cog.groupby(["taxid"]).count()

    tree = Tree(ref_tree)
    R = tree.get_midpoint_outgroup()
    if not R is None:
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)
    e_tree.rename_leaves(ref_names)

    face_params = {"color": EteTree.RED}
    e_tree.add_column(SimpleColorColumn.fromSeries(df_cog_count.cog,
        header=format_cog(cog_id), face_params=face_params))

    for og in df_og_count:
        og_serie = df_og_count[og]
        color_chooser = FamCogColorFunc(og, red_color)
        col_column = SimpleColorColumn(og_serie.to_dict(), header=format_orthogroup(og),
                col_func=color_chooser.get_color)
        e_tree.add_column(col_column)

    asset_path = f"/temp/fam_tree_{cog_id}.svg"
    path = settings.BASE_DIR+"/assets/"+asset_path
    e_tree.render(path, dpi=500)

    group_count = [format_orthogroup(og, to_url=True) for og in group_count]
    func, cog_description = cog_info.loc[cog_id]
    info_func = "<br>".join((cog_func[code] for code in func))
    type = "cog"

    info = [info_func, cog_description]
    menu = True
    envoi = True
    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def format_pathway(pat_id):
    return f"map{pat_id:05d}"

def format_module(mod_id):
    return f"M{mod_id:05d}"

def fam_ko(request, ko_str):
    ko_id = int(ko_str[len("K"):])
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb_path, settings.BIODB_CONF)

    df_ko_hits  = db.get_ko_hits([ko_id], search_on="ko", indexing="seqid", keep_taxid=True)
    seqids      = df_ko_hits.index.tolist()
    seqid_to_og = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    red_color   = set(tuple(entry) for entry in seqid_to_og.to_numpy())

    pathways = db.get_ko_pathways([ko_id])
    modules  = db.get_ko_modules([ko_id])
    modules_id = [mod_id for key, values in modules.items() for mod_id, desc in values]
    modules_data = db.get_modules_info(modules_id)
    ko_desc      = db.get_ko_desc([ko_id])[ko_id]
    all_locus_data, group_count = get_all_prot_infos(db, seqids, seqid_to_og)

    pathway_data = format_ko_path(pathways, ko_id, as_list=True)
    module_data = [(format_ko_module(mod_id), cat, mod_desc)
            for mod_id, mod_desc, mod_def, path, cat in modules_data]

    ref_tree     = db.get_reference_phylogeny()
    leaf_to_name = db.get_genomes_description().description.to_dict()
    df_og_count  = db.get_og_count([int(og) for og in group_count], search_on="orthogroup").T
    df_ko_count  = df_ko_hits.groupby(["taxid"]).count()

    tree = Tree(ref_tree)
    R = tree.get_midpoint_outgroup()
    if not R is None:
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)
    e_tree.rename_leaves(leaf_to_name)

    face_params = {"color": EteTree.RED}
    e_tree.add_column(SimpleColorColumn.fromSeries(df_ko_count.ko,
        header=format_ko(ko_id), face_params=face_params))

    for og in df_og_count:
        og_serie = df_og_count[og]
        color_chooser = FamCogColorFunc(og, red_color)
        col_column = SimpleColorColumn(og_serie.to_dict(), header=format_orthogroup(og),
                col_func=color_chooser.get_color)
        e_tree.add_column(col_column)

    fam  = format_ko(ko_id)
    group_count = [format_orthogroup(og, to_url=True) for og in group_count]
    asset_path = f"/temp/fam_tree_{fam}.svg"
    path = settings.BASE_DIR + f"/assets/{asset_path}"

    e_tree.render(path, dpi=500)
    type = "ko"
    menu = True
    envoi = True
    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def fam(request, fam, type):
    biodb = settings.BIODB
    if request.method == 'GET': 

        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)
        # should have be redirected in fam_cog instead
        assert(type != "cog")

        print ('-- family request: biodb %s -- type %s -- name %s' % (biodb, type, fam))
        if type =='pfam':
            sql1 =   'select seqfeature_id from interpro_signature t1 ' \
                     ' inner join interpro_interpro t2 on t1.signature_id=t2.signature_id ' \
                     ' where t1.signature_accession="%s" group by seqfeature_id;' % (fam)
            sql2 = 'select signature_description from interpro where signature_accession="%s" limit 1' % (fam)
            try:
                info = server.adaptor.execute_and_fetchall(sql2, )[0]
            except:
                valid_id = False
        elif type == 'cog':
            sql1 = 'select seqfeature_id from COG_seqfeature_id2best_COG_hit t1 ' \
                   ' inner join COG_cog_names_2014 t2 on t1.hit_cog_id=t2.cog_id ' \
                   ' where COG_name="%s"' % (fam)
                   
            sql2 = 'select t2.description from COG_seqfeature_id2best_COG_hit t1 ' \
                   ' inner join COG_cog_names_2014 t2 on t1.hit_cog_id=t2.cog_id where t2.COG_name="%s";' % (fam)
            try:
                info = server.adaptor.execute_and_fetchall(sql2, )[0]
            except:
                valid_id = False
        elif type == 'interpro':
            sql1 =   'select seqfeature_id from interpro_entry t1 inner join interpro_signature t2 on t1.interpro_id=t2.interpro_id ' \
                     ' inner join interpro_interpro t3 on t2.signature_id=t3.signature_id ' \
                     ' where name="%s" group by seqfeature_id;' % (fam)
            sql2 = 'select interpro_description from interpro where interpro_accession="%s" limit 1' % (fam)
            
            try:
                info = server.adaptor.execute_and_fetchall(sql2, )[0]
            except IndexError:
                valid_id = False
        elif type == 'EC':
            sql1 = f'select distinct seqfeature_id from enzyme_seqfeature_id2ec t1' \
                   f' inner join enzyme_enzymes t2 on t1.ec_id=t2.enzyme_id where ec="{fam}";' 

            sql2 = 'select line,value from (select * from enzyme_enzymes where ec="%s") t1 ' \
                   ' inner join enzyme_enzymes_dat as t2 on t1.enzyme_id=t2.enzyme_dat_id;' % (fam)
            path = fam.split('.')

            external_link = 'https://www.qmul.ac.uk/sbcs/iubmb/enzyme/EC%s/%s/%s/%s.html' % (path[0], path[1], path[2], path[3])

            sql_pathways = 'select distinct pathway_name,pathway_category,description ' \
                           ' from (select * from enzyme_enzymes where ec = "%s") t1 ' \
                           ' inner join enzyme_kegg2ec as t2 on t2.ec_id=t1.enzyme_id ' \
                           ' inner join enzyme_kegg_pathway as t3 on t2.pathway_id=t3.pathway_id' \
                           ' where pathway_category !="1.0 Global and overview maps";' % (fam)

            pathway_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pathways, )]
            
            try:
                info =  server.adaptor.execute_and_fetchall(sql2, )
            except:
                valid_id = False

            sql_ko = f'select distinct t3.ko_accession,t3.definition from enzyme_ko2ec t1' \
                     f' inner join enzyme_enzymes t2 on t1.enzyme_id=t2.enzyme_id ' \
                     f' inner join enzyme_ko_annotation t3 on t1.ko_id=t3.ko_accession where t2.ec="{fam}";'
            associated_ko = [list(i) for i in server.adaptor.execute_and_fetchall(sql_ko, )]
            sql_ko_freq = f'select ko_accession,count(*) as n from (select distinct t3.ko_accession,seqfeature_id from enzyme_ko2ec t1 ' \
                          f' inner join enzyme_enzymes t2 on t1.enzyme_id=t2.enzyme_id' \
                          f' inner join enzyme_ko_annotation t3 on t1.ko_id=t3.ko_accession ' \
                          f' inner join enzyme_seqfeature_id2ko t4 on t3.ko_id=t4.ko_id ' \
                          f'where t2.ec="{fam}") A group by ko_accession;'
            ko2freq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_ko_freq,))

        elif type == 'ko':
            sql1 = 'select distinct seqfeature_id from enzyme_seqfeature_id2ko t1 ' \
                   ' inner join enzyme_ko_annotation t2 on t1.ko_id=t2.ko_id where t2.ko_accession="%s";' % (fam)

            sql2 = 'select * from enzyme_ko_annotation where ko_accession="%s"' % (fam)
            try:
                ko_data = server.adaptor.execute_and_fetchall(sql2,)[0]
                external_link = 'http://www.genome.jp/dbget-bin/www_bget?%s' % (fam)

                sql_modules = 'select pathways, modules from enzyme_ko_annotation where ko_accession="%s";' % (fam)
                data = server.adaptor.execute_and_fetchall(sql_modules,)[0]

                if data[0] != '-':
                    import re
                    pathway_list = [re.sub('ko', 'map',i) for i in data[0].split(',')]
                    pathway_list = '("' + '","'.join(pathway_list) + '")'
                    sql = 'select distinct pathway_name,pathway_category,description from enzyme_kegg_pathway where pathway_name in %s' % pathway_list
                    pathway_data = server.adaptor.execute_and_fetchall(sql,)
                if data[1] != '-':
                    module_list = '("' + '","'.join(data[1].split(',')) + '")'
                    sql = 'select distinct module_name,module_sub_sub_cat,description from enzyme_kegg_module where module_name in %s' % module_list
                    module_data = server.adaptor.execute_and_fetchall(sql,)
            except:
                valid_id = False
        else:
            valid_id = False
            return render(request, 'chlamdb/fam.html', my_locals(locals()))
        try:
            seqfeature_id_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql1, )]
            seqfeature_list_form = '"' + '","'.join(seqfeature_id_list) + '"'
            
        except IndexError:
            valid_id = False
            return render(request, 'chlamdb/fam.html', my_locals(locals()))
        else:
            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql2 = 'select %s from orthology_detail where seqfeature_id in (%s)' % (columns, seqfeature_list_form)

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
        from chlamdb.phylo_tree_display import ete_motifs
        if type =='pfam':
            taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb,
                                                                              [fam],
                                                                              'Pfam')

            sql3= 'select distinct taxon_id,orthogroup,signature_accession from interpro ' \
                  ' where analysis="Pfam" and orthogroup in (%s);' % ('"'+'","'.join(set(orthogroup_list))+'"')

        elif type == 'cog':
            taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb,
                                                                              [fam],
                                                                              'COG')

            sql3='select distinct taxon_id,orthogroup,COG_id from (select taxon_id,locus_tag,orthogroup ' \
                 ' from orthology_detail where orthogroup in (%s)) A ' \
                 ' inner join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag;' % ('"'+'","'.join(set(orthogroup_list))+'"')

        elif type == 'interpro':
            taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb,
                                                                              [fam],
                                                                              'interpro')

            sql3 = 'select distinct taxon_id,orthogroup,interpro_accession from ' \
                   ' interpro where orthogroup in (%s);' % ('"'+'","'.join(set(orthogroup_list))+'"')

        elif type == 'EC':
            taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb,
                                                                              [fam],
                                                                              'EC')

            group_filter = '"'+'","'.join(set(orthogroup_list))+'"'   

            sql3 = f'select distinct taxon_id,orthogroup_name,ec from enzyme_seqfeature_id2ec t1 ' \
                   f' inner join enzyme_enzymes t2 on t1.ec_id=t2.enzyme_id ' \
                   f' inner join annotation_seqfeature_id2locus t3 on t1.seqfeature_id=t3.seqfeature_id ' \
                   f' inner join orthology_seqfeature_id2orthogroup t4 on t1.seqfeature_id=t4.seqfeature_id ' \
                   f' inner join orthology_orthogroup as t5 on t4.orthogroup_id=t5.orthogroup_id ' \
                   f' where orthogroup_name in ({group_filter});'

        elif type == 'ko':
            taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb,
                                                                              [fam],
                                                                              'ko')

            sql3 = 'select t1.taxon_id, t1.orthogroup, t3.ko_accession ' \
                   ' from orthology_detail t1 ' \
                   ' inner join enzyme_seqfeature_id2ko t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                   ' inner join enzyme_ko_annotation t3 on t2.ko_id=t3.ko_id ' \
                   ' where t1.orthogroup in (%s);' % ('"'+'","'.join(set(orthogroup_list))+'"')

        else:
            taxon2orthogroup2count_reference = ete_motifs.get_taxon2name2count(biodb,
                                                                              [fam],
                                                                              'ko')
            sql3 = 'select distinct A.taxon_id,A.orthogroup,B.ko_id from (' \
                   ' select locus_tag,orthogroup,taxon_id from orthology_detail ' \
                   ' where orthogroup in (%s)) A inner join enzyme_locus2ko as B ' \
                   ' on A.locus_tag=B.locus_tag;' % ('"'+'","'.join(set(orthogroup_list))+'"')

        print(sql3)
        data = server.adaptor.execute_and_fetchall(sql3,)


        taxon2orthogroup2ec = {}
        # dico of the form:
        # taxid: {group: {domain}}
        # 58: {'group_414': ['PF06723']}, 216: {'group_414': ['PF06723']}, 269: {'group_414': ['PF06723']},
        for one_row in data:
            taxon = str(one_row[0])
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

        if len(taxon2orthogroup2ec) == 0:
            no_match = True
        else:
            taxon2orthogroup2count = ete_motifs.get_taxon2orthogroup2count(biodb, group_count)
            merged_dico = taxon2orthogroup2count
            for i in taxon2orthogroup2count_reference:
                merged_dico[i] = taxon2orthogroup2count_reference[i]

            labels = [fam] + group_count

            print("plotting")

            #print("merged_dico", merged_dico)
            #print("taxon2orthogroup2ec", taxon2orthogroup2ec)
            print("taxon2orthogroup2count_reference", taxon2orthogroup2count_reference)
            tree, style = ete_motifs.multiple_profiles_heatmap(biodb,
                                                               labels,
                                                               merged_dico,
                                                               taxon2group2value=taxon2orthogroup2ec,
                                                               highlight_first_column=True)


            if len(labels) > 30:
                big = True
                path = settings.BASE_DIR + '/assets/temp/fam_tree_%s.png' % fam
                asset_path = '/temp/fam_tree_%s.png' % fam
                tree.render(path, dpi=300, tree_style=style)
            else:
                big = False
                path = settings.BASE_DIR + '/assets/temp/fam_tree_%s.svg' % fam
                asset_path = '/temp/fam_tree_%s.svg' % fam

                tree.render(path, dpi=300, tree_style=style)

    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def fam_interpro(request, fam, type):

    # flag for template rendering
    fam_interpro = True

    if type == "Gene3D":
        Gene3D_acc = fam.split(":")[1]

    '''
    Same as fam but don't rely on comparative matrix
    retrieve pattern of presence/absence from interproscan result table
    '''
    biodb = settings.BIODB
    if request.method == 'GET': 
        from chlamdb.phylo_tree_display import ete_motifs

        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        print ('-- family request: biodb %s -- type %s -- name %s' % (biodb, type, fam))

        sql1 = f'select distinct t1.seqfeature_id from interpro.interpro_{biodb} t1 ' \
               f' inner join interpro.signature t2 on t1.signature_id=t2.signature_id ' \
               f' inner join interpro.analysis t3 on t2.analysis_id=t3.analysis_id ' \
               f' where analysis_name="{type}" and signature_accession="{fam}";'

        sql2 =   f'select signature_accession,signature_description from interpro.interpro_{biodb} t1 ' \
                 f' inner join interpro.signature t2 on t1.signature_id=t2.signature_id ' \
                 f' inner join interpro.analysis t3 on t2.analysis_id=t3.analysis_id ' \
                 f' where analysis_name="{type}" and signature_accession="{fam}" limit 1'

        try:
            info = server.adaptor.execute_and_fetchall(sql2, )[0]
        except IndexError:
            #    valid_id = False
            no_match = True
            menu = True
            return render(request, 'chlamdb/fam.html', my_locals(locals()))

        try:
            seqfeature_id_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql1, )]

            seqfeature_list_form = '"' + '","'.join(seqfeature_id_list) + '"'

        except IndexError:
            valid_id = False
            return render(request, 'chlamdb/fam.html', my_locals(locals()))
        else:

            # retrieve locus list and their annotations
            columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                      'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'
            sql2 = 'select %s from orthology_detail_%s where seqfeature_id in (%s)' % (columns, biodb, seqfeature_list_form)

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
        
        

        signature2taxon_id2count = ete_motifs.get_interpro2taxon_id2count(biodb,
                                                                          [fam],
                                                                          type)


        sql3 = 'select distinct taxon_id,orthogroup,signature_accession from ' \
                ' interpro where orthogroup in (%s) and signature_accession="%s";' % ('"'+'","'.join(set(orthogroup_list))+'"',
                                                                                       fam)


        data = server.adaptor.execute_and_fetchall(sql3,)

        # dico of the form:
        # taxid: {group: {domain}}
        # one group can be associated to multiple domains
        # 58: {'group_414': ['PF06723']}, 216: {'group_414': ['PF06723']}, 269: {'group_414': ['PF06723']},

        taxon2orthogroup2target_domain = {}
        for one_row in data:
            taxon = int(one_row[0])
            group = one_row[1]
            signature = one_row[2]
            if taxon not in taxon2orthogroup2target_domain:
                taxon2orthogroup2target_domain[taxon] = {}
                taxon2orthogroup2target_domain[taxon][group] = [signature]
            else:
                if group not in taxon2orthogroup2target_domain[taxon]:
                    taxon2orthogroup2target_domain[taxon][group] = [signature]
                else:
                    if signature not in taxon2orthogroup2target_domain[taxon][group]:
                        taxon2orthogroup2target_domain[taxon][group].append(signature)

        if len(taxon2orthogroup2target_domain) == 0:
            no_match = True
        else:
            taxon2orthogroup2count = ete_motifs.get_taxon2orthogroup2count(biodb, group_count)
            merged_dico = taxon2orthogroup2count
            for i in signature2taxon_id2count:
                merged_dico[i] = signature2taxon_id2count[i]

            labels = [fam] + group_count


            print("plotting")

            tree, style = ete_motifs.multiple_profiles_heatmap(biodb,
                                                               labels,
                                                               merged_dico, # group2taxon2count
                                                               taxon2group2value=taxon2orthogroup2target_domain,
                                                               highlight_first_column=True)


            if len(labels) > 30:
                big = True
                path = settings.BASE_DIR + '/assets/temp/fam_tree_%s.png' % fam
                asset_path = '/temp/fam_tree_%s.png' % fam
                tree.render(path, dpi=300, tree_style=style)
            else:
                big = False
                path = settings.BASE_DIR + '/assets/temp/fam_tree_%s.svg' % fam
                asset_path = '/temp/fam_tree_%s.svg' % fam

                tree.render(path, dpi=300, tree_style=style)

    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def COG_phylo_heatmap(request, frequency):
    biodb = settings.BIODB_DB_PATH

    if request.method != "GET":
        return render(request, 'chlamdb/COG_phylo_heatmap.html', my_locals(locals()))
    freq = frequency!="False"

    from ete3 import Tree
    from chlamdb.phylo_tree_display import ete_motifs
    from chlamdb.plots import cog_heatmap

    db = db_utils.DB.load_db_from_name(biodb)
    tree = db.get_reference_phylogeny()
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()
    tree, style = cog_heatmap.plot_cog_heatmap(db, t1, frequency=freq)

    freq = frequency
    path = settings.BASE_DIR + f"/assets/temp/COG_tree_{freq}.svg"
    asset_path = f"/temp/COG_tree_{freq}.svg"
    tree.render(path, dpi=600, tree_style=style)
    envoi = True
    return render(request, 'chlamdb/COG_phylo_heatmap.html', my_locals(locals()))


def plot_hist(data, xlab, ylab, title, abline, div_id):

    import plotly.graph_objects as go
    from collections import Counter
    import plotly.figure_factory as ff

    plot_data = [go.Histogram(x=data)]

    layout = go.Layout(
        title=title,
        yaxis=go.layout.YAxis(
            title=go.layout.yaxis.Title(
                text=ylab,
                font=dict(
                    family='Courier New, monospace',
                    size=15,
                    color='#7f7f7f'
                )
            )
        ),
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text=xlab,
                font=dict(
                    family='Courier New, monospace',
                    size=15,
                    color='#7f7f7f'
                )
            )
        )
    )

    fig = go.Figure()

    fig = go.Figure(data=plot_data, 
                    layout=layout)

    if abline:
        # Add abline
        fig.add_shape(
                # Line Vertical
                go.layout.Shape(
                    type="line",
                    x0=abline,
                    y0=0,
                    x1=abline,
                    y1=200,
                    line=dict(
                        color="RoyalBlue",
                        width=3
                    )
        ))

    fig.layout.margin.update({"l": 80,
                              "r": 20,
                              "b": 40,
                              "t": 80,
                              "pad": 10,
                                })

    html_plot = manipulate_biosqldb.make_div(fig, div_id=div_id)
    return html_plot


def effector_predictions(request, genome):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = f'select taxon_id from biosqldb.bioentry where accession="{genome}"'
    taxid_list = [str(server.adaptor.execute_and_fetchall(sql,)[0][0])]

    print(taxid_list)

    sql_locus_tag_effectiveT3 = 'select distinct t2.locus_tag from effectors.predicted_effectiveT3_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    locus_tag_list_effective_T3 = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql_locus_tag_effectiveT3, )]

    sql_locus_tag_T3MM = 'select distinct t2.locus_tag from effectors.predicted_T3MM_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s) and probability > 0.5;' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    locus_tag_list_effective_T3MM = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql_locus_tag_T3MM, )]

    sql_locus_tag_BPBAac = 'select distinct t2.locus_tag from effectors.predicted_BPBAac_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    locus_tag_list_effective_BPBAac = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql_locus_tag_BPBAac, )]


    sql_locus_tag_DeepT3 = 'select distinct t2.locus_tag from effectors.predicted_DeepT3_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    locus_tag_list_effective_DeepT3 = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql_locus_tag_DeepT3, )]


    serie_effective_T3 = '{' + f'name: "effectiveT3", data: {locus_tag_list_effective_T3}'  + '}' 
    serie_effective_T3MM = '{' + f'name: "T3MM", data: {locus_tag_list_effective_T3MM}'  + '}' 
    serie_effective_BPBAac = '{' + f'name: "BPBAac", data: {locus_tag_list_effective_BPBAac}'  + '}' 
    serie_effective_DeepT3 = '{' + f'name: "DeepT3", data: {locus_tag_list_effective_DeepT3}' + '}' 

    #series = [serie_pfam, serie_interpro, serie_blastnr, serie_effective_ELD] # serie_effective_T3, serie_effective_T3MM

    series = [serie_effective_T3, serie_effective_T3MM, serie_effective_BPBAac, serie_effective_DeepT3]
    series_string = '[%s]' % ','.join(series)

    envoi_venn = True



    sql_scores_effectiveT3 = 'select score from effectors.predicted_effectiveT3_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    score_list_effective_T3 = [float(i[0]) for i in server.adaptor.execute_and_fetchall(sql_scores_effectiveT3, )]

    
    sql_scores_T3MM = 'select probability from effectors.predicted_T3MM_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    score_list_T3MM = [float(i[0]) for i in server.adaptor.execute_and_fetchall(sql_scores_T3MM, )]

    
    sql_scores_BPBAac = 'select SVM_value from effectors.predicted_BPBAac_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    score_list_BPBAac = [float(i[0]) for i in server.adaptor.execute_and_fetchall(sql_scores_BPBAac, )]

    
    sql_scores_DeepT3 = 'select score from effectors.predicted_DeepT3_%s t1 ' \
                                           'inner join annotation.seqfeature_id2locus_%s t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (biodb,
                                                                           biodb,
                                                                           ','.join(taxid_list))

    # score_list_DeepT3 = [float(i[0]) for i in server.adaptor.execute_and_fetchall(sql_scores_DeepT3, )]

    
    title = 'Score distribution: effective_T3 (%s hits)' % len(score_list_effective_T3)
    xlab = 'Score'
    ylab = 'Density'
    plot_effective_T3 = plot_hist(score_list_effective_T3, xlab, ylab, title, False, "plot1")
    title = 'Score distribution: T3MM (%s hits)' % len(score_list_T3MM)
    plot_T3MM = plot_hist(score_list_T3MM, xlab, ylab, title, False, "plot2")
    title = 'Score distribution: BPBAac (%s hits)' % len(score_list_BPBAac)
    plot_BPBAac = plot_hist(score_list_BPBAac, xlab, ylab, title, False, "plot3")    
    # title = 'Score distribution: DeepT3 (%s hits)' % len(score_list_DeepT3)
    # plot_DeepT3 = plot_hist(score_list_DeepT3, xlab, ylab, title, False, "plot4")      


    

    return render(request, 'chlamdb/venn_effectors.html', my_locals(locals()))


def venn_candidate_effectors(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    taxid_list = [1279839,1279496,1035343,314,886707,804807,48,283,55,1279822,46,1279815,49,87925,52,1137444,67,1172028,1069693,1172027,307,59,60,313,1069694,62,1143376,293,1279767,1279497,1279774,64,66]
    taxid_list = [str(i) for i in taxid_list]

    '''
    taxid_list = ['67',
'1279774',
'1279496',
'1280030',
'1280034',
'1279969',
'48',
'46',
'55',
'87925',
'1279815',
'1280079',
'1279822',
'66',
'1280085',
'52',
'49',
'64',
'60',
'804807',
'886707',
'283',
'314',
'1069693',
'1280091',
'1137444',
'1280044',
'288',
'290',
'1172027',
'1172028',
'1035343',
'315',
'293',
'1117985',
'1280098',
'1280035',
'1280065',
'1280069',
'1280073',
'1279839',
'1279972',
'1279975',
'1279978',
'1279981',
'1279984',
'1279987',
'1279990',
'1279993',
'1279996',
'1279999',
'1280002',
'1280005',
'1280008',
'1280011',
'1280014',
'1280017',
'1280020',
'1279497']
    '''

    pfam_bacteria_freq_cutoff = '0.01'
    pfam_eukaryota_freq_cutoff = '0'
    interpro_euk_cutoff='95'

    sql_locus_tag_pfam = 'select distinct t5.locus_tag from interpro_interpro_signature2pfam_id t1  ' \
                         ' inner join interpro_interpro t2 on t1.signature_id=t2.signature_id ' \
                         ' inner join pfam.pfam2superkingdom_frequency_31 t3 on t1.pfam_id=t3.pfam_id ' \
                         ' inner join interpro_signature t4 on t1.signature_id=t4.signature_id ' \
                         ' inner join annotation_seqfeature_id2locus t5 on t2.seqfeature_id=t5.seqfeature_id' \
                         ' where bacteria_freq<=%s and eukaryota_freq>%s and taxon_id in (%s);' % (pam_bacteria_freq_cutoff,
                                                                                                   pfam_eukaryota_freq_cutoff,
                                                                                                   ','.join(taxid_list))

    locus_tag_list_pfam = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_pfam,)]


    sql_locus_tag_interpro = 'select distinct locus_tag from interpro_interpro t1 ' \
                             ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
                             ' inner join interpro_interpro_taxonomy_v_60 t3 on t2.interpro_id=t3.interpro_id ' \
                             ' inner join annotation_seqfeature_id2locus t5 on t1.seqfeature_id=t5.seqfeature_id ' \
                             ' where p_eukaryote>%s and taxon_id in (%s);' % (interpro_euk_cutoff,
                                                                              ','.join(taxid_list))

    locus_tag_list_interpro = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_interpro, )]



    sql_locus_tag_blast_refseq = 'select distinct locus_tag from blastnr_blastnr_best_non_self_phylum t1' \
                      ' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                      ' where superkingdom="Eukaryota" and t1.query_taxon_id in (%s);' % (','.join(taxid_list))

    locus_tag_list_blastnr = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_blast_refseq, )]

    sql_locus_tag_effectiveT3 = 'select distinct t2.locus_tag from effectors_predicted_effectiveT3 t1 ' \
                                           'inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (','.join(taxid_list))

    locus_tag_list_effective_T3 = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_effectiveT3, )]

    sql_locus_tag_T3MM = 'select distinct t2.locus_tag from effectors_predicted_T3MM t1 ' \
                                           'inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s) and probability >0.5;' % (','.join(taxid_list))

    locus_tag_list_effective_T3MM = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_T3MM, )]

    sql_locus_tag_BPBAac = 'select distinct t2.locus_tag from effectors_predicted_BPBAac t1 ' \
                                           'inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s);' % (','.join(taxid_list))

    locus_tag_list_effective_BPBAac = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_BPBAac, )]

    sql_locus_tag_ELD = 'select distinct t2.locus_tag from effectors_predicted_ELD t1 ' \
                                           'inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                           'where t1.taxon_id in (%s) and score >=10;' % (','.join(taxid_list))

    locus_tag_list_effective_ELD = [i[0] for i in server.adaptor.execute_and_fetchall(sql_locus_tag_ELD, )]


    interpro_or_pfam_not_in_blastnr = []
    for locus in locus_tag_list_interpro:
        if locus not in locus_tag_list_blastnr and locus not in interpro_or_pfam_not_in_blastnr:
            interpro_or_pfam_not_in_blastnr.append(locus)
    for locus in locus_tag_list_pfam:
        if locus not in locus_tag_list_blastnr and locus not in interpro_or_pfam_not_in_blastnr:
            interpro_or_pfam_not_in_blastnr.append(locus)

    sql = 'select species, count(*) as n from (' \
          ' select t1.locus_tag,t2.subject_taxon_id from annotation_seqfeature_id2locus t1 ' \
          ' left join blastnr_blastnr_best_non_self_phylum t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' where locus_tag in (%s)) A left join blastnr_blastnr_taxonomy B on  A.subject_taxon_id=taxon_id ' \
          ' group by family order by n DESC;' % ('"'+'","'.join(interpro_or_pfam_not_in_blastnr)+'"')

    data = server.adaptor.execute_and_fetchall(sql, )
    for row in data:
        print('%s\t%s' % (row[0], row[1]))
    serie_pfam = '{name: "pfam", data: %s}' % locus_tag_list_pfam
    serie_interpro = '{name: "interpro", data: %s}' % locus_tag_list_interpro
    serie_blastnr = '{name: "blastnr", data: %s}' % locus_tag_list_blastnr
    serie_effective_T3 = '{name: "effectiveT3", data: %s}' % locus_tag_list_effective_T3
    serie_effective_T3MM = '{name: "T3MM", data: %s}' % locus_tag_list_effective_T3MM
    serie_effective_BPBAac = '{name: "BPBAac", data: %s}' % locus_tag_list_effective_BPBAac
    serie_effective_ELD = '{name: "ELD", data: %s}' % locus_tag_list_effective_ELD

    #series = [serie_pfam, serie_interpro, serie_blastnr, serie_effective_ELD] # serie_effective_T3, serie_effective_T3MM

    series = [serie_effective_T3, serie_effective_T3MM, serie_effective_BPBAac]
    series = [serie_pfam, serie_interpro,serie_blastnr]
    series_string = '[%s]' % ','.join(series)
    pfam2description = {}
    envoi_venn = True
    return render(request, 'chlamdb/venn_euk_domains.html', my_locals(locals()))


def pfam_taxonomy_with_homologs(request, bacteria_freq, eukaryote_freq):


    from ete3 import TreeStyle

    biodb = settings.BIODB

    if request.method == 'GET': 
        from chlamdb.phylo_tree_display import ete_motifs
        server, db = manipulate_biosqldb.load_db(biodb)

        taxid_list = ['1279839', '1279496', '1035343', '314', '886707', '804807', '48', '283', '55', '1279822', '46', '1279815', '49', '87925', '52', '1137444', '67', '1172028', '1069693', '1172027', '307', '59', '60', '313', '1069694', '62', '1143376', '293', '1279767', '1279497', '1279774', '64', '66']
        taxid_list = ['314', '886707', '804807', '48', '283', '55', '1279822', '46', '1279815', '49', '87925', '52']

        '''
        taxid_list = ['67',
                      '1279774',
                      '1279496',
                      '1280030',
                      '1280034',
                      '1279969',
                      '48',
                      '46',
                      '55',
                      '87925',
                      '1279815',
                      '1280079',
                      '1279822',
                      '66',
                      '1280085',
                      '52',
                      '49',
                      '64',
                      '60',
                      '804807',
                      '886707',
                      '283',
                      '314',
                      '1069693',
                      '1280091',
                      '1137444',
                      '1280044',
                      '288',
                      '290',
                      '1172027',
                      '1172028',
                      '1035343',
                      '315',
                      '293',
                      '1117985',
                      '1280098',
                      '1280035',
                      '1280065',
                      '1280069',
                      '1280073',
                      '1279839',
                      '1279972',
                      '1279975',
                      '1279978',
                      '1279981',
                      '1279984',
                      '1279987',
                      '1279990',
                      '1279993',
                      '1279996',
                      '1279999',
                      '1280002',
                      '1280005',
                      '1280008',
                      '1280011',
                      '1280014',
                      '1280017',
                      '1280020',
                      '1279497']
        '''
        sql = 'select distinct signature_accession,signature_description from ' \
              ' (select t4.*,t5.* from interpro_interpro_signature2pfam_id t1 ' \
              ' inner join interpro_interpro t2 on t1.signature_id=t2.signature_id ' \
              ' inner join pfam.pfam2superkingdom_frequency_31 t3 on t1.pfam_id=t3.pfam_id ' \
              ' inner join interpro_signature t4 on t1.signature_id=t4.signature_id ' \
              ' inner join annotation_seqfeature_id2locus t5 on t2.seqfeature_id=t5.seqfeature_id ' \
              ' where bacteria_freq<=%s and eukaryota_freq>=%s and taxon_id in (%s)) A;' % (bacteria_freq,
                                                                                            eukaryote_freq,
                                                                                            ','.join(
                                                                                                taxid_list))

        data = server.adaptor.execute_and_fetchall(sql, )
        pfam2description = {}
        for row in data:
            pfam2description[row[0]] = row[1]

        '''
        # number of groups with identified signature domains
        sql = 'select name,n from (select AA.interpro_id, count(*) as n from ' \
              ' (select distinct A.interpro_id, B.orthogroup_id from ' \
              ' (select distinct seqfeature_id,t2.interpro_id from interpro_interpro t1 ' \
              ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
              ' inner join interpro_interpro_taxonomy_v_60 t3 on t2.interpro_id=t3.interpro_id where %s>=%s) A ' \
              ' inner join orthology_seqfeature_id2orthogroup B on A.seqfeature_id=B.seqfeature_id ' \
              ' inner join orthology_orthogroup C on B.orthogroup_id=C.orthogroup_id) AA ' \
              ' group by AA.interpro_id) BB inner join interpro_entry CC on BB.interpro_id=CC.interpro_id;' % (domain,
                                                                                                               percentage)
        '''


        # get list of all orthogroups with corresponding interpro entry
        sql = 'select signature_accession, t7.orthogroup_name from interpro_interpro_signature2pfam_id t1 ' \
              ' inner join interpro_interpro t2 on t1.signature_id=t2.signature_id ' \
              ' inner join pfam.pfam2superkingdom_frequency_31 t3 on t1.pfam_id=t3.pfam_id ' \
              ' inner join interpro_signature t4 on t1.signature_id=t4.signature_id ' \
              ' inner join annotation_seqfeature_id2locus t5 on t2.seqfeature_id=t5.seqfeature_id' \
              ' inner join orthology_seqfeature_id2orthogroup t6 on t2.seqfeature_id=t6.seqfeature_id ' \
              ' inner join orthology_orthogroup t7 on t6.orthogroup_id=t7.orthogroup_id' \
              ' where bacteria_freq<=%s and eukaryota_freq>%s and taxon_id in (%s) ' \
              ' group by signature_description, t7.orthogroup_name;' % (bacteria_freq,
                                                                        eukaryote_freq,
                                                                        ','.join(
                                                                        taxid_list))


        orthogroup_data = server.adaptor.execute_and_fetchall(sql, )

        pfam2orthogroups = {}
        orthogroup_list = []
        for i in orthogroup_data:
            if i[0] not in pfam2orthogroups:
                pfam2orthogroups[i[0]] = [i[1]]
            else:
                pfam2orthogroups[i[0]].append(i[1])
            orthogroup_list.append(i[1])

        pfam_accession2n_groups = {}
        for pfam in pfam2orthogroups:
            pfam_accession2n_groups[pfam] = len(pfam2orthogroups[pfam])


        pfam_list = list(set(pfam2orthogroups.keys()))

        taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup", taxon_filter=taxid_list)

        taxon2pfam2count = ete_motifs.get_taxon2name2count(biodb, pfam_list, type="Pfam", taxon_filter=taxid_list)

        # get total euk domain per taxon_id
        taxon_id2total = {}
        for taxon_id in taxid_list:
            taxon_id2total[taxon_id] = 0
        for accession in taxon2pfam2count:
            for taxon_id in taxid_list:
                if int(taxon2pfam2count[accession][taxon_id]) > 0 :
                    taxon_id2total[taxon_id]+=1

        taxon2pfam2count['TOTAL'] = {}
        for taxon in taxon_id2total:
            taxon2pfam2count['TOTAL'][taxon] = taxon_id2total[taxon]

        # rename interpro accession
        label_list = []
        for pfam_accession in list(taxon2pfam2count.keys()):
            if pfam_accession == 'TOTAL':
                continue
            try:
                pfam_des = "%s: %s (%s groups)" % (
                    pfam_accession,
                    pfam2description[pfam_accession],
                    pfam_accession2n_groups[pfam_accession])
            except:
                interpro_des = "%s: %s (? groups)" % (pfam_accession,
                                                      taxon2pfam2count[pfam_accession])
            taxon2pfam2count[pfam_des] = taxon2pfam2count[pfam_accession]
            pfam2orthogroups[pfam_des] = pfam2orthogroups[pfam_accession]
            if pfam_des not in label_list:
                label_list.append(pfam_des)

        # reporter interpro_accessions based on their frequency
        accession2count =  {}
        for i in label_list:
            accession2count[i] = sum([taxon2pfam2count[i][n] for n in taxon2pfam2count[i]])
        import pandas
        freq_table =  pandas.DataFrame.from_dict(accession2count, orient='index')
        sort = freq_table.sort_values(freq_table.columns[0], ascending=False)

        label_list = ['TOTAL'] + list(sort.index)

        tree, style = ete_motifs.multiple_profiles_heatmap(biodb,
                                                           label_list,
                                                           taxon2pfam2count,
                                                           rotate=True, column_scale=True)
        style.rotation=90



        tree2, tss = ete_motifs.combined_profiles_heatmap(biodb,
                                                     label_list,
                                                     taxon2orthogroup2count,
                                                     taxon2pfam2count,
                                                     pfam2orthogroups,
                                                     rotate=True,
                                                          column_scale=True)


        module_name = "taxonomy"

        big = False
        path = settings.BASE_DIR + '/assets/temp/pfam_tree_%s.svg' % module_name
        asset_path = '/temp/pfam_tree_%s.svg' % module_name
        tree.render(path, dpi=800, w=600, tree_style=style)

        path2 = settings.BASE_DIR + '/assets/temp/pfam_tree_%s_complete.svg' % module_name
        asset_path2 = '/temp/pfam_tree_%s_complete.svg' % module_name

        tree2.render(path2, dpi=800, w=600, tree_style=tss)

        envoi = True
        menu = True
        valid_id = True

    return render(request, 'chlamdb/interpro_taxonomy_homologs.html', my_locals(locals()))


def interpro_taxonomy_with_homologs(request, domain, percentage):


    from ete3 import TreeStyle

    biodb = settings.BIODB

    if request.method == 'GET': 
        from chlamdb.phylo_tree_display import ete_motifs
        server, db = manipulate_biosqldb.load_db(biodb)

        taxid_list = ['1279839', '1279496', '1035343', '314', '886707', '804807', '48', '283', '55', '1279822', '46', '1279815', '49', '87925', '52', '1137444', '67', '1172028', '1069693', '1172027', '307', '59', '60', '313', '1069694', '62', '1143376', '293', '1279767', '1279497', '1279774', '64', '66']
        #taxid_list = ['314', '886707', '804807', '48', '283', '55', '1279822', '46', '1279815', '49', '87925', '52']
        '''
        taxid_list = ['67',
                      '1279774',
                      '1279496',
                      '1280030',
                      '1280034',
                      '1279969',
                      '48',
                      '46',
                      '55',
                      '87925',
                      '1279815',
                      '1280079',
                      '1279822',
                      '66',
                      '1280085',
                      '52',
                      '49',
                      '64',
                      '60',
                      '804807',
                      '886707',
                      '283',
                      '314',
                      '1069693',
                      '1280091',
                      '1137444',
                      '1280044',
                      '288',
                      '290',
                      '1172027',
                      '1172028',
                      '1035343',
                      '315',
                      '293',
                      '1117985',
                      '1280098',
                      '1280035',
                      '1280065',
                      '1280069',
                      '1280073',
                      '1279839',
                      '1279972',
                      '1279975',
                      '1279978',
                      '1279981',
                      '1279984',
                      '1279987',
                      '1279990',
                      '1279993',
                      '1279996',
                      '1279999',
                      '1280002',
                      '1280005',
                      '1280008',
                      '1280011',
                      '1280014',
                      '1280017',
                      '1280020',
                      '1279497']
        '''
        sql = 'select * from (select distinct interpro_accession from interpro where ' \
              ' interpro_accession!="0" and taxon_id in (%s))A inner join interpro_entry B on A.interpro_accession=B.name ' \
              ' inner join interpro_interpro_taxonomy_v_60 C on B.interpro_id=C.interpro_id where %s>=%s;' % (','.join(
                                                                                                                  taxid_list),
                                                                                                              domain,
                                                                                                              percentage)

        data = server.adaptor.execute_and_fetchall(sql, )
        interpro2description = {}
        for row in data:
            interpro2description[row[0]] = row[3]

        '''
        # number of groups with identified signature domains
        sql = 'select name,n from (select AA.interpro_id, count(*) as n from ' \
              ' (select distinct A.interpro_id, B.orthogroup_id from ' \
              ' (select distinct seqfeature_id,t2.interpro_id from interpro_interpro t1 ' \
              ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
              ' inner join interpro_interpro_taxonomy_v_60 t3 on t2.interpro_id=t3.interpro_id where %s>=%s) A ' \
              ' inner join orthology_seqfeature_id2orthogroup B on A.seqfeature_id=B.seqfeature_id ' \
              ' inner join orthology_orthogroup C on B.orthogroup_id=C.orthogroup_id) AA ' \
              ' group by AA.interpro_id) BB inner join interpro_entry CC on BB.interpro_id=CC.interpro_id;' % (domain,
                                                                                                               percentage)
        '''


        # get list of all orthogroups with corresponding interpro entry
        sql = 'select name,orthogroup_name from (select distinct A.interpro_id, C.orthogroup_name from ' \
              '(select distinct seqfeature_id,t2.interpro_id from interpro_interpro t1 ' \
              ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
              ' inner join interpro_interpro_taxonomy_v_60 t3 on t2.interpro_id=t3.interpro_id where %s>=%s) A ' \
              ' inner join orthology_seqfeature_id2orthogroup B on A.seqfeature_id=B.seqfeature_id ' \
              ' inner join annotation_seqfeature_id2locus D on A.seqfeature_id=D.seqfeature_id' \
              ' inner join orthology_orthogroup C on B.orthogroup_id=C.orthogroup_id where D.taxon_id in (%s)) BB ' \
              ' inner join interpro_entry CC on BB.interpro_id=CC.interpro_id;' % (domain,
                                                                                   percentage,
                                                                                   ','.join(taxid_list))


        orthogroup_data = server.adaptor.execute_and_fetchall(sql, )

        interpro2orthogroups = {}
        orthogroup_list = []
        for i in orthogroup_data:
            if i[0] not in interpro2orthogroups:
                interpro2orthogroups[i[0]] = [i[1]]
            else:
                interpro2orthogroups[i[0]].append(i[1])
            orthogroup_list.append(i[1])

        interpro_accession2n_groups = {}
        for interpro in interpro2orthogroups:
            interpro_accession2n_groups[interpro] = len(interpro2orthogroups[interpro])


        interpro_list = list(set(interpro2orthogroups.keys()))

        taxon2orthogroup2count = ete_motifs.get_taxon2name2count(biodb, orthogroup_list, type="orthogroup", taxon_filter=taxid_list)

        taxon2interpro2count = ete_motifs.get_taxon2name2count(biodb, interpro_list, type="interpro", taxon_filter=taxid_list)

        # get total euk domain per taxon_id
        taxon_id2total = {}
        for taxon_id in taxid_list:
            taxon_id2total[taxon_id] = 0
        for accession in taxon2interpro2count:
            for taxon_id in taxid_list:
                if int(taxon2interpro2count[accession][taxon_id]) > 0 :
                    taxon_id2total[taxon_id]+=1
        taxon2interpro2count['TOTAL'] = {}
        for taxon in taxon_id2total:
            taxon2interpro2count['TOTAL'][taxon] = taxon_id2total[taxon]

        # rename interpro accession
        label_list = []
        for interpro_accession in list(taxon2interpro2count.keys()):
            if interpro_accession == 'TOTAL':
                continue
            try:
                interpro_des = "%s: %s (%s groups)" % (
                    interpro_accession,
                    interpro2description[interpro_accession],
                    interpro_accession2n_groups[interpro_accession])
            except:
                interpro_des = "%s: %s (? groups)" % (interpro_accession,
                                                      interpro2description[interpro_accession])
            taxon2interpro2count[interpro_des] = taxon2interpro2count[interpro_accession]
            interpro2orthogroups[interpro_des] = interpro2orthogroups[interpro_accession]
            if interpro_des not in label_list:
                label_list.append(interpro_des)





        # reporter interpro_accessions based on their frequency
        accession2count =  {}
        for i in label_list:
            accession2count[i] = sum([taxon2interpro2count[i][n] for n in taxon2interpro2count[i]])
        import pandas
        freq_table =  pandas.DataFrame.from_dict(accession2count, orient='index')
        sort = freq_table.sort_values(freq_table.columns[0], ascending=False)

        label_list = ['TOTAL']+list(sort.index)

        tree, style = ete_motifs.multiple_profiles_heatmap(biodb,
                                                           label_list,
                                                           taxon2interpro2count,
                                                           rotate=True, column_scale=True)
        style.rotation=90



        tree2, tss = ete_motifs.combined_profiles_heatmap(biodb,
                                                     label_list,
                                                     taxon2orthogroup2count,
                                                     taxon2interpro2count,
                                                     interpro2orthogroups,
                                                     rotate=True,
                                                          column_scale=True)


        module_name = "taxonomy"

        big = False
        path = settings.BASE_DIR + '/assets/temp/interpro_tree_%s.svg' % module_name
        asset_path = '/temp/interpro_tree_%s.svg' % module_name
        tree.render(path, dpi=800, w=600, tree_style=style)

        path2 = settings.BASE_DIR + '/assets/temp/interpro_tree_%s_complete.svg' % module_name
        asset_path2 = '/temp/interpro_tree_%s_complete.svg' % module_name

        tree2.render(path2, dpi=800, w=600, tree_style=tss)

        ko_url = '+' + '+'.join(interpro_list)
        envoi = True
        menu = True
        valid_id = True

    return render(request, 'chlamdb/interpro_taxonomy_homologs.html', my_locals(locals()))


def KEGG_module_map(request, module_name):
    from metagenlab_libs.ete_phylo import EteTree, SimpleColorColumn, ModuleCompletenessColumn 
    from metagenlab_libs.KO_module import ModuleParser
    from ete3 import Tree

    if request.method != "GET":
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    try:
        module_id = int(module_name[len("M"):])
    except:
        # add error message: module not formated correctly
        valid_id = False
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))

    module_infos = db.get_modules_info([module_id])
    if len(module_infos) != 1:
        # add error message
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))
    else:
        mod_id, module_descr, module_def, cat, sub_cat = module_infos[0]

    parser = ModuleParser(module_def)
    expr_tree = parser.parse()

    ko_ids = db.get_module_kos(module_id)
    map_data = [(format_ko(ko_id), ko_desc) for ko_id, ko_desc in db.get_ko_desc(ko_ids).items()]
    leaf_to_name = db.get_genomes_description().description.to_dict()
    mat = db.get_ko_count_for_ko(ko_ids)
    if len(mat.index) == 0:
        # should add an error message: no gene was associated for any
        # of the KO of the current module
        envoi = True
        menu = True
        valid_id = True
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))

    mat = mat.set_index(["bioentry", "ko_id"]).unstack(level=1, fill_value=0)
    mat.columns = [col for col in mat["count"].columns.values]
    tree = Tree(db.get_reference_phylogeny())
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    for column in mat.columns:
        values = mat[column].to_dict()
        new_col = SimpleColorColumn(values)
        new_col.header = format_ko(column)
        e_tree.add_column(new_col)

    hsh_n_missing = {}
    for bioentry, _ in leaf_to_name.items():
        index = int(bioentry)
        if index not in mat.index:
            n_missing = expr_tree.get_n_missing({})
        else:
            n_missing = expr_tree.get_n_missing(mat.loc[index].to_dict())
        hsh_n_missing[index] = n_missing

    completeness = ModuleCompletenessColumn(hsh_n_missing, "Missing")
    e_tree.add_column(completeness)
    e_tree.rename_leaves(leaf_to_name)

    big = len(mat.columns) >= 40
    dpi = 800 if big else 1200
    path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.svg' % module_name
    asset_path = '/temp/KEGG_tree_%s.svg' % module_name
    e_tree.render(path, dpi=dpi)
    envoi = True
    menu = True
    valid_id = True
    return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))


def kegg_multi(request, map_name, ko_name):
    biodb = settings.BIODB


    #cache.clear()

    if request.method == 'GET': 
        from chlamdb.phylo_tree_display import ete_motifs
        from chlamdb.plots import kegg_maps

        server, db = manipulate_biosqldb.load_db(biodb)


        sql = 'select bb.ko_accession, count(*) from (select C.ko_accession, E.seqfeature_id from  ' \
              ' (select * from enzyme_kegg_pathway  where pathway_name="%s") A ' \
              ' inner join enzyme_pathway2ko as B  on A.pathway_id=B.pathway_id ' \
              ' inner join enzyme_ko_annotation as C on B.ko_id=C.ko_id ' \
              ' inner join enzyme_seqfeature_id2ko E on B.ko_id=E.ko_id ' \
              ' group by B.ko_id,seqfeature_id) bb group by ko_accession;' % (map_name)


        ko2freq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql = 'select node_id from enzyme_pathway2ortholog_associations t1 ' \
              ' inner join enzyme_kegg_pathway  t2 on t1.pathway_id=t2.pathway_id ' \
              ' where pathway_name="%s" and ko_id="%s";' % (map_name, ko_name)

        node_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

        sql2 = 'select ko_id from enzyme_pathway2ortholog_associations t1 ' \
               ' inner join enzyme_kegg_pathway  t2 on t1.pathway_id=t2.pathway_id ' \
               ' where pathway_name="%s" and node_id=%s;' % (map_name, node_id)

        ko_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql2,)]
        ko_filter = '"'+ '","'.join(ko_list)+'"'
        sql = 'select ko_accession,name,definition,modules,pathways from enzyme_ko_annotation where ko_accession in (%s)' % (ko_filter)
        ko_data = server.adaptor.execute_and_fetchall(sql,)


    return render(request, 'chlamdb/KEGG_map_multi.html', my_locals(locals()))


def KEGG_mapp_ko(request, map_name):
    biodb = settings.BIODB

    if request.method == 'GET': 

        task = KEGG_map_ko_task.delay(biodb, map_name)
        task_id = task.id

    return render(request, 'chlamdb/KEGG_map_ko.html', my_locals(locals()))



def KEGG_mapp_ko_organism(request, map_name, taxon_id):
    biodb = settings.BIODB

    if request.method == 'GET': 
 
        task = KEGG_map_ko_organism_task.delay(biodb, map_name, taxon_id)
        task_id = task.id

    return render(request, 'chlamdb/KEGG_map_ko.html', my_locals(locals()))


def KEGG_mapp(request, map_name):
    biodb = settings.BIODB
    if request.method == 'GET': 
        from chlamdb.phylo_tree_display import ete_motifs
        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select pathway_name,pathway_category,description,ec,value from ' \
              '(select pathway_name,pathway_category,description,ec_id from enzyme_kegg_pathway as t1 ' \
              'inner join enzyme_kegg2ec as t2 on t1.pathway_id=t2.pathway_id where pathway_name="%s") A ' \
              'inner join enzyme_enzymes as B on A.ec_id=B.enzyme_id inner join enzyme_enzymes_dat on enzymes_dat.enzyme_dat_id=enzyme_id ' \
              'where line="description";' % (map_name)

        map_data = server.adaptor.execute_and_fetchall(sql,)

        if len(map_data) == 0:
            return KEGG_mapp_ko(request, map_name)

        enzyme_list = [i[3] for i in map_data]

        sql = 'select id from comparative_tables_EC where id in (%s);' % ('"' + '","'.join(enzyme_list) + '"')
        enzyme_list_found_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

        # get list of all orthogroups with corresponding EC
        sql = 'select distinct ec,orthogroup from enzyme_locus2ec as t1 ' \
              ' inner join enzyme_enzymes as t2 on t1.ec_id=t2.enzyme_id where ec in (%s);' % ('"' + '","'.join(enzyme_list_found_in_db) + '"')
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

        labels = enzyme_list_found_in_db
        tree, style = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2enzyme2count)

        tree2, style2 = ete_motifs.combined_profiles_heatmap(biodb,
                                                     labels,
                                                     taxon2orthogroup2count,
                                                     taxon2enzyme2count,
                                                     ec2orthogroups)


        if len(labels) > 40:
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
        envoi = True
        menu = True
        valid_id = True


    return render(request, 'chlamdb/KEGG_map.html', my_locals(locals()))


def sunburst(request, locus):
    biodb = settings.BIODB


    #cache.clear()

    if request.method == 'GET': 
        valid_id = True

        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession from orthology_detail where locus_tag = "%s"' % (locus)
        accession = server.adaptor.execute_and_fetchall(sql,)[0][0]

        sql1 = 'select t3.superkingdom,  t3.phylum,  t3.order,  t3.family,  t3.genus,  t3.species  from ' \
                   ' blastnr_blastnr_hits_%s as t1' \
                   ' inner join blastnr_blastnr_taxonomy as t3 on ' \
                   ' t1.subject_taxid = t3.taxon_id inner join blastnr_blastnr_hsps_%s as t4 ' \
                   ' on t1.nr_hit_id=t4.nr_hit_id where t1.locus_tag="%s"' % (accession, accession, locus)

        try:
            '''
            sql1 = 'select t3.superkingdom,  t3.phylum,  t3.order,  t3.family,  t3.genus,  t3.species  from ' \
                   ' blastnr_blastnr_hits_%s as t1  inner join blastnr_blastnr_hits_taxonomy_filtered_%s ' \
                   ' as t2 on t1.nr_hit_id = t2.nr_hit_id  inner join blastnr_blastnr_taxonomy as t3 on ' \
                   ' t2.subject_taxon_id = t3.taxon_id inner join blastnr_blastnr_hsps_%s as t4 ' \
                   ' on t1.nr_hit_id=t4.nr_hit_id where t1.locus_tag="%s"' % (accession, accession, accession, locus)
            '''

            raw_data = server.adaptor.execute_and_fetchall(sql1,)

        except:
            valid_id = False
            return render(request, 'chlamdb/sunburst.html', my_locals(locals()))

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

                                #print '\'%s, %s, %s, %s\\n\' +' % (i, 1, superkingdom, 0)
                                #print "\'%s, %s, %s, %s\\n\' +" % (i, 2, phylum, 0)
                                #print "\'%s, %s, %s, %s\\n\' +" % (i, 3, order, 0)
                                #print "\'%s, %s, %s, %s\\n\' +" % (i, 4, family, 0)
                                #print "\'%s, %s, %s, %s\\n\' +" % (i, 5, genus, 0)
                                #print "\'%s, %s, %s, %s\\n\' +" % (i, 6, species, dico[superkingdom][phylum][order][family][genus][species])
        out.close()

        envoi = True
        menu = True


    return render(request, 'chlamdb/sunburst.html', my_locals(locals()))


def get_cog(request, taxon_id, category):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)

    cog_hits = db.get_cog_hits([int(taxon_id)], indexing="seqid", search_on="taxid")
    cog_ids = cog_hits.cog.unique().tolist()
    cog_summaries = db.get_cog_summaries(cog_ids, only_cog_desc=True)
    prot_infos = db.get_proteins_info(cog_hits.index.unique().to_list())
    organisms = db.get_genomes_description().description.to_dict()
    functions = db.get_cog_code_description()

    data = []
    organism = organisms[int(taxon_id)]

    for seqid, cog_hit_data in cog_hits.iterrows():
        cog_id = cog_hit_data.cog
        if cog_id not in cog_summaries:
            print("Cog id unknown: ", cog_id)
            continue
        cog_func, cog_desc = cog_summaries[cog_id]
        if category not in cog_func:
            continue
        locus_tag = prot_infos[seqid][0]
        product = prot_infos[seqid][3]
        data.append([organism, locus_tag, format_cog(cog_id), cog_desc, product])

    data_type = "cog"
    description = functions[category]
    return render(request, 'chlamdb/cog_info.html', my_locals(locals()))


def get_cog_multiple(request, category, accessions=False):
    biodb = settings.BIODB
    if accessions == 'False' or accessions == 'F':
        accessions = False

    '''
    idem as get_cog but possibility to get more complex requests:
    - one or multiple include taxons
    - one or multiple explude taxons
    return the list of match COGs with their annotations
    '''
    from chlamdb.biosqldb import biosql_own_sql_tables

    server, db = manipulate_biosqldb.load_db(biodb)
    include = [i for i in request.GET.getlist('i')]
    exclude = [i for i in request.GET.getlist('o')]
    n_missing = request.GET.getlist('m')[0]
    if exclude[0] == '':
        exclude = []
        target_taxons = include
    else:
        target_taxons = include + exclude
    freq_missing = (len(include)-float(n_missing))/len(include)

    biodb_id_sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    biodb_id = server.adaptor.execute_and_fetchall(biodb_id_sql,)[0][0]

    # get sub matrix and complete matrix
    mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                  "COG",
                                                                  "id",
                                                                  include,
                                                                  exclude,
                                                                  freq_missing,
                                                                  accessions=accessions,
                                                                  cache=cache)

    match_groups_subset = mat.index.tolist()
    filter = '"' + '","'.join(match_groups_subset) + '"'
    sql = 'select COG_name, t3.code,t3.description,A.description from (select * from COG_cog_names_2014 where COG_name in (%s)) A ' \
          ' inner join COG_cog_id2cog_category t2 on A.COG_id=t2.COG_id ' \
          ' inner join COG_code2category t3 on t2.category_id=t3.category_id where t3.code="%s";' % (filter, category)

    data = server.adaptor.execute_and_fetchall(sql,)

    data_type = 'cog'

    return render(request, 'chlamdb/cog_info_multiple.html', my_locals(locals()))


def get_orthogroup_multiple_cog(request, category):
    biodb = settings.BIODB
    '''
    idem as get_cog but possibility to get more complex requests:
    - one or multiple include taxons
    - one or multiple explude taxons
    return the list of match COGs with their annotations
    '''
    from chlamdb.biosqldb import biosql_own_sql_tables

    server, db = manipulate_biosqldb.load_db(biodb)
    match_groups_subset = [i for i in request.GET.getlist('h')]

    # get list of all orthogroup with at least one hit in the specified category
    filter = '"' + '","'.join(match_groups_subset) + '"'
    sql = 'select orthogroup from (select orthogroup,locus_tag from orthology_detail ' \
          ' where orthogroup in (%s)) A left join COG_locus_tag2gi_hit as B ' \
          ' on A.locus_tag=B.locus_tag left join COG_cog_names_2014 as C on B.COG_id=C.COG_id ' \
          'where function="%s" group by orthogroup;' % (filter,
                                                       category)

    orthogroup_subset = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    filter2 = '"' + '","'.join(orthogroup_subset) + '"'
    # get detailed COG annotation of all match groups
    annot_grp = ' select A.*,B.COG_id,C.* from (select orthogroup,locus_tag ' \
                ' from orthology_detail where orthogroup in (%s)) A left join COG_locus_tag2gi_hit ' \
                ' as B on A.locus_tag=B.locus_tag left join COG_cog_names_2014 as C on B.COG_id=C.COG_id;' % (filter2)

    sql2 = 'select * from COG_code2category;'
    code2category = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    code2category[None] = '-'

    sql3 = 'select orthogroup, count(*) from orthology_detail where orthogroup in (%s) group by orthogroup' % (filter2)

    orthogroup2size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))


    # orthogroup | locus_tag   | COG_id  | COG_id  | function | name
    # for each group, count the number of matches in each category
    data = server.adaptor.execute_and_fetchall(annot_grp,)

    orthogroup2category2count = {}
    locus2count = {}

    for i in data:
        if i[1] not in locus2count:
            locus2count[i[1]] = 1
        else:
            locus2count[i[1]] += 1

    for row in data:
        # new group
        if row[0] not in orthogroup2category2count:
            orthogroup2category2count[row[0]] = {}
            orthogroup2category2count[row[0]][row[4]] = (1/float(locus2count[row[1]]))
        else:
            # existing category
            if row[4] in orthogroup2category2count[row[0]]:
                orthogroup2category2count[row[0]][row[4]] += (1/float(locus2count[row[1]]))
            else:
                # new category
                orthogroup2category2count[row[0]][row[4]] = (1/float(locus2count[row[1]]))

    data_type = 'cog'

    return render(request, 'chlamdb/get_orthogroup_multiple_cog.html', my_locals(locals()))



def get_ko_multiple(request, type, category):
    biodb = settings.BIODB
    '''
    idem as module_cat_info but possibility to get more complex requests:
    - one or multiple include taxons
    - one or multiple explude taxons
    return the list of match ko with their annotations
    '''
    from chlamdb.biosqldb import biosql_own_sql_tables
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    category = re.sub('\+', ' ', category)

    include = [i for i in request.GET.getlist('i')]
    exclude = [i for i in request.GET.getlist('o')]
    n_missing = request.GET.getlist('m')[0]
    if exclude[0] == '':
        exclude = []
        target_taxons = include
    else:
        target_taxons = include + exclude
    freq_missing = (len(include)-float(n_missing))/len(include)


    # get sub matrix and complete matrix
    mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                  "ko",
                                                                  "id",
                                                                  include,
                                                                  exclude,
                                                                  freq_missing,
                                                                              cache=cache)

    match_groups_subset = mat.index.tolist()
    filter = '"' + '","'.join(match_groups_subset) + '"'
    if type == 'module':
        sql = 'select A.ko_accession,name,definition,pathways,modules,module_name, module_sub_cat,description ' \
              ' from (select * from enzyme_ko_annotation where ko_accession in (%s)) A inner join enzyme_module2ko as B ' \
              ' on A.ko_id=B.ko_id inner join enzyme_kegg_module as C on B.module_id=C.module_id where module_sub_sub_cat="%s";' % (filter, category)
    if type == 'pathway':
        sql = 'select A.ko_accession,name,definition,pathway_name,pathway_category,description from (select * from enzyme_ko_annotation ' \
              'where ko_accession in  (%s)) A inner join enzyme_pathway2ko as B on A.ko_id=B.ko_id  ' \
              ' inner join enzyme_kegg_pathway as C on B.pathway_id=C.pathway_id' \
              ' where description="%s";' % (filter, category)

    data = list(server.adaptor.execute_and_fetchall(sql,))
    if type == 'module':
        for i, info in enumerate(data):
            data[i] = list(data[i])
            data[i][7] = info[7].split('[')[0]




    data_type = 'ko'

    return render(request, 'chlamdb/ko_info_multiple.html', my_locals(locals()))


def cog_venn_subset(request, category):
    # Note: add error handling

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    targets = [int(i) for i in request.GET.getlist('h')]
    if len(targets)>5:
        targets = targets[0:6]

    cog_hits    = db.get_cog_hits(targets, indexing="taxid", search_on="taxid")
    genome_desc = db.get_genomes_description().description.to_dict()
    cog_description = db.get_cog_summaries(cog_hits.index.tolist(), only_cog_desc=True, as_df=True)
    selected_cogs   = cog_description[cog_description.function.str.contains(category)]
    cog_codes       = db.get_cog_code_description()

    cog2description_l = []
    for cog, data in selected_cogs.iterrows():
        name = data.description
        func = data.function
        functions = ",".join(f"\"{abbr}\"" for abbr in func)
        cog2description_l.append(f"h[\"{format_cog(cog)}\"] = [[{functions}], \"{name}\"]")
    cog2description = ";".join(cog2description_l)

    sel_cog_ids = selected_cogs.index
    cog_hits = cog_hits.reindex(sel_cog_ids)

    series_tab = []
    for target in targets:
        cogs = cog_hits[target]
        non_zero_cogs = cogs[cogs > 0]
        data = ",".join(f"\"{format_cog(cog)}\"" for cog, count in non_zero_cogs.iteritems())
        series_tab.append( f"{{name: \"{genome_desc[target]}\", data: [{data}]}}" )
    series = "[" + ",".join(series_tab) + "]"

    cog_func_dict = (f"\"{func}\": \"{descr}\"" for func, descr in cog_codes.items())
    cog_func_dict = "{"+",".join(cog_func_dict)+"}"
    display_form = False
    envoi_venn = True
    return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))


def ko_venn_subset(request, category):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    category = category.replace("+", " ")
    try:
        targets = [int(i) for i in request.GET.getlist('h')]
    except:
        # add an error message
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    if len(targets)>5:
        targets = targets[0:6]

    genomes   = db.get_genomes_description(targets).description.to_dict()
    ko_counts = db.get_ko_count_cat(taxon_ids=targets, category_name=category, index=False)
    ko_count  = ko_counts.groupby("taxon_id")["KO"].apply(list)

    # shameful copy/cape from venn_ko
    fmt_data = []
    ko_set = set()
    for taxid in targets:
        if not taxid in ko_count.index:
            continue
        kos = ko_count.loc[taxid]
        kos_str = ",".join(f"{to_s(format_ko(ko))}" for ko in kos)
        genome = genomes[taxid]
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{kos_str}] }}")
        ko_set = ko_set.union({ko for ko in kos })
    series = "[" + ",".join(fmt_data) + "]"

    ko_list = list(ko_set)
    ko_descriptions = db.get_ko_desc(ko_list)
    ko2description = []
    for ko, ko_desc in ko_descriptions.items():
        forbidden = "\""
        ko_item = f"h[{to_s(format_ko(ko))}] = {forbidden}{ko_desc}{forbidden};"
        ko2description.append(ko_item)

    display_form = False
    envoi_venn = True
    return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))


def module_cat_info(request, taxid, category):
    # Not really efficient code: it would be better
    # to get the cat_id associated with category to select
    # in the category list to avoid multiple string comparison

    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb_path, settings.BIODB_CONF)

    organisms = db.get_genomes_description().description.to_dict()
    taxid = int(taxid)
    if len(organisms) == 0 or taxid not in organisms:
        return render(request, 'chlamdb/cog_info.html', my_locals(locals()))
    else:
        organism = organisms[taxid]

    category = category.replace("+", " ")
    ko_counts = db.get_ko_count([taxid], keep_seqids=True, as_multi=False)
    ko_modules = db.get_ko_modules(ko_counts["KO"].values.tolist(), as_pandas=True, compact=True)
    ko_modules_info = db.get_modules_info(ko_modules["module_id"].unique().tolist(), as_pandas=True)
    filtered_modules = ko_modules_info[ko_modules_info.subcat == category]
    selected_kos = filtered_modules.merge(ko_modules, left_on="module_id",
            right_on="module_id", how="inner")["ko_id"].unique()
    selected_seqids = ko_counts[ko_counts.KO.isin(selected_kos)]
    seqids = selected_seqids["seqid"].unique().tolist()
    hsh_to_prot = db.get_proteins_info(seqids)
    hsh_ko_desc = db.get_ko_desc(selected_kos.tolist())

    # description, locus, KO, KO name, KO description
    data = []
    for index, row in selected_seqids[["KO", "seqid"]].iterrows():
        seqid, ko_id = row.seqid, row.KO
        if seqid not in hsh_to_prot:
            continue
        locus, prot_id, gene, product = hsh_to_prot[seqid]
        ko_desc = hsh_ko_desc[ko_id]
        piece = [organism, locus, format_ko(ko_id), ko_desc, product]
        data.append(piece)
    description = category
    data_type = 'ko'
    return render(request, 'chlamdb/cog_info.html', my_locals(locals()))


def to_s(f):
    return "\"" + str(f) + "\""

def js_bioentries_to_description(hsh):
    taxon_map = 'var taxon2description = { '
    mid = ",".join(f"{to_s(bioentry)}:{to_s(description)}" for bioentry, description in hsh.items())
    return taxon_map + mid + "};"


def module_barchart(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    venn_form_class = make_venn_from(db)
    if request.method != "POST":
        form = venn_form_class()
        return render(request, 'chlamdb/module_barplot.html', my_locals(locals()))

    form = venn_form_class(request.POST)
    if not form.is_valid():
        form = venn_form_class()
        return render(request, 'chlamdb/module_barplot.html', my_locals(locals()))

    taxids = form.get_taxids()
    taxon2description = db.get_genomes_description().description.to_dict()

    ko_counts = db.get_ko_count(taxids, keep_seqids=True, as_multi=False)
    ko_ids = ko_counts.KO.unique()
    ko_module_ids = db.get_ko_modules(ko_ids.tolist(), as_pandas=True, compact=True)
    ko_modules_info = db.get_modules_info(ko_module_ids["module_id"].unique().tolist(), as_pandas=True)

    merged = ko_counts.merge(ko_module_ids, left_on="KO", right_on="ko_id", how="inner")
    merged = merged.merge(ko_modules_info, left_on="module_id", right_on="module_id", how="inner")
    cat_count = merged[["taxid", "subcat", "seqid"]].groupby(["taxid", "subcat"]).nunique()
    subcategories_list = cat_count.index.unique(level="subcat").to_list()
    subcategories = ",".join(f"{to_s(cat)}" for cat in subcategories_list)
    labels = f"[{subcategories}]"

    taxon_map = js_bioentries_to_description(taxon2description)

    series_data = []

    # not ideal, but I'm really fed up with multi-indices. Be my guest
    # if you want to improve on this.
    cat_count_dict = cat_count["seqid"].to_dict()
    taxids = cat_count.index.unique(level="taxid")
    for taxid in taxids:
        entry_data = []
        for subcat in subcategories_list:
            if (taxid, subcat) in cat_count_dict:
                entry_data.append(cat_count_dict[(taxid, subcat)])
            else:
                entry_data.append(0)
        str_entry_data = (str(entry) for entry in entry_data)
        string = f"{{ label: {to_s(taxid)}, values : [" + ",".join(str_entry_data) + "]}"
        series_data.append(string)

    taxids = "?" + "&".join((f"h={i}" for i in taxids))
    series = "[" + ",".join(series_data) + "]"
    envoi = True
    form = venn_form_class()
    return render(request, 'chlamdb/module_barplot.html', my_locals(locals()))


def add_comment(request, locus_tag):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    comment_form_class = make_comment_from(locus_tag)


    if request.method == 'POST':
        form = comment_form_class(request.POST)
        if form.is_valid():
            envoi = True
            locus_tag = form.cleaned_data['locus']
            comment =  form.cleaned_data['comment']

            if '%' in comment:
                comment = comment.replace('%', '%%')

            sql = 'select * from manual_annotation where locus_tag="%s"' % locus_tag
            data = server.adaptor.execute_and_fetchall(sql,)
            if len(data) == 0:
                sql = 'insert into manual_annotation (locus_tag, annotation) values("%s", "%s")' % (locus_tag, comment)
            else:
                sql = 'update manual_annotation set annotation="%s" where locus_tag="%s"' % (comment, locus_tag)
            server.adaptor.execute(sql,)
            server.commit()
    else:
        form = comment_form_class()
    return render(request, 'chlamdb/comment_form.html', my_locals(locals()))


def add_locus_int(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST':
        form = LocusInt(request.POST)
        if form.is_valid():
            from datetime import datetime

            now = datetime.now()
            str_date = "%s-%s-%s" % (now.year, now.month, now.day)

            envoi = True

            category = form.cleaned_data['category']
            gene = form.cleaned_data['gene']
            locus_tag = form.cleaned_data['locus_tag']
            description = form.cleaned_data['description']
            reference = form.cleaned_data['reference']

            if '%' in description:
                comment = description.replace('%', '%%')
            if '%' in reference:
                comment = reference.replace('%', '%%')

            sql = 'select * from manual_annotation where locus_tag="%s"' % locus_tag
            data = server.adaptor.execute_and_fetchall(sql,)
            if len(data) > 0:
                raise("locus already present")

            else:
                sql = 'insert into custom_tables_annot_table (category,gene,locus_tag, description, ' \
                      ' reference, date) values("%s", "%s", "%s","%s","%s", "%s")' % (category,
                                                                                     gene,
                                                                                     locus_tag,
                                                                                     description,
                                                                                     reference,
                                                                                     str_date)

            server.adaptor.execute(sql,)
            server.commit()
    else:
        form = LocusInt()
    return render(request, 'chlamdb/add_inter_form.html', my_locals(locals()))

def ko_subset_barchart(request, type):
    biodb = settings.BIODB
    '''

    create KO sub sub category barchart of selected KO
    url parameters:
                    i: taxons to include
                    o: taxons to exclude
                    m: frequ missing (allow COG to miss in m incided genomes)

    GET THE LIST OF cogS PRESENT IN MINIMUM LEN(i)-freq-missing AND NOT PRESENT IN ALL(o)

    create dictionnary of counts for each category
    construct the barchart with javascript

    :param request:
    :param biodb:
    :return:
    '''

    from chlamdb.biosqldb import biosql_own_sql_tables
    server, db = manipulate_biosqldb.load_db(biodb)

    include = [i for i in request.GET.getlist('i')]
    exclude = [i for i in request.GET.getlist('o')]
    # if not exclude taxon
    if exclude[0] == '':
        exclude = []
    n_missing = request.GET.getlist('m')[0]

    freq_missing = (len(include)-float(n_missing))/len(include)

    # get sub matrix and complete matrix
    mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                  "ko",
                                                                  "id",
                                                                  include,
                                                                  exclude,
                                                                  freq_missing,
                                                                              cache=cache)

    match_groups_subset = mat.index.tolist()

    if type == 'module':
        sql = 'select module_sub_sub_cat,count(*) as n from (select * from enzyme_ko_annotation_v1 where ko_id in ' \
              ' (%s)) A inner join enzyme_module2ko_v1 as B on A.ko_id=B.ko_id ' \
              ' inner join enzyme_kegg_module_v1 as C on B.module_id=C.module_id group by module_sub_sub_cat;' % ('"'+'","'.join(match_groups_subset)+'"')
    if type == 'pathway':
        sql = 'select description,count(*) as n from (select * from enzyme_ko_annotation_v1 where ko_id in  ' \
              ' (%s)) A inner join enzyme_pathway2ko_v1 as B on A.ko_id=B.ko_id  inner join enzyme_kegg_pathway as C' \
              ' on B.pathway_id=C.pathway_id where pathway_category != "1.0 Global and overview maps"' \
              ' group by description;' % ('"'+'","'.join(match_groups_subset)+'"')

    data_subset = server.adaptor.execute_and_fetchall(sql,)


    # on récupère tous les KO de tous les génomes inclus pour faire une comparaison
    filter = '`' + '`>0 or `'.join(include) + '`>0'
    sql = 'select id from comparative_tables_ko where (%s)' % (filter)

    match_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    if type == 'module':
        sql = 'select module_sub_sub_cat,count(*) as n from (select * from enzyme_ko_annotation_v1 where ko_id in ' \
              ' (%s)) A inner join enzyme_module2ko_v1 as B on A.ko_id=B.ko_id ' \
              ' inner join enzyme_kegg_module_v1 as C on B.module_id=C.module_id group by module_sub_sub_cat;' % ('"'+'","'.join(match_groups)+'"')
    if type == 'pathway':
        sql = 'select description,count(*) as n from (select * from enzyme_ko_annotation_v1 where ko_id in ' \
              ' (%s)) A inner join enzyme_pathway2ko_v1 as B on A.ko_id=B.ko_id ' \
              ' inner join enzyme_kegg_pathway as C on B.pathway_id=C.pathway_id where pathway_category != "1.0 Global and overview maps"' \
              ' group by description;' % ('"'+'","'.join(match_groups)+'"')

    data_all_include = server.adaptor.execute_and_fetchall(sql,)

    # count total
    total_subset = sum([i[1] for i in data_subset])
    total_all_include = sum([i[1] for i in data_all_include])

    # calculate freq of each category, add missing categories in the 2 dictionnaries
    category2freq_subset = {}
    category2count_subset = {}
    for i in data_subset:
        category2freq_subset[i[0]] = float(i[1])/total_subset
        category2count_subset[i[0]] = i[1]
    category2freq_all = {}
    category2count_all = {}
    for i in data_all_include:
        category2freq_all[i[0]] = float(i[1])/total_all_include
        category2count_all[i[0]] = i[1]
        if i[0] not in category2freq_subset:
            category2freq_subset[i[0]] = 0
            category2count_subset[i[0]] = 0
    for cat in category2freq_subset:
       if cat not in category2freq_all:
           category2freq_all[cat] = 0
           category2count_all[cat] = 0

    labels_template = '[\n' \
                      '%s\n' \
                      ']\n'

    serie_template = '[%s\n' \
                     ']\n'

    one_serie_template = '{\n' \
                         'label: "%s",\n' \
                         'values: [%s]\n' \
                         '},\n'
    ref_serie = []
    all_serie = []
    for cat in category2freq_all.keys():

        ref_serie.append(str(round(category2freq_subset[cat],4)*100))
        all_serie.append(str(round(category2freq_all[cat],4)*100))

    serie_all = one_serie_template % ("complete genomes", ','.join(all_serie))
    serie_target = one_serie_template % ("selection", ','.join(ref_serie))

    series = serie_template % ''.join([serie_all, serie_target])
    labels = labels_template % ('"'+'","'.join([str(i) for i in category2freq_all.keys()]) + '"')

    category_count_complete = 'var category_count_complete = {'
    for i in category2count_all:
        category_count_complete+='"%s":["%s", "%s"],' % (i, category2count_all[i], category2count_subset[i])
    category_count_complete = category_count_complete[0:-1] + '};'

    taxons_in_url = "?i="+("&i=").join(include) + '&m=%s' % str(n_missing)
    taxon_out_url = "&o="+("&o=").join(exclude)
    return render(request, 'chlamdb/ko_subset_barchart.html', my_locals(locals()))


def compare_homologs(request):
    biodb = settings.BIODB

    '''
    plot protein length of closest homologs between 2 genomes
    :param request:
    :param biodb:
    :return:
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    venn_form_class = make_venn_from(biodb)

    if request.method == 'POST':
        form = venn_form_class(request.POST)

        if form.is_valid():
            target_taxons = form.cleaned_data['targets']
            target_taxons = target_taxons[0:2]
            sql = 'select locus_1,locus_2 from comparative_tables_identity_closest_homolog ' \
                  ' where (taxon_1=%s and taxon_2=%s) ' \
                  ' UNION select locus_2,locus_1 from comparative_tables_identity_closest_homolog ' \
                  ' where (taxon_1=%s and taxon_2=%s)' % (target_taxons[0],
                                                           target_taxons[1],
                                                           target_taxons[1],
                                                           target_taxons[0])
            sql = 'select locus_1,locus_2 from comparative_tables_identity_closest_homolog ' \
                  ' where (taxon_1=%s and taxon_2=%s)' % (target_taxons[0],
                                                           target_taxons[1])

            locus_list = list(server.adaptor.execute_and_fetchall(sql,))

            locus2length = {}
            locus2orthogroup = {}
            for taxon in target_taxons:
                sql1 = 'select locus_tag, orthogroup from orthology_detail where taxon_id=%s' % (taxon)
                if db_driver == 'mysql':
                    sql2 = 'select locus_tag, char_length(translation) as len from orthology_detail where taxon_id=%s' % (taxon)
                if db_driver == 'sqlite':
                    sql2 = 'select locus_tag, length(translation) as len from orthology_detail where taxon_id=%s' % (taxon)
                tmp1 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
                tmp2 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
                locus2length.update(tmp2)
                locus2orthogroup.update(tmp1)
            '''
            var dataset = [
            [5, 20], [480, 90], [250, 50], [100, 33], [330, 95],
            [410, 12], [475, 44], [25, 67], [85, 21], [220, 88],
            [600, 150]
            ];
            '''
            dataset_template = '''
            var dataset = [
            %s
            ];
            '''
            datastring = ''
            for locus_1, locus_2 in locus_list:

                datastring+='[%s,%s,"%s"], ' % (locus2length[locus_1],
                                            locus2length[locus_2],
                                            locus2orthogroup[locus_1])
            data_string = datastring[0:-2]

            dataset = dataset_template % data_string

    else:  
        form = venn_form_class()

    return render(request, 'chlamdb/prot_length_scatter.html', my_locals(locals()))


def orthogroup2cog_series(orthogroup_list, reference_taxon=None, accessions=False):

    if accessions=='False' or accessions == 'F':
        accessions=False

    server, db = manipulate_biosqldb.load_db(biodb)

    '''
    for each orthogroup, get the list of associated COGs
    ponderate the count value by the frequency of the COG within the group
    TODO: inclue locus without cog hits within the counts?
    '''


    if not accessions:
        sql = 'select A.orthogroup,C.function, count(*) from (select * from orthology_detail as t1 ' \
              ' where orthogroup in (%s) and taxon_id=%s) A left join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag ' \
              ' inner join COG_cog_names_2014 as C on B.COG_id=C.COG_id group by orthogroup,function;' % ('"' + '","'.join(orthogroup_list) + '"',
                                                                                                         reference_taxon)

    else:
        sql = 'select A.orthogroup,C.function, count(*) from (select * from orthology_detail as t1 ' \
              ' where orthogroup in (%s) and accession="%s") A left join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag ' \
              ' inner join COG_cog_names_2014 as C on B.COG_id=C.COG_id group by orthogroup,function;' % ('"' + '","'.join(orthogroup_list) + '"',
                                                                                                         reference_taxon)


    data = server.adaptor.execute_and_fetchall(sql,)

    # reference taxon allow to compare a subset of group with a whole genome. Not mandatory.
    if not reference_taxon:
        reference_taxon = data[0][-1]


    # get genome description
    if not accessions:
        sql = 'select t2.description from biodatabase as t1 inner join bioentry as t2 ' \
              ' on t1.biodatabase_id=t2.biodatabase_id where t1.name="%s" ' \
              ' and taxon_id=%s and t2.description not like "%%%%plasmid%%%%";' % (biodb,
                                                                                    reference_taxon)

    else:
        sql = 'select t2.description from biodatabase as t1 inner join bioentry as t2 ' \
              ' on t1.biodatabase_id=t2.biodatabase_id where t1.name="%s" ' \
              ' and accession="%s";' % (biodb,
                                        reference_taxon)

    genome_reference = server.adaptor.execute_and_fetchall(sql,)[0]


    # count cog categories for each orthogroup
    # count the total number of locus with a COG to ponderate the cog categorie counts
    orthogroup2counts = {}
    orthogroup2total_count = {}
    for row in data:
        if row[0] not in orthogroup2counts:
            orthogroup2counts[row[0]] = {}
            orthogroup2total_count[row[0]] = float(row[2])
            orthogroup2counts[row[0]][row[1]] = float(row[2])
        else:
            orthogroup2counts[row[0]][row[1]] = float(row[2])
            orthogroup2total_count[row[0]] += float(row[2])

    # calculate fraction of each category for each group
    # sum of all categories fractions = 1 for each group
    cog_category2count = {}
    for group in orthogroup2counts:
        for category in orthogroup2counts[group]:
            if category not in cog_category2count:
                cog_category2count[category] = orthogroup2counts[group][category]/orthogroup2total_count[group]
            else:
                cog_category2count[category] += orthogroup2counts[group][category]/orthogroup2total_count[group]
    # count total count, then fraction of the total for each category
    category2fraction = {}
    total = sum([cog_category2count[i] for i in cog_category2count])
    for category in cog_category2count:
        category2fraction[category] = (cog_category2count[category]/total)*100

    # get data for the whole reference genome
    if not reference_taxon:
        sql = 'select A.orthogroup,C.function, count(*) as n from (select * from orthology_detail as t1) A ' \
              ' left join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag ' \
              ' inner join COG_cog_names_2014 as C on B.COG_id=C.COG_id group by orthogroup,function;'

    else:
        if not accessions:
            sql = 'select A.orthogroup,C.function, count(*) as n from (select * from orthology_detail as t1 where taxon_id=%s) A ' \
                  ' left join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag ' \
                  ' inner join COG_cog_names_2014 as C on B.COG_id=C.COG_id group by orthogroup,function;' % (reference_taxon)
        else:
            sql = 'select A.orthogroup,C.function, count(*) as n from (select * from orthology_detail as t1 where accession="%s") A ' \
                  ' left join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag ' \
                  ' inner join COG_cog_names_2014 as C on B.COG_id=C.COG_id group by orthogroup,function;' % (reference_taxon)
    # same counts as previously, but with the whole genome
    # TODO seperate function do do those counts
    data_all = server.adaptor.execute_and_fetchall(sql,)
    orthogroup2counts_all = {}
    orthogroup2total_count_all = {}
    for row in data_all:
        if row[0] not in orthogroup2counts_all:
            orthogroup2counts_all[row[0]] = {}
            orthogroup2total_count_all[row[0]] = float(row[2])
            orthogroup2counts_all[row[0]][row[1]] = float(row[2])
        else:
            orthogroup2counts_all[row[0]][row[1]] = float(row[2])
            orthogroup2total_count_all[row[0]] += float(row[2])

    # count fraction of each category for each group
    # sum of all category/ies = 1 for each group
    cog_category2count_all = {}
    for group in orthogroup2counts_all:
        for category in orthogroup2counts_all[group]:
            if category not in cog_category2count_all:
                cog_category2count_all[category] = orthogroup2counts_all[group][category]#/orthogroup2total_count_all[group]
            else:
                cog_category2count_all[category] += orthogroup2counts_all[group][category]#/orthogroup2total_count_all[group]
    # count total count, then fraction of the total for each category
    category2fraction_all = {}
    total = sum([cog_category2count_all[i] for i in cog_category2count_all])

    # number and list of of groups without any cog hit
    n_missing_cog = len(orthogroup_list) - len(orthogroup2counts)
    missing_cog_list = []
    for group in orthogroup_list:
        if group not in orthogroup2counts:
            missing_cog_list.append(group)

    for category in cog_category2count_all:
        # get count in percent
        category2fraction_all[category] = (cog_category2count_all[category]/total)*100

    for category in category2fraction_all:
        if category not in category2fraction:
            category2fraction[category] = 0
            cog_category2count[category] = 0

    # prepare javascript code to make barcharts
    labels_template = '[\n' \
                      '%s\n' \
                      ']\n'

    serie_template = '[%s\n' \
                     ']\n'

    one_serie_template = '{\n' \
                         'label: "%s",\n' \
                         'values: [%s]\n' \
                         '},\n'

    # percentage of each category
    # comparison with counts of a complete reference genome
    ref_serie = []
    all_serie = []
    for cat in category2fraction_all.keys():
        ref_serie.append(str(round(category2fraction[cat],2)))
        all_serie.append(str(round(category2fraction_all[cat],2)))

    serie_all = one_serie_template % ("%s" % genome_reference, ','.join(all_serie))
    serie_target = one_serie_template % ("selection", ','.join(ref_serie))

    series = serie_template % ''.join([serie_all, serie_target])
    labels = labels_template % ('"'+'","'.join([str(i) for i in category2fraction_all.keys()]) + '"')

    # counts in absolute numbers
    ref_serie_counts = []
    all_serie_counts = []
    for cat in cog_category2count.keys():
        ref_serie_counts.append(str(cog_category2count[cat]))
        all_serie_counts.append(str(cog_category2count_all[cat]))

    serie_all_counts = one_serie_template % ("%s" % genome_reference, ','.join(all_serie_counts))
    serie_target_counts = one_serie_template % ("selection", ','.join(ref_serie_counts))

    series_counts = serie_template % ''.join([serie_target_counts])
    labels_counts = labels_template % ('"'+'","'.join([str(i) for i in cog_category2count.keys()]) + '"')

    sql = 'select * from COG_code2category'
    category_description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    category_map = 'var category_description = {'
    for i in category_description:
        category_map+='"%s":"%s",' % (i, category_description[i])
    category_map = category_map[0:-1] + '};'

    return series, labels, serie_all_counts, serie_target_counts, series_counts, labels_counts, category_description, category_map, n_missing_cog, missing_cog_list


def locus_tag2cog_series(locus_tag_list, reference_taxon=None):

    server, db = manipulate_biosqldb.load_db(biodb)



    # get one cog/locus tag
    sql = ' select A.*,B.taxon_id from (select t1.locus_tag,t2.function, count(*) as n from COG_locus_tag2gi_hit as t1  ' \
          ' inner join COG_cog_names_2014 as t2 on t1.COG_id=t2.COG_id where locus_tag in (%s) group by locus_tag,function) A ' \
          ' inner join orthology_detail as B on A.locus_tag=B.locus_tag;' % ('"' + '","'.join(locus_tag_list) + '"')

    sql = 'select locus_tag, taxon_id from COG_seqfeature_id2best_COG_hit t1 ' \
          ' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' where locus_tag in (%s) group by locus_tag;' % ('"' + '","'.join(locus_tag_list) + '"')

    data = server.adaptor.execute_and_fetchall(sql,)
    # non redundant list of locus with associated COG
    locus_with_cog_list = set([i[0] for i in data])

    if not reference_taxon:
        reference_taxon = data[0][-1]


    sql = 'select t2.description from biodatabase as t1 inner join bioentry as t2 ' \
          ' on t1.biodatabase_id=t2.biodatabase_id where t1.name="%s" ' \
          ' and taxon_id=%s and t2.description not like "%%%%plasmid%%%%";' % (biodb,
                                                                                reference_taxon)
    genome_reference = server.adaptor.execute_and_fetchall(sql,)[0]

    # counting COG categories
    cog_category2count = {}
    for row in data:
        if row[1] not in cog_category2count:
            cog_category2count[row[1]] = 1
        else:
            cog_category2count[row[1]] += 1

    # calculating COG categories as percentages of total COG
    total = sum([cog_category2count[i] for i in cog_category2count])

    cog_category2fraction = {}
    for category in cog_category2count:
        cog_category2fraction[category] = (cog_category2count[category]/float(total))*100

    sql = 'select A.locus_tag,C.function, count(*) as n from (select * from orthology_detail as t1 where taxon_id=%s) A ' \
          ' left join COG_locus_tag2gi_hit as B on A.locus_tag=B.locus_tag ' \
          ' inner join COG_cog_names_2014 as C on B.COG_id=C.COG_id group by orthogroup,function;' % (reference_taxon)
    # attention!!!! Pourquoi group by orthogroup ci-dessus???????!!!!!!!
    sql = 'select A.locus_tag,code, count(*) as n from (select * from orthology_detail as t1 where taxon_id=%s) A ' \
          ' left join COG_seqfeature_id2best_COG_hit as B on A.seqfeature_id=B.seqfeature_id ' \
          ' inner join COG_cog_names_2014 as C on B.hit_cog_id=C.COG_id ' \
          ' inner join COG_cog_id2cog_category D on C.COG_id=D.COG_id ' \
          ' inner join COG_code2category E on D.category_id=E.category_id group by locus_tag,code;' % (eference_taxon)




    data_all = server.adaptor.execute_and_fetchall(sql,)
    # counting COG categories for the whole genome
    cog_category2count_all = {}
    for row in data_all:
        if row[1] not in cog_category2count_all:
            cog_category2count_all[row[1]] = 1

        else:
            cog_category2count_all[row[1]] += 1

    # getting list of locus without COGs
    n_missing_cog = len(locus_tag_list) - len(locus_with_cog_list)
    missing_cog_list = list(set(locus_tag_list) - locus_with_cog_list)

    # claculating fraction of each category
    total = sum([cog_category2count_all[i] for i in cog_category2count_all])
    cog_category2fraction_all = {}
    for category in cog_category2count_all:
        # get count in percent
        cog_category2fraction_all[category] = (cog_category2count_all[category]/float(total))*100

    for category in cog_category2fraction_all:
        if category not in cog_category2count:
            cog_category2fraction[category] = 0
            cog_category2count[category] = 0

    labels_template = '[\n' \
                      '%s\n' \
                      ']\n'

    serie_template = '[%s\n' \
                     ']\n'

    one_serie_template = '{\n' \
                         'label: "%s",\n' \
                         'values: [%s]\n' \
                         '},\n'

    # preparing series for locus list and complete genome
    ref_serie = []
    all_serie = []
    for cat in cog_category2fraction_all.keys():
        ref_serie.append(str(round(cog_category2fraction[cat],2)))
        all_serie.append(str(round(cog_category2fraction_all[cat],2)))

    serie_all = one_serie_template % ("%s" % genome_reference, ','.join(all_serie))
    serie_target = one_serie_template % ("selection", ','.join(ref_serie))

    series = serie_template % ''.join([serie_all, serie_target])
    labels = labels_template % ('"'+'","'.join([str(i) for i in cog_category2fraction_all.keys()]) + '"')

    # counts in absolute numbers
    ref_serie_counts = []
    all_serie_counts = []
    for cat in cog_category2count.keys():
        ref_serie_counts.append(str(cog_category2count[cat]))
        all_serie_counts.append(str(cog_category2count_all[cat]))


    serie_all_counts = one_serie_template % ("%s" % genome_reference, ','.join(all_serie_counts))
    serie_target_counts = one_serie_template % ("selection", ','.join(ref_serie_counts))

    series_counts = serie_template % ''.join([serie_target_counts])
    labels_counts = labels_template % ('"'+'","'.join([str(i) for i in cog_category2count.keys()]) + '"')

    sql = 'select * from COG_code2category'
    category_description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    category_map = 'var category_description = {'
    for i in category_description:
        category_map+='"%s":"%s",' % (i, category_description[i])
    category_map = category_map[0:-1] + '};'

    return series, labels, serie_all_counts, serie_target_counts, series_counts, labels_counts, category_description, category_map, n_missing_cog, missing_cog_list


def orthogroup_list_cog_barchart(request, accessions=False):
    biodb = settings.BIODB
    if accessions == 'False' or accessions == 'F':
        accessions = False

    orthogroup_list = [i for i in request.GET.getlist('h')]

    reference = request.GET.getlist('ref')[0]

    series, \
    labels, \
    serie_all_counts, \
    serie_target_counts, \
    series_counts, \
    labels_counts, \
    category_description, \
    category_map, \
    n_missing_cog, \
    missing_cog_list = orthogroup2cog_series(orthogroup_list, reference_taxon=reference, accessions=accessions)


    no_cogs_url = "?g=" + ('&g=').join(missing_cog_list)
    orthogroups_url = '?h=' + ('&h=').join(orthogroup_list)

    return render(request, 'chlamdb/orthogroup_list_cog_barchart.html', my_locals(locals()))


def cog_barchart(request):
    # Not too happy with this code
    # should be reviewed and the old hashes replaced by pandas
    # operations

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)
    venn_form_class = make_venn_from(db)

    if request.method != 'POST':
        form = venn_form_class()
        return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))

    form = venn_form_class(request.POST)

    if not form.is_valid():
        # add error message
        return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))

    target_bioentries = form.get_taxids()

    hsh_counts = db.get_cog_counts_per_category(target_bioentries)
    taxon2description = db.get_genomes_description().description.to_dict()
    category_dico = db.get_cog_code_description()

    # create a dictionnary to convert cog category description and one letter code
    category_map = 'var description2category = {'
    map_lst = (f"\"{func_descr}\" : \"{func}\"" for func, func_descr in category_dico.items())
    category_map = category_map + ",".join(map_lst) + '};'

    taxon_map = 'var taxon2description = {'
    taxon_map_lst = (f"\"{target}\":\"{taxon2description[target]}\"" for target in target_bioentries)
    taxon_map = taxon_map + ",".join(taxon_map_lst) + '};'

    # Not too happy with this code and its level of indentation
    # Could also be made faster by avoiding string comparisons and list lookup
    taxon2category2count = {}
    all_categories = []
    for bioentry, hsh_cnt in hsh_counts.items():
        bioentry_str = str(bioentry)
        if bioentry_str not in taxon2category2count:
            taxon2category2count[bioentry_str] = {}

        for func, cnt in hsh_cnt.items():
            # a cog can have multiple functions
            for i in range(0, len(func)):
                f = func[i]
                category = category_dico[f]
                if category in taxon2category2count[bioentry_str]:
                    taxon2category2count[bioentry_str][category] += cnt
                else:
                    taxon2category2count[bioentry_str][category] = cnt
                    if category not in all_categories:
                        all_categories.append(category)

    labels_template = '[\n' \
                      '%s\n' \
                      ']\n'

    serie_template = '[%s\n' \
                     ']\n'

    one_serie_template = '{\n' \
                         'label: "%s",\n' \
                         'values: [%s]\n' \
                         '},\n'


    all_series_templates = []
    for taxon in taxon2category2count:
        one_category_list = []
        for category in all_categories:
            count = taxon2category2count[taxon].get(category, 0)
            one_category_list.append(count)
        one_category_list = [str(i) for i in one_category_list]

        all_series_templates.append(one_serie_template % (taxon, ','.join(one_category_list)))

    taxids = "?" + "&".join((f"h={i}" for i in target_bioentries))
    series = serie_template % ''.join(all_series_templates)
    labels = labels_template % ('"'+'","'.join(all_categories) + '"')
    envoi = True
    return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))


def get_locus_annotations(biodb, locus_list):


    '''
    get annotation from a serie of locus
    - genbank annot
    - ko annot
    - modules
    - pathways
    - cogs


    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    import re
    from string import digits
    server, db = manipulate_biosqldb.load_db(biodb)

    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
              'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, taxon_id'
    sql = 'select %s from orthology_detail where locus_tag in (%s)' % (columns, '"' + '","'.join(locus_list) + '"')

    all_data = server.adaptor.execute_and_fetchall(sql,)

    locus_list = [i[1] for i in all_data]

    locus2annot = []
    for i, data in enumerate(all_data):
        locus2annot.append((i,) + data)

    
    sql2 = 'select A.locus_tag, B.COG_name from (select locus_tag, COG_id from COG_locus_tag2gi_hit ' \
          ' where locus_tag in (%s)) A inner JOIN ' \
          ' COG_cog_names_2014 as B on A.COG_id=B.COG_name' % ('"' + '","'.join(locus_list) + '"')


    sql = 'select A.locus_tag, t3.code from (select locus_tag, COG_id from COG_locus_tag2gi_hit ' \
           ' where locus_tag in (%s)) A inner JOIN COG_cog_names_2014 B on A.COG_id=B.COG_name ' \
           ' inner join COG_cog_id2cog_category t2 on B.COG_id=t2.COG_id ' \
           ' inner join COG_code2category t3 on t2.category_id=t3.category_id;' % ('"' + '","'.join(locus_list) + '"')

    sql3 = 'select locus_tag,ko_id from enzyme_locus2ko where locus_tag in (%s) ' % ('"' + '","'.join(locus_list) + '"')
    sql4 = 'select pathway_name,pathway_category from enzyme_kegg_pathway'
    sql5 = 'select module_name,description from enzyme_kegg_module_v1'

    sql6 = 'select * from (select distinct locus_tag,interpro_accession,interpro_description ' \
           ' from interpro where locus_tag in (%s)) A where interpro_accession!="0"' % ('"' + '","'.join(locus_list) + '"')

    try:
        locus_tag2cog_catego = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        locus_tag2cog_name = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    # TODO: deal with missing data (hide columns in tables)
    # case when missing tables
    except:
        locus_tag2cog_catego = {}
        locus_tag2cog_name = {}
    try:
        locus_tag2ko = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))
        pathway2category = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))
        module2category = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql5,))
    except:
        locus_tag2ko = {}
        pathway2category = {}
        module2category = {}

    try:
        interpro_data = server.adaptor.execute_and_fetchall(sql6,)
    except:
        # TODO: deal with missing data (hide columns in tables)
        # case when missing tables
        interpro_data = []
    locus2interpro = {}
    for row in interpro_data:
        if row[0] not in locus2interpro:
            locus2interpro[row[0]] = [row[1:]]
        else:
            locus2interpro[row[0]].append(row[1:])

    for locus in locus_list:
        if locus not in locus_tag2cog_name:
            locus_tag2cog_name[locus] = '-'
            locus_tag2cog_catego[locus] = '-'
        if locus not in locus_tag2ko:
            locus_tag2ko[locus] = '-'
        if locus not in locus2interpro:
            locus2interpro[locus] = [('-', '-')]

    sql4 = 'select ko_id,pathways,modules from enzyme_ko_annotation_v1 where ko_id in (%s); ' % ('"' + '","'.join(locus_tag2ko.values()) + '"')

    try:
        ko_data = server.adaptor.execute_and_fetchall(sql4,)
    # TODO: deal with missing data (hide columns in tables)
    # case when missing tables
    except:
        ko_data = []

    ko2ko_pathways = {}
    ko2ko_modules = {}
    ko2pathway_categories = {}
    for one_ko in ko_data:
        if one_ko[1] != '-':
            ko2ko_pathways[one_ko[0]] = ''
            ko2pathway_categories[one_ko[0]] = ''
            for one_pathway in one_ko[1].split(','):
                one_pathway = one_pathway.replace('ko', 'map')
                try:
                    ko2ko_pathways[one_ko[0]]+='''<a href="/chlamdb/KEGG_mapp_ko/%s" target="_top">%s / %s</a></br>''' % (one_pathway,
                                                                                                           one_pathway,
                                                                                                           pathway2category[one_pathway]) # .translate(None, digits+'\.')
                except:
                    ko2ko_pathways[one_ko[0]]+='''<a href="/chlamdb/KEGG_mapp_ko/%s" target="_top">%s / %s</a></br>''' % (one_pathway,
                                                                                                           one_pathway,
                                                                                                           '?')

        if one_ko[2] != '-':
            ko2ko_modules[one_ko[0]] = ''
            for one_module in one_ko[2].split(','):
                #one_pathway = one_pathway.replace('ko', 'map')
                ko2ko_modules[one_ko[0]]+='''<a href="/chlamdb/KEGG_module_map/%s" target="_top">%s / %s</a></br>''' % (one_module,
                                                                                                           one_module,
                                                                                                           module2category[one_module])
    for ko in locus_tag2ko.values():
        if ko not in ko2ko_pathways:
            ko2ko_pathways[ko] = '-'
        if ko not in ko2ko_modules:
            ko2ko_modules[ko] = '-'

    return locus2annot, \
           locus_tag2cog_catego, \
           locus_tag2cog_name, \
           locus_tag2ko, \
           pathway2category, \
           module2category, \
           ko2ko_pathways, \
           ko2ko_modules, \
           locus2interpro



def genome_annotation(request, accession):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = f'select t1.seqfeature_id,locus_tag,start,stop,t2.name,product from annotation_seqfeature_id2locus t1 ' \
          f' inner join term t2 on t1.feature_type_id=t2.term_id ' \
          f' inner join annotation_seqfeature_id2RNA_annotation t3 on t1.seqfeature_id=t3.seqfeature_id ' \
          f' inner join bioentry t4 on t1.bioentry_id=t4.bioentry_id where t4.accession="{accession}" ' \
          f' union select t1.seqfeature_id,locus_tag,start,stop,t2.name,product from annotation_seqfeature_id2locus t1 ' \
          f' inner join term t2 on t1.feature_type_id=t2.term_id ' \
          f' inner join annotation_seqfeature_id2CDS_annotation t3 on t1.seqfeature_id=t3.seqfeature_id ' \
          f' inner join bioentry t4 on t1.bioentry_id=t4.bioentry_id ' \
          f' where t4.accession="{accession}" order by start ASC;'

    ordered_data = server.adaptor.execute_and_fetchall(sql,)

    locus_list = []
    for i in ordered_data:
        if i[4] == 'CDS':
            locus_list.append(i[1])

    biodb_id_sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb
    biodb_id = server.adaptor.execute_and_fetchall(biodb_id_sql,)[0][0]

    locus2annot, \
    locus_tag2cog_catego, \
    locus_tag2cog_name, \
    locus_tag2ko, \
    pathway2category, \
    module2category, \
    ko2ko_pathways, \
    ko2ko_modules, \
    locus2interpro = get_locus_annotations(biodb, locus_list)

    locus2annot_dico = {}
    for i in locus2annot:
        locus2annot_dico[i[2]] = i

    accession2taxon = manipulate_biosqldb.accession2taxon_id(server, biodb)
    sql = 'select bioentry_id, taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s"' % biodb
    accession2taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    '''
    series, \
    labels, \
    serie_all_counts, \
    serie_target_counts, \
    series_counts, \
    labels_counts, \
    category_description, \
    category_map, \
    n_missing_cog, \
    missing_cog_list = locus_tag2cog_series(biodb, locus_list, reference_taxon=None)
    '''

    return render(request, 'chlamdb/genome_annotation.html', my_locals(locals()))



def blastnr_cat_info(request, accession, rank, taxon):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    target_accessions = [i for i in request.GET.getlist('h')]
    counttype = request.GET.getlist('t')[0]
    top_n = request.GET.getlist('n')[0]

    biodb_id_sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb
    biodb_id = server.adaptor.execute_and_fetchall(biodb_id_sql,)[0][0]

    if counttype == 'Majority':
        sql = 'select B.locus_tag, A.%s ,A.n from (select seqfeature_id,%s, count(*) as n from blastnr_blastnr A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id where hit_number<=%s and query_bioentry_id=%s ' \
              ' group by seqfeature_id, %s order by seqfeature_id,n DESC) A ' \
              ' inner join custom_tables_locus2seqfeature_id B on A.seqfeature_id=B.seqfeature_id' % (rank,
                                                                                                      rank,
                                                                                                      top_n,
                                                                                                      accession,
                                                                                                      rank)

        data = server.adaptor.execute_and_fetchall(sql,)
        category2count = {}
        all_query_locus_list = []
        majority_locus_list = []

        for i in data:
            # keep only the majoritary taxon
            if i[0] not in all_query_locus_list:
                # keep only data for taxon of interest
                if i[1] == taxon:
                    majority_locus_list.append(i[0])
                all_query_locus_list.append(i[0])
        locus_list = majority_locus_list

    elif counttype == 'BBH':
        sql = ' select locus_tag from (select t2.*,t1.locus_tag from custom_tables_locus2seqfeature_id t1 ' \
              ' inner join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id' \
              ' where hit_number=1) A inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id ' \
              ' where %s="%s"  and query_bioentry_id=%s;' % (rank, taxon, accession)

        locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    else:
        raise 'invalide type'

    # get number of hits for each kingdom
    sql_superkingdom = 'select locus_tag,superkingdom, count(*) as n from' \
    ' (select t2.*,t1.locus_tag from custom_tables_locus2seqfeature_id t1 ' \
    ' inner join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id' \
    ' where locus_tag in ("%s")) A' \
    ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id  group by locus_tag,superkingdom;' % (  '","'.join(locus_list))

    superkingdom_data = server.adaptor.execute_and_fetchall(sql_superkingdom,)
    locus2superkingdom_counts = {}
    for row in superkingdom_data:
        if row[0] not in locus2superkingdom_counts:
            locus2superkingdom_counts[row[0]] = {}
            locus2superkingdom_counts[row[0]][row[1]] = row[2]
        else:
            locus2superkingdom_counts[row[0]][row[1]] = row[2]
    # calculate sum of hits
    for locus_tag in locus2superkingdom_counts:
        total = sum([int(i) for i in locus2superkingdom_counts[locus_tag].values()])
        lst = []
        for superkingdom in locus2superkingdom_counts[locus_tag]:
            lst.append(["%s_percent" % superkingdom, round(float(locus2superkingdom_counts[locus_tag][superkingdom])/total, 2)])
        for kd in lst:
            locus2superkingdom_counts[locus_tag][kd[0]] = kd[1]
        locus2superkingdom_counts[locus_tag]["TOTAL"] = total


    locus2annot, \
    locus_tag2cog_catego, \
    locus_tag2cog_name, \
    locus_tag2ko, \
    pathway2category, \
    module2category, \
    ko2ko_pathways, \
    ko2ko_modules, \
    locus2interpro = get_locus_annotations( locus_list)

    accession2taxon = manipulate_biosqldb.accession2taxon_id(server, biodb)
    sql = 'select bioentry_id, taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s"' % biodb
    accession2taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    circos_url = '?ref=%s&' % locus2annot[0][-1]
    target_taxons = [str(accession2taxon[i]) for i in target_accessions]
    reference_taxon = str(accession2taxon[accession])
    target_taxons.pop(target_taxons.index(reference_taxon))
    circos_url += "t="+('&t=').join((target_taxons)) + '&h=' + ('&h=').join(locus_list)

    locus_filter = '"' + '","'.join(locus_list) + '"'
    sql = 'select locus_tag,subject_accession,subject_kingdom,subject_scientific_name,subject_taxid,evalue,percent_identity ' \
          ' from custom_tables. locus2seqfeature_id_%s t1 ' \
          ' inner join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' where t1.locus_tag in (%s) and hit_number=1;' % (  locus_filter)

    locus_tag2blastnr_BBH = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    '''
    series, \
    labels, \
    serie_all_counts, \
    serie_target_counts, \
    series_counts, \
    labels_counts, \
    category_description, \
    category_map, \
    n_missing_cog, \
    missing_cog_list = locus_tag2cog_series( locus_list, reference_taxon=None)
    '''


    return render(request, 'chlamdb/blastnr_info.html', my_locals(locals()))

def identity_heatmap(request):
    biodb = settings.BIODB
    import numpy as np
    from chlamdb.plots import pairwiseid_plots

    server, db = manipulate_biosqldb.load_db(biodb)

    form_class = heatmap_form(biodb)

    if request.method == 'POST':

        form_venn = form_class(request.POST)

        if 'venn' in request.POST and form_venn.is_valid():
            taxon_list = form_venn.cleaned_data['targets']

            if len(taxon_list) < 3:
                wrong_count = True
                return render(request, 'chlamdb/identity_heatmap.html', my_locals(locals()))

            plot_type = form_venn.cleaned_data['plot']
            taxon_filter = '"'+'","'.join(taxon_list)+'"'


            if plot_type == 'blast_identity':
                sql2 = 'select taxon_1, taxon_2, median_identity from comparative_tables_reciprocal_BBH_average_identity ' \
                       'where taxon_1 in (%s) and taxon_2 in (%s);' % (taxon_filter,
                                                                       taxon_filter)
            if plot_type == 'core_msa_identity':
                sql2 = 'select taxon_1, taxon_2, identity from comparative_tables_core_orthogroups_identity_msa ' \
                       'where taxon_1 in (%s) and taxon_2 in (%s);' % (taxon_filter, taxon_filter)
            if plot_type == 'n_RBBH':
                sql2 = 'select taxon_1, taxon_2, n_pairs from comparative_tables_reciprocal_BBH_average_identity ' \
                      ' where taxon_1 in (%s) and taxon_2 in (%s) UNION select taxon_1, taxon_2, n_pairs ' \
                      ' from comparative_tables_reciprocal_BBH_average_identity' \
                      ' where taxon_2 in (%s) and taxon_1 in (%s)' % (taxon_filter,
                                                                      taxon_filter,taxon_filter,
                                                                      taxon_filter)
            if plot_type == 'n_shared_orthogroups':
                sql2 = 'select taxon_1, taxon_2, n_shared_orthogroups from comparative_tables_shared_orthogroups ' \
                       'where taxon_1 in (%s) and taxon_2 in (%s);' % (taxon_filter,
                                                                      taxon_filter)
            data = server.adaptor.execute_and_fetchall(sql2,)
            m = np.empty([len(taxon_list), len(taxon_list)], dtype=float)
            for row in data:
                index1 = taxon_list.index(str(row[0]))
                index2 = taxon_list.index(str(row[1]))

                m[index1][index2] = float(row[2])
                m[index2][index1] = float(row[2])
            for n in range(0, len(taxon_list)):
                m[n][n] = None

            sql3 = 'select taxon_id, t2.description from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
                   ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % (biodb)

            taxon_id2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

            labels = [taxon_id2description[i] for i in taxon_list]
            if len(labels) < 9:
                small = True
            path = settings.BASE_DIR + '/assets/temp/heatmap.svg'
            asset_path = '/temp/heatmap.svg'
            pairwiseid_plots.identity_heatmap_plot(m, output_path=path, labels=labels, reverse=True)

            envoi = True

    else:
        form_venn = form_class()
    return render(request, 'chlamdb/identity_heatmap.html', my_locals(locals()))


# TODO: implement with plotly to avoid having to save a picture at 
# every request + weird errors regarding Qt
def pan_genome(request, type):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)

    venn_form_class = make_venn_from(db, plasmid=False)

    if request.method != "POST":
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    form = venn_form_class(request.POST)
    if not form.is_valid():
        # should add an error message
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    import seaborn as sn
    import matplotlib.pyplot as plt
    import pandas as pd
    import time

    taxids = form.get_taxids()
    print( "taxids",taxids)

    if type == "COG":
        df_hits = db.get_cog_hits(taxids, search_on="taxid", indexing="taxid")
        type_txt = "COG"
    elif type == "orthology":
        df_hits = db.get_og_count(taxids, search_on="taxid")
        type_txt = "orthologs"
    elif type == "ko":
        df_hits = db.get_ko_hits(taxids, search_on="taxid")
        type_txt = "KO"
    else:
        # should add an error message
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    target2description = db.get_genomes_description().description.to_dict()
    df_hits.columns = [target2description[i] for i in df_hits.columns.values]
    print(" df_hits.columns",  df_hits )
    path2 = settings.BASE_DIR + '/assets/temp/pangenome_barplot.svg'
    asset_path2 = '/temp/pangenome_barplot.svg'
    n_entries_per_genome = df_hits[df_hits>0].count(axis=0).sort_values()

    fig, ax = plt.subplots()
    barplot = ax.bar(list(n_entries_per_genome.index), list(n_entries_per_genome.values))
    ax.set_ylabel(f"{type} count")
    ax.set_xticklabels(list(n_entries_per_genome.index), rotation=45, horizontalalignment="right")
    fig.tight_layout()
    fig.savefig(path2)

    hsh_shared_count = {}
    hsh_total_count = {}
    for entry, pres in df_hits[n_entries_per_genome.index[0]].items():
        if pres==0:
            continue
        hsh_shared_count[entry] = 1
        hsh_total_count[entry] = 1

    total_entries = [n_entries_per_genome.values[0]]
    shared_entries = []
    for i in range(1, len(n_entries_per_genome)):
        cur_total = total_entries[-1]
        cur_shared = 0
        for entry, pres in df_hits[n_entries_per_genome.index[i]].items():
            if pres==0:
                continue
            v = hsh_shared_count.get(entry, 0)
            if v==i:
                cur_shared += 1
                hsh_shared_count[entry] = i+1
            if entry not in hsh_total_count:
                cur_total += 1
                hsh_total_count[entry] = 1
        total_entries.append(cur_total)
        shared_entries.append(cur_shared)

    path = settings.BASE_DIR + '/assets/temp/pangenome.svg'
    asset_path = '/temp/pangenome.svg'
    fig, ax = plt.subplots()

    ax.plot(total_entries)
    ax2 = ax.twinx()
    ax2.plot([i for i in range(1, len(n_entries_per_genome))], shared_entries, color="red")

    ax.set_xticks([i for i in range(0, len(n_entries_per_genome))])
    ax.set_xticklabels(n_entries_per_genome.index, rotation=45, horizontalalignment="right")
    ax2.set_ylabel(f"Number of shared {type_txt}")
    ax.set_ylabel(f"Number of {type_txt} in pangenome")

    fig.tight_layout()
    fig.savefig(path)

    # Serie with the number of genomes having a given cog
    missing_entries = df_hits[df_hits == 0].count(axis=1)
    core = len(missing_entries[missing_entries == 0])

    # cogs present in all but one genome
    core_minus1 = len(missing_entries[missing_entries == 1])
    total = len(df_hits.index)

    genome_count_list = []
    for i, count in enumerate(n_entries_per_genome):
        genome_count_list.append([i+1, count])
    envoi = True
    return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))


def core_genome_missing(request, type):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    venn_form_class = make_venn_from(biodb)

    if request.method == 'POST':
        form = venn_form_class(request.POST)
        if form.is_valid():
            from chlamdb.plots import core_pan_genome_plots
            import numpy
            import pandas
            from chlamdb.phylo_tree_display import ete_motifs
            from chlamdb.biosqldb import biosql_own_sql_tables


            taxon_list = form.cleaned_data['targets']

            filter = '`'+'`,`'.join(taxon_list)+'`'

            sql = 'select  orthogroup from comparative_tables_%s' % (type)
            group_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


            sql = 'select %s from comparative_tables_%s' % (filter, type)

            data = numpy.array([list(i) for i in server.adaptor.execute_and_fetchall(sql,)])

            count_df = pandas.DataFrame(data, columns=taxon_list, index=group_list)

            groups_with_paralogs = count_df[(count_df > 1).sum(axis=1) > 0].index
            count_df = count_df.drop(groups_with_paralogs)

            t0 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)]
            t1 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)-1]
            t2 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)-2]
            t3 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)-3]
            t4 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)-4]
            t5 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)-5]
            t6 = count_df[(count_df == 1).sum(axis=1) >= len(taxon_list)-6]

            core0 = t0.index.tolist()
            core1 = t1.index.tolist()
            core2 = t2.index.tolist()
            core3 = t3.index.tolist()
            core4 = t4.index.tolist()
            core5 = t5.index.tolist()
            core6 = t6.index.tolist()

            count_list = t1.sum(axis=1)
            group2count = zip(count_list.index.tolist(), list(count_list))

            taxon2n_missing1 = {}
            taxon2n_missing2 = {}
            taxon2n_missing3 = {}
            taxon2n_missing4 = {}
            taxon2n_missing5 = {}
            taxon2n_missing6 = {}

            for group in core1:
                for taxon in taxon_list:
                    if count_df.loc[group,taxon] != 1:
                        if taxon not in taxon2n_missing1:
                            taxon2n_missing1[taxon] = 1
                        else:
                            taxon2n_missing1[taxon] += 1

            for group in core2:
                for taxon in taxon_list:
                    if count_df.loc[group,taxon] != 1:
                        if taxon not in taxon2n_missing2:
                            taxon2n_missing2[taxon] = 1
                        else:
                            taxon2n_missing2[taxon] += 1
            for group in core3:
                for taxon in taxon_list:
                    if count_df.loc[group, taxon] != 1:
                        if taxon not in taxon2n_missing3:
                            taxon2n_missing3[taxon] = 1
                        else:
                            taxon2n_missing3[taxon] += 1
            for group in core4:
                for taxon in taxon_list:
                    if count_df.loc[group, taxon] != 1:
                        if taxon not in taxon2n_missing4:
                            taxon2n_missing4[taxon] = 1
                        else:
                            taxon2n_missing4[taxon] += 1
            for group in core5:
                for taxon in taxon_list:
                    if count_df.loc[group,taxon] != 1:
                        if taxon not in taxon2n_missing5:
                            taxon2n_missing5[taxon] = 1
                        else:
                            taxon2n_missing5[taxon] += 1
            for group in core6:
                for taxon in taxon_list:
                    if count_df.loc[group,taxon] != 1:
                        if taxon not in taxon2n_missing6:
                            taxon2n_missing6[taxon] = 1
                        else:
                            taxon2n_missing6[taxon] += 1

            # group2taxon2count

            n_missing2taxon2count = {}
            n_missing2taxon2count['1'] = taxon2n_missing1
            n_missing2taxon2count['2'] = taxon2n_missing2
            n_missing2taxon2count['3'] = taxon2n_missing3
            n_missing2taxon2count['4'] = taxon2n_missing4
            n_missing2taxon2count['5'] = taxon2n_missing5
            n_missing2taxon2count['6'] = taxon2n_missing6

            labels = ['1','2','3','4', '5', '6']
            tree1, style1 = ete_motifs.multiple_profiles_heatmap(labels,n_missing2taxon2count, column_scale=True, as_float=False)

            path2 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
            asset_path2 = '/temp/ortho_tree1.svg'
            tree1.render(path2, dpi=800, tree_style=style1)

            envoi = True

            envoi = True
    else:  
        form = venn_form_class()
    return render(request, 'chlamdb/core_genome_missing.html', my_locals(locals()))





def pairwiseid(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    pairwiseid_form_class = make_pairwiseid_form(biodb)

    if request.method == 'POST':
        from chlamdb.plots import pairwiseid_plots
        from chlamdb.biosqldb import biosql_own_sql_tables
        from chlamdb.phylo_tree_display import ete_motifs

        form = pairwiseid_form_class(request.POST)

        if form.is_valid():
            plot_type = form.cleaned_data['plot']
            genome_1 = form.cleaned_data['genome_1']
            genome_2 = form.cleaned_data['genome_2']
            genome_3 = form.cleaned_data['genome_3']
            genome_4 = form.cleaned_data['genome_4']
            genome_5 = form.cleaned_data['genome_5']
            genome_6 = form.cleaned_data['genome_6']

            taxid2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            sql = 'select blast_identity_a_vs_b from comparative_tables_reciprocal_BBH where taxon_1 in (%s,%s) ' \
                  'and taxon_2 in (%s,%s);' % (genome_1,
                                                genome_2,
                                                genome_1,
                                                genome_2)

            data1 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            data_list = [data1]
            label_list = ["%s" %(taxid2description[genome_2])]
            if genome_3 != 'None':
                sql = 'select blast_identity_a_vs_b from comparative_tables_reciprocal_BBH where taxon_1 in (%s,%s) ' \
                      'and taxon_2 in (%s,%s);' % (genome_1,
                                                    genome_3,
                                                    genome_1,
                                                    genome_3)

                data2 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data2)
                label_list.append("%s" %(taxid2description[genome_3]))

            if genome_4 != 'None':
                sql = 'select blast_identity_a_vs_b from comparative_tables_reciprocal_BBH where taxon_1 in (%s,%s) ' \
                      'and taxon_2 in (%s,%s);' % (
                                    genome_1,
                                    genome_4,
                                    genome_1,
                                    genome_4)

                data3 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data3)
                label_list.append("%s" %(taxid2description[genome_4]))

            if genome_5 != 'None':
                sql = 'select blast_identity_a_vs_b from comparative_tables_reciprocal_BBH where taxon_1 in (%s,%s) ' \
                      'and taxon_2 in (%s,%s);' % (genome_1,
                                                    genome_5,
                                                    genome_1,
                                                    genome_5)

                data3 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data3)
                label_list.append("%s" %(taxid2description[genome_5]))

            if genome_6 != 'None':
                sql = 'select blast_identity_a_vs_b from comparative_tables_reciprocal_BBH where taxon_1 in (%s,%s) ' \
                      'and taxon_2 in (%s,%s);' % (genome_1,
                                                    genome_6,
                                                    genome_1,
                                                    genome_6)

                data3 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data3)
                label_list.append("%s" %(taxid2description[genome_6]))

            path1 = settings.BASE_DIR + '/assets/temp/plot.svg'

            asset_path1 = '/temp/plot.svg'
            pairwiseid_plots.density_plot(data_list,
                                          label_list,
                                          header=taxid2description[genome_1],
                                          xlab="identity (%)",
                                          ylab="density",
                                          output_path=path1,
                                          min_value=10,
                                          max_value=100,
                                          show_mode=True,
                                          show_median=True
                                          )

            sql = 'SELECT orthogroup, count(*) as n FROM (select  orthogroup,taxon_id from orthology_detail ' \
                  ' group by orthogroup,taxon_id) A  GROUP BY orthogroup' % biodb
            group2n_organisms = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select orthogroup_1,blast_identity_a_vs_b from (select seqfeature_id_1,blast_identity_a_vs_b,orthogroup_1 ' \
                  ' from comparative_tables_reciprocal_BBH where orthogroup_1=orthogroup_2 and taxon_1 in (%s,%s) ' \
                  ' and taxon_2 in (%s,%s)) A inner join custom_tables_locus2seqfeature_id t2 ' \
                  ' on A.seqfeature_id_1=t2.seqfeature_id;' % (genome_1,
                                                               genome_2,
                                                               genome_1,
                                                               genome_2,
                                                               biodb)
            identity_list = []
            genome_count_list = []
            data = server.adaptor.execute_and_fetchall(sql,)
            for row in data:
                identity_list.append(row[1])

                genome_count_list.append(group2n_organisms[row[0]])

            sql = 'select taxon_2, median_identity, average_identity,n_pairs from comparative_tables_reciprocal_BBH_average_identity ' \
                  ' where taxon_1=%s UNION select taxon_1, median_identity, average_identity,n_pairs ' \
                  ' from comparative_tables_reciprocal_BBH_average_identity' \
                  ' where taxon_2=%s' % (genome_1,genome_1)
            data = server.adaptor.execute_and_fetchall(sql,)

            sql2 = 'select taxon_2, identity from comparative_tables_core_orthogroups_identity_msa ' \
                   'where taxon_1=%s;' % (genome_1)
            #   taxon2core_msa_identity = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

            sql3 = 'select taxon_2, n_shared_orthogroups from comparative_tables_shared_orthogroups ' \
                   'where taxon_1=%s;' % (genome_1)

            taxon2n_shared_orthogroups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

            sql = 'select SUM(n_CDS) from biodatabase t1 ' \
                  ' inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                  ' inner join genomes_info t3 on t2.ACCESSION=t3.ACCESSION ' \
                  ' where t1.name="%s" and taxon_id=%s group by taxon_id;' % (genome_1)
            n_CDS = server.adaptor.execute_and_fetchall(sql,)[0][0]

            sql = 'select count(*) from (select taxon_id from orthology_detail ' \
                  ' WHERE taxon_id=%s group by orthogroup) A;' % (genome_1)
            n_orthogroups = server.adaptor.execute_and_fetchall(sql,)[0][0]

            taxon2category2id = {}
            taxon2category2id['mean'] = {}
            taxon2category2id['median'] = {}
            taxon2category2id['n_RBBH'] = {}
            taxon2category2id['n_RBBH'][str(genome_1)] = n_CDS
            taxon2category2id['n_shared_orthogroups'] = {}
            taxon2category2id['n_shared_orthogroups'][str(genome_1)] = n_orthogroups
            #taxon2category2id['core_align'] = {}

            for row in data:
                taxon2category2id['mean'][str(row[0])] = row[2]
                taxon2category2id['median'][str(row[0])] = row[1]
                taxon2category2id['n_RBBH'][str(row[0])] = int(row[3])
                taxon2category2id['n_shared_orthogroups'][str(row[0])] = int(taxon2n_shared_orthogroups[str(row[0])])
                #taxon2category2id['core_align'][str(row[0])] = taxon2core_msa_identity[str(row[0])]

            tree1, style1 = ete_motifs.multiple_profiles_heatmap(
                                                        column_labels=['mean', 'median', 'n_RBBH', 'n_shared_orthogroups'],
                                                        group2taxon2count=taxon2category2id,
                                                        identity_scale=False,
                                                        show_labels=True,
                                                        column_scale=True,
                                                                 reference_taxon=genome_1)

            path2 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
            asset_path2 = '/temp/ortho_tree1.svg'
            tree1.render(path2, dpi=800, tree_style=style1)
            envoi = True
    else:  
        form = pairwiseid_form_class()
    return render(request, 'chlamdb/pairwise_id.html', my_locals(locals()))

def multiple_codon_usage(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    pairwiseCDS_length_form_class = make_pairwiseCDS_length_form(biodb)

    if request.method == 'POST':
        from chlamdb.plots import pairwiseid_plots
        from chlamdb.biosqldb import biosql_own_sql_tables
        from chlamdb.phylo_tree_display import ete_motifs
        import numpy

        form = pairwiseCDS_length_form_class(request.POST)

        genome_median_cds_length = []

        if form.is_valid():
            from chlamdb.plots import pca_seq_composition

            genome_1 = form.cleaned_data['genome_1']
            genome_2 = form.cleaned_data['genome_2']
            genome_3 = form.cleaned_data['genome_3']
            genome_4 = form.cleaned_data['genome_4']

            taxid2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            taxon_id_list = [genome_1]

            if genome_2 != 'None':
                taxon_id_list.append(genome_2)

            if genome_3 != 'None':
                taxon_id_list.append(genome_3)
            if genome_4 != 'None':
                taxon_id_list.append(genome_4)

            taxon_id_filrter = '"'+'","'.join(taxon_id_list)+'"'
            sql2 = 'select t2.description, t1.* from custom_tables_codon_usage_percent t1  inner join bioentry t2 on t1.taxon_id=t2.taxon_id' \
                   ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                   ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" and t1.taxon_id in (%s)' % (taxon_id_filrter)
            data = [i for i in server.adaptor.execute_and_fetchall(sql2,)]


            mat = numpy.array([[list(i)[0]] + list(i)[4:68] for i in data])
            path = settings.BASE_DIR + '/assets/temp/hydro.png'
            asset_path = '/temp/hydro.png'
            pca_seq_composition.multiple_aa_composition_pca(mat, path)

            envoi = True
    else:  
        form = pairwiseCDS_length_form_class()
    return render(request, 'chlamdb/codons_multiple_pca.html', my_locals(locals()))

def get_BBH_non_chlamydiae_taxonomy(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)


    acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
    filter = '"' + '","'.join(acc_list) + '"'
    sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
    taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]

    filter = ','.join(taxon_list)

    sql = 'select seqfeature_id,query_taxon_id,hit_number,subject_taxid,t2.superkingdom, t2.phylum ' \
          ' from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 on t1.subject_taxid=taxon_id ' \
          ' where t1.query_taxon_id!=127 and phylum!="Chlamydiae" order by hit_number'

    sql = 'select seqfeature_id,query_taxon_id,hit_number,subject_taxid,' \
          ' t2.superkingdom, t2.kingdom,t2.order, t2.phylum ' \
          ' from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 on t1.subject_taxid=taxon_id ' \
          ' where t1.query_taxon_id in (%s) and phylum!="Chlamydiae" order by hit_number' % (filter)

    data = server.adaptor.execute_and_fetchall(sql,)

    '''
    phylum2count = {}
    seqfeature_ok = []
    for n,row in enumerate(data):
        if n % 1000 == 0:
            print "%s / %s" % (n, len(data))
        if row[0] in seqfeature_ok:
            continue
        else:
            if row[5] not in phylum2count:
                phylum2count[row[5]] = 1
            else:
                phylum2count[row[5]] += 1
            seqfeature_ok.append(row[0])
    '''
    superkingdom2count = {}
    kingdom2count = {}
    order2count = {}
    seqfeature_ok = []
    phylum2count= {}
    for n,row in enumerate(data):
        if n % 1000 == 0:
            print ("%s / %s" % (n, len(data)))
        if row[0] in seqfeature_ok:
            continue
        else:
            if row[4] not in superkingdom2count:
                superkingdom2count[row[4]] = 1
            else:
                superkingdom2count[row[4]] += 1
            if row[5] not in kingdom2count:
                kingdom2count[row[5]] = 1
            else:
                kingdom2count[row[5]] += 1
            if row[6] not in order2count:
                order2count[row[6]] = 1
            else:
                order2count[row[6]] += 1

            if row[7] not in phylum2count:
                phylum2count[row[7]] = 1
            else:
                phylum2count[row[7]] += 1

            seqfeature_ok.append(row[0])

    with open('/home/trestan/superkingdom_count.tab', 'w') as f:
        for superkingdom in superkingdom2count:
            f.write('%s\t%s\n' % (superkingdom, str(superkingdom2count[superkingdom])))
    with open('/home/trestan/kingdom_count.tab', 'w') as f:
        for kingdom in kingdom2count:
            f.write('%s\t%s\n' % (kingdom, str(kingdom2count[kingdom])))
    with open('/home/trestan/order_count.tab', 'w') as f:
        for order in order2count:
            f.write('%s\t%s\n' % (order, str(order2count[order])))
    with open('/home/trestan/phylum_count.tab', 'w') as f:
        for phylum in phylum2count:
            f.write('%s\t%s\n' % (phylum, str(phylum2count[phylum])))

def multipleGC(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    pairwiseCDS_length_form_class = make_pairwiseCDS_length_form(biodb)

    if request.method == 'POST':
        from chlamdb.plots import pairwiseid_plots
        from chlamdb.biosqldb import biosql_own_sql_tables
        from chlamdb.phylo_tree_display import ete_motifs
        import numpy

        form = pairwiseCDS_length_form_class(request.POST)

        genome_median_cds_length = []

        if form.is_valid():
            from chlamdb.plots import pca_seq_composition

            genome_1 = form.cleaned_data['genome_1']
            genome_2 = form.cleaned_data['genome_2']
            genome_3 = form.cleaned_data['genome_3']
            genome_4 = form.cleaned_data['genome_4']

            taxid2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            taxon_id_list = [genome_1]

            if genome_2 != 'None':
                taxon_id_list.append(genome_2)

            if genome_3 != 'None':
                taxon_id_list.append(genome_3)
            if genome_4 != 'None':
                taxon_id_list.append(genome_4)

            taxon_id_filrter = '"'+'","'.join(taxon_id_list)+'"'
            sql2 = 'select t2.description, t1.* from custom_tables_aa_usage_count t1  inner join bioentry t2 on t1.taxon_id=t2.taxon_id' \
                   ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                   ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" and t1.taxon_id in (%s)' % (biodb, taxon_id_filrter)
            data = [i for i in server.adaptor.execute_and_fetchall(sql2,)]


            mat = numpy.array([[list(i)[0]] + list(i)[4:23] for i in data])
            path = settings.BASE_DIR + '/assets/temp/hydro.png'
            asset_path = '/temp/hydro.png'
            pca_seq_composition.multiple_aa_composition_pca(mat, path)

            envoi = True
    else:  
        form = pairwiseCDS_length_form_class()
    return render(request, 'chlamdb/aa_multiple_pca.html', my_locals(locals()))

def pairwiseCDS_length(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    pairwiseCDS_length_form_class = make_pairwiseCDS_length_form(biodb)

    if request.method == 'POST':
        from chlamdb.plots import pairwiseid_plots
        from chlamdb.biosqldb import biosql_own_sql_tables
        from chlamdb.phylo_tree_display import ete_motifs
        import numpy

        form = pairwiseCDS_length_form_class(request.POST)

        genome_median_cds_length = []

        if form.is_valid():
            genome_1 = form.cleaned_data['genome_1']
            genome_2 = form.cleaned_data['genome_2']
            genome_3 = form.cleaned_data['genome_3']
            genome_4 = form.cleaned_data['genome_4']

            taxid2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
            if db_driver == 'mysql':
                sql = 'select CHAR_LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_1)
            if db_driver == 'sqlite':
                sql = 'select LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_1)

            data1 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            genome_median_cds_length.append([taxid2description[genome_1],
                                             numpy.mean(data1),
                                             round((sum(i <= 400 for i in data1)/float(len(data1)))*100,2)])

            data_list = [data1]
            label_list = ["%s" % taxid2description[genome_1]]
            if genome_2 != 'None':
                if db_driver == 'mysql':
                    sql = 'select CHAR_LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_2)
                if db_driver == 'sqlite':
                    sql = 'select LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_2)

                data2 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data2)
                label_list.append("%s" % taxid2description[genome_2])
                genome_median_cds_length.append([taxid2description[genome_2], numpy.mean(data2),
                                                 round((sum(i <= 400 for i in data2)/float(len(data2)))*100,2)])

            if genome_3 != 'None':
                sql = 'select CHAR_LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_3)

                data3 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data3)
                label_list.append("%s" % taxid2description[genome_3])
                genome_median_cds_length.append([taxid2description[genome_3],
                                                 numpy.mean(data3),
                                                 round((sum(i <= 400 for i in data3)/float(len(data3)))*100,2)])
            if genome_4 != 'None':
                if db_driver == 'mysql':
                    sql = 'select CHAR_LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_4)
                if db_driver == 'sqlite':
                    sql = 'select LENGTH(translation)*3 from orthology_detail where taxon_id =%s;' % (genome_4)

                data4 = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                data_list.append(data4)
                label_list.append("%s" % taxid2description[genome_4])
                genome_median_cds_length.append([taxid2description[genome_4],
                                                 numpy.mean(data4),
                                                 round((sum(i <= 400 for i in data4)/float(len(data4)))*100,2)])


            m = 0
            for one_list in data_list:
                if max(one_list) > m:
                    m=max(one_list)
            path1 = settings.BASE_DIR + '/assets/temp/plot.svg'

            asset_path1 = '/temp/plot.svg'
            pairwiseid_plots.density_plot(data_list,
                                          label_list,
                                          output_path=path1,
                                          show_median=False,
                                          min_value=0,
                                          max_value=m,
                                          header="Distribution of CDS length",
                                          xlab="length (bp)",
                                          ylab="density",
                                          breaks_manual=m/200)


            envoi = True
    else:  
        form = pairwiseCDS_length_form_class()
    return render(request, 'chlamdb/pairwise_CDS_length.html', my_locals(locals()))





def blastnr_euk(request):
    biodb = settings.BIODB

    from chlamdb.phylo_tree_display import phylo_tree_bar
    from chlamdb.plots import pairwiseid_plots

    server, db = manipulate_biosqldb.load_db(biodb)


    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;'

    data = server.adaptor.execute_and_fetchall(sql,)

    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]

    sql = 'select seqfeature_id, locus_tag from custom_tables_locus2seqfeature_id'

    seqfeature_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag,n_species from custom_tables_seqfeature_id2n_species t1 ' \
          ' inner join custom_tables_locus2seqfeature_id t2 on t1.seqfeature_id=t2.seqfeature_id;'

    locus_tag2n_species = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 ' \
               ' on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    sql_taxon = 'select taxon_id,t2.description from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % biodb
                
    taxon_id2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_taxon,))

    taxon_list = taxon_id2description.keys()

    sql_genome_size = 'select taxon_id, n_CDS from biodatabase t1 ' \
                      ' inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                      ' inner join genomes_info t3 on t2.accession=t3.accession ' \
                      ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%";' % (biodb)

    taxon_id2n_CDS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_genome_size,))

    '''
    sql_best_hit_euk = 'select t1.query_taxon_id, count(*) from blastnr_blastnr t1 ' \
                       ' inner join blastnr_blastnr_taxonomy t2 on t1.subject_taxid=t2.taxon_id ' \
                       ' where t2.superkingdom = "Eukaryota" and t1.hit_number=1 group by t1.query_taxon_id;'
    '''
    sql_any_hit_euk = 'select A.query_taxon_id, count(*) from (select t1.query_taxon_id, t1.seqfeature_id ' \
                      ' from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 on t1.subject_taxid=t2.taxon_id ' \
                      ' where t2.superkingdom = "Eukaryota" group by t1.query_taxon_id, t1.seqfeature_id) A ' \
                      ' group by A.query_taxon_id;'

    #taxon_id2n_best_hits_euk = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_best_hit_euk,))

    '''
    print 'taxon_id2any_hit_euk...'
    taxon_id2n_any_hits_euk = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_any_hit_euk,))
    '''
    sql_best_non_self_euk = 'select query_taxon_id, count(*) as n from blastnr_blastnr_best_non_self_phylum ' \
                      ' where superkingdom="Eukaryota" group by query_taxon_id;'
    '''
    sql_best_non_self_euk_50 = 'select query_taxon_id, count(*) as n from blastnr_blastnr_best_non_self_phylum ' \
                      ' where superkingdom="Eukaryota" and percent_identity>=50 group by query_taxon_id;'
    '''

    sql_best_non_self_euk_species_specific = 'select query_taxon_id, count(*) as n from blastnr_blastnr_best_non_self_phylum t1' \
                                             ' inner join custom_tables_seqfeature_id2n_species t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                                             ' where superkingdom="Eukaryota" and n_species=1 group by query_taxon_id;'

    taxon_id2n_best_non_chlamydial_euk = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_best_non_self_euk,))

    '''
    taxon_id2n_best_non_chlamydial_euk_50 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_best_non_self_euk_50,))
    '''

    taxon_id2n_species_specific = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_best_non_self_euk_species_specific,))

    taxon2values = {}
    set2taxon2values = {}
    taxon_id2n_any_hits_euk_perc = {}
    taxon_id2n_best_non_chlamydial_euk_perc = {}
    for taxon in taxon_list:
        '''
        try:
            n_best_euk = taxon_id2n_best_hits_euk[taxon]
        except:
            n_best_euk = 0

        try:
            n_euk = taxon_id2n_any_hits_euk[taxon]
            percent_euk = round((taxon_id2n_any_hits_euk[taxon]/float(taxon_id2n_CDS[taxon]))*100,2)
            taxon_id2n_any_hits_euk_perc[taxon] = percent_euk
        except:
            n_euk = 0
            percent_euk = 0
            taxon_id2n_any_hits_euk_perc[taxon] = 0

        try:
            perc_best_euk = round((taxon_id2n_best_non_chlamydial_euk[taxon]/float(taxon_id2n_CDS[taxon]))*100,2)
        except:
            perc_best_euk = 0
            taxon_id2n_best_non_chlamydial_euk[taxon] = 0

        try:
            print taxon_id2n_best_non_chlamydial_euk_50[taxon]
        except:
            taxon_id2n_best_non_chlamydial_euk_50[taxon] = 0
        '''
        try:
            print (taxon_id2n_species_specific[taxon])
        except:
            taxon_id2n_species_specific[taxon] = 0
        try:
            taxon_id2n_best_non_chlamydial_euk_perc[taxon] = round((taxon_id2n_best_non_chlamydial_euk[taxon]/float(taxon_id2n_CDS[taxon]))*100,2)
        except:

            taxon_id2n_best_non_chlamydial_euk_perc[taxon] = 0
        try:
            print (taxon_id2n_best_non_chlamydial_euk[taxon])
        except:
            taxon_id2n_best_non_chlamydial_euk[taxon] = 0

        '''
        taxon2values[taxon] = [n_euk,
                               percent_euk,
                               taxon_id2n_best_non_chlamydial_euk[taxon],
                               perc_best_euk,
                               taxon_id2n_best_non_chlamydial_euk_50[taxon],
                               taxon_id2n_species_specific[taxon]]
        print taxon2values[taxon]
        '''
    #set2taxon2values['n any euk'] = taxon_id2n_any_hits_euk
    #set2taxon2values['percent euk'] = taxon_id2n_any_hits_euk_perc
    set2taxon2values['best other other phylum euk'] = taxon_id2n_best_non_chlamydial_euk
    #set2taxon2values['best other other phylum euk (%)'] = taxon_id2n_best_non_chlamydial_euk

    #set2taxon2values['n_best_non_self_euk_50'] = taxon_id2n_best_non_chlamydial_euk_50
    set2taxon2values['best other other phylum euk (%)'] = taxon_id2n_best_non_chlamydial_euk_perc
    set2taxon2values['n best_non_self species specific'] = taxon_id2n_species_specific

    '''
    tree1, style1 = phylo_tree_bar.plot_tree_barplot(tree,
                                                    taxon2values,
                                                    ["n any euk",
                                                     'percent euk',
                                                     'n_best_non_self_euk',
                                                     'perc_best_non_self_euk',
                                                     'n_best_non_self_euk_50',
                                                     'n specific'],
                                                    taxon2set2value_heatmap=False,
                                                    header_list2=False,
                                                    biodb=biodb,
                                                    general_max=False)
    '''

    tree1, style1 = phylo_tree_bar.plot_tree_stacked_barplot(tree,
                                 taxon2value_list_barplot=False,
                                 header_list=False, # header stackedbarplots
                                 taxon2set2value_heatmap=False,
                                 taxon2label=taxon2most_freq_phylum,
                                 header_list2=False, # header counts columns
                                 biodb=biodb,
                                 column_scale=True,
                                 general_max=False,
                                 header_list3 =['best other other phylum euk',
                                                'best other other phylum euk (%)',
                                                'n best_non_self species specific'],
                                 set2taxon2value_list_simple_barplot=set2taxon2values,
                                 set2taxon2value_list_simple_barplot_counts=True)


    path1 = settings.BASE_DIR + '/assets/temp/interpro_tree2.svg'
    asset_path1 = '/temp/interpro_tree2.svg'
    tree1.render(path1, dpi=600, tree_style=style1)

    #pairwiseid_plots.density_plot([identity_values],["identity hits Eukaryota"])
    #pairwiseid_plots.basic_plot(identity_values, count_n_species, output_path="~/tata.svg")
    #pairwiseid_plots.basic_plot(count_n_species, output_path="~/tata2.svg")

    return render(request, 'chlamdb/blastnr_euk.html', my_locals(locals()))

def prot_length_barchart(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 ' \
               ' on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb


    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    from chlamdb.phylo_tree_display import phylo_tree_bar

    if db_driver == 'mysql':
        sql = 'select taxon_id, CHAR_LENGTH(translation) from orthology_detail;'
    if db_driver == 'sqlite':
        sql = 'select taxon_id, LENGTH(translation) from orthology_detail;'
    # blast_hits_taxonomy_overview

    CDS_length_data = server.adaptor.execute_and_fetchall(sql,)

    taxon_id2values = {}


    for row in CDS_length_data:
        if row[0] not in taxon_id2values:
            taxon_id2values[row[0]] = [row[1]]
        else:
            taxon_id2values[row[0]].append(row[1])

    taxon_id2CDS_length_counts = {}
    for taxon in taxon_id2values:
        taxon_data = taxon_id2values[taxon]

        len_50_99 = 0
        len_100_199 = 0
        len_200_299 = 0
        len_300_399 = 0
        len_400_499 = 0
        len_500_599 = 0
        len_more_600 = 0

        for CDS in taxon_data:
            if int(CDS) >49 and int(CDS) < 100:
                len_50_99+=1
            elif int(CDS) >99 and int(CDS) < 200:
                len_100_199+=1
            elif int(CDS) >199 and int(CDS) < 300:
                len_200_299+=1
            elif int(CDS) >299 and int(CDS) < 400:
                len_300_399+=1
            elif int(CDS) >399 and int(CDS) < 500:
                len_400_499+=1
            elif int(CDS) >499 and int(CDS) < 600:
                len_500_599+=1
            elif int(CDS) >599:
                len_more_600+=1
            else:
                pass
        taxon_id2CDS_length_counts[taxon] = [[len_50_99, len_100_199, len_200_299, len_300_399, len_400_499, len_500_599, len_more_600]]

    header_list = ['size']

    tree1, style1 = phylo_tree_bar.plot_tree_stacked_barplot(tree,
                                                    taxon_id2CDS_length_counts,
                                                    header_list,
                                                    biodb=biodb)

    path1 = settings.BASE_DIR + '/assets/temp/CDS_length.svg'
    asset_path1 = '/temp/CDS_length.svg'
    tree1.render(path1, dpi=600, tree_style=style1)
    return render(request, 'chlamdb/CDS_length.html', my_locals(locals()))


def blastnr_overview(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 ' \
               ' on t1.biodatabase_id=t2.biodatabase_id where name="%s";'



    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]


    header_list2 = ['effectiveT3',  'BPBAac', 'T3MM', 'T4SEpre_bpbAac', 'T4SEpre_psAac', 'chaperones', 'intesect']

    from chlamdb.phylo_tree_display import phylo_tree_bar
    from chlamdb.plots import hmm_heatmap
    from chlamdb.plots import pathway_heatmap
    from chlamdb.plots import module_heatmap

    set2taxon2value = {} #, column_names = hmm_heatmap.get_set_data(biodb, score_cutoff=20)

    sql_checkm_completeness = 'select taxon_id, completeness from custom_tables_checkm;'
    sql_checkm_n_duplicated = 'select taxon_id,n_total from custom_tables_checkm;'
    sql_n_contigs = 'select taxon_id,n_contigs from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)
    sql_n_no_CDS = 'select taxon_id,n_no_CDS from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)
    sql_n_no_BBH_chlamydiae = 'select taxon_id,n_no_BBH_chlamydiae from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    sql_n_no_BBH_chlamydiae = 'select taxon_id,n_no_BBH_chlamydiae from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    sql_GC = 'select taxon_id,GC from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    sql_genome_size = 'select taxon_id,genome_size/1000000 from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    taxon_id2completeness = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_checkm_completeness))
    taxon_id2duplicate = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_checkm_n_duplicated))
    taxon_id2n_no_CDS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_n_no_CDS))
    taxon_id2n_contigs = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_n_contigs))
    taxon_id2n_contigs_no_BBH_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_n_no_BBH_chlamydiae))
    taxon_id2GC = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_GC))
    taxon_id2genome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_genome_size))

    taxon_id2GC['127'] = 50
    # , 'T6SSi' 'rinke_et_al_2013', 'dupont_et_al_2012', 'eisen_et_al_2013'
    my_sets = ['duplicates', 'contigs','contigs_no_CDS', 'contigs_no_BBH', 'ADP/ATP carrier protein']

    #
    set2taxon2value['duplicates'] = taxon_id2duplicate
    set2taxon2value['contigs'] = taxon_id2n_contigs
    set2taxon2value['contigs_no_CDS'] = taxon_id2n_no_CDS
    set2taxon2value['contigs_no_BBH'] = taxon_id2n_contigs_no_BBH_chlamydiae


    set2taxon2value_list_simple_barplot = {}
    set2taxon2value_list_simple_barplot['gc'] = taxon_id2GC
    set2taxon2value_list_simple_barplot['complet.'] = taxon_id2completeness
    set2taxon2value_list_simple_barplot['size'] = taxon_id2genome_size

    header3 = ['gc', 'size', 'complet.' ]

    sql_ntt_transporters = 'select * from comparative_tables_interpro where id="IPR004667";'
    ntt_data = list(server.adaptor.execute_and_fetchall(sql_ntt_transporters)[0])
    sql_headers = 'show columns from  comparative_tables_interpro'
    ntt_taxid = [i[0] for i in server.adaptor.execute_and_fetchall(sql_headers)]

    taxon_id2count_ntt = dict(zip(ntt_taxid[1:], ntt_data[1:]))

    set2taxon2value['ADP/ATP carrier protein'] = taxon_id2count_ntt

    flagellum_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map02040"])
    peptidoglycan_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map00550"])
    purines_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map00230"])
    pyrim_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map00240"])

    set2taxon2value['Flagellum'] = flagellum_data[1]['Flagellar assembly']
    set2taxon2value['Peptidoglycan biosynthesis'] = peptidoglycan_data[1]['Peptidoglycan biosynthesis']
    set2taxon2value['Purine metabolism'] = purines_data[1]['Purine metabolism']
    set2taxon2value['Pyrimidine metabolism'] = pyrim_data[1]['Pyrimidine metabolism']

    module_list, code2taxon2count_modules = module_heatmap.module_list2profile_dico(biodb, ["M00001","M00004","M00009"])

    set2taxon2value['gylcolysis'] = code2taxon2count_modules['Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate [PATH:map01200 map00010]']
    set2taxon2value['TCA'] = code2taxon2count_modules['Citrate cycle (TCA cycle, Krebs cycle) [PATH:map01200 map00020]']
    set2taxon2value['PPP'] = code2taxon2count_modules['Pentose phosphate pathway (Pentose phosphate cycle) [PATH:map01200 map00030]']

    my_sets+= ["T3SS",
               "T4SS",
               'Flagellum',
                  "chemosensory",
                  'Peptidoglycan biosynthesis',
                  "Menaquinone",
                  "Menaquinone/futalosine pathway",
                  "NADH-quinone oxidoreductase",
                  "NQR",
                  "NDH-2",
                  "SDH",
                  "Cytochrome bc1 complex",
                  "Cbb3-type cytochrome c oxidase",
                  "Cytochrome bo(3) ubiquinol oxidase",
                  "Cytochrome bd-I ubiquinol oxidase",
                  "F-type ATPase",
                  "V-type ATPase",
                   ]

    my_sets+=['TCA',
              'gylcolysis',
              'PPP', "glycogen metabolism",'Purine metabolism','Pyrimidine metabolism']

    ######################
    # taxonomy project: locus int
    ######################
    #my_sets.append('ntt')
    try:
        # group by orthogroup ==> even if one prot has multiple homologs, count only one
        sql_complex = 'select C.taxon_id, count(*) as n from (select distinct A.category, A.orthogroup, B.taxon_id ' \
              ' from (select distinct t1.category,t2.orthogroup from custom_tables_annot_table t1 ' \
              ' inner join orthology_detail t2 on t1.locus_tag=t2.locus_tag ' \
              ' where category="%s") A inner join orthology_detail B ' \
              ' on A.orthogroup=B.orthogroup) C group by C.category,C.taxon_id;'
        # not group by orthogroup
        sql_simple = 'select C.taxon_id, count(*) as n from (select A.category, A.orthogroup, B.taxon_id ' \
              ' from (select distinct t1.category,t2.orthogroup from custom_tables_annot_table t1 ' \
              ' inner join orthology_detail t2 on t1.locus_tag=t2.locus_tag ' \
              ' where category="%s") A inner join orthology_detail B ' \
              ' on A.orthogroup=B.orthogroup) C group by C.category,C.taxon_id;'

        c_list = ["chemosensory",
                  "T4SS",
                  "T3SS",
                  "F-type ATPase",
                  "V-type ATPase",
                  "glycogen metabolism",
                  "Menaquinone",
                  "Menaquinone/futalosine pathway",
                  "NADH-quinone oxidoreductase",
                  "SDH",
                  "NQR",
                  "Cytochrome bc1 complex",
                  "Cbb3-type cytochrome c oxidase",
                  "Cytochrome bo(3) ubiquinol oxidase",
                  "Cytochrome bd-I ubiquinol oxidase"
                   ]
        for category in c_list:
            sql = sql_complex % ("2017_06_29b_motile_chlamydiae", biodb, category, biodb)
            #sql_ftype = sql_complex % (biodb, biodb, biodb)
            #sql_T3SS = sql_complex % (biodb, biodb, biodb)
            #sql_T4SS = sql_complex % (biodb, biodb, biodb)
            #sql_chemo = sql_complex % (biodb, biodb, biodb)
            #sql_futa = sql_complex % (biodb, biodb, biodb)
            #sql_mena = sql_complex % (biodb, biodb, biodb)

            #taxon_id2euo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_euo,))
            #taxon_id2vtype = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_vtype,))
            #taxon_id2ftype = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_ftype,))
            #taxon_id2T3SS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_T3SS,))
            #taxon_id2T4SS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_T4SS,))
            #taxon_id2chemo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_chemo,))
            #taxon_id2futa = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_futa,))
            taxon_id2counts = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            set2taxon2value[category] = taxon_id2counts
            #my_sets.append(category)
            #set2taxon2value['euo'] = taxon_id2euo
            #set2taxon2value['V-type ATPase'] = taxon_id2vtype
            #set2taxon2value['F-type ATPase'] = taxon_id2ftype
            #set2taxon2value['T3SS'] = taxon_id2T3SS
            #set2taxon2value['T4SS'] = taxon_id2T4SS
            #set2taxon2value['chemosensory'] = taxon_id2chemo
            #set2taxon2value['Menaquinone/futalosine pathway'] = taxon_id2futa
            #set2taxon2value['Menaquinone'] = taxon_id2mena

        #my_sets+=['chemosensory', 'V-type ATPase', 'F-type ATPase', 'T4SS', 'T3SS', 'Menaquinone/futalosine pathway',
        #          'Menaquinone']

        s_list = ["NDH-2"]#

        s_list = ["NDH-2",
                  "uhpC",
                  "CPAF",
                  "CopN",
                  "NUE",
                  "pkn5",
                  "mip",
                  "Lda3",
                  "Tarp",
                  "TepP",
                  "CT_847",
                  "CT_868",
                  "CT_867",
                  "CT_082",
                  "capN",
                  "CT_695",
                  "CT_620",
                  "CT_621",
                  "CT_622",
                  "CT_711",
                  "CT_694",
                  "CT_163",
                  "CT_365",
                  "CT_610",
                  "CT_157",
                  "CT_119",
                  "CT_115",
                  "CT_118",
                  "CT_228",
                  "CT_229",
                  "CT_813",
                  "CT_850",
                  "CPn0585",
                  "CPn0517",
                  "euo",
                  "OmcB",
                  "OmpA",
                  "OmpA1",
                  "pmp1",
                  "pmp2",
                  "Ctad1",
                  "FtsH",
                  "FtsI/Pbp3",
                  "FtsY",
                  "FtsK",
                  "FtsQ",
                  "FtsL",
                  "FtsW",
                  "MreB",
                  "RodZ",
                  "Ddl",
                  "DacB",
                  "DacC",
                  "DapA",
                  "DapB",
                  "DapL(AspC3)",
                  "DapF",
                  "MraY",
                  "MurJ (MviN)",
                  "Pbp2",
                  "AmiA",
                  "AmiB",
                  "NlpD",
                  "Pal",
                  "TolQ",
                  "TolR",
                  "TolA",
                  "TolA",
                  "TolB"]

        for locus in s_list:
            if locus == 'ntt':
                continue
            sql_locus = sql_simple % ("2017_06_29b_motile_chlamydiae",locus)
            taxid2locus_count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_locus,))
            set2taxon2value[locus] = taxid2locus_count
            my_sets.append(locus)
    except:
        pass




    #######################
    # metabolism aa, cofactors, nucleotides

    sql = 'select taxon_id,count(*) from ' \
          ' (select distinct t1.taxon_id,t1.ko_id,pathway_category from enzyme_locus2ko t1 ' \
          ' inner join enzyme_pathway2ko_v1 t2 on t1.ko_id=t2.ko_id ' \
          ' inner join enzyme_kegg_pathway t3 on t2.pathway_id=t3.pathway_id) A ' \
          ' where pathway_category="%s"  group by A.taxon_id,pathway_category;'

    #taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % (biodb,"1.4 Nucleotide metabolism",)))
    #set2taxon2value['Nucleotide metabolism'] = taxon_id2nucleo
    #my_sets.append('Nucleotide metabolism')
    taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.8 Metabolism of cofactors and vitamins",)))
    set2taxon2value['Cofactors and vitamins'] = taxon_id2nucleo
    my_sets.append('Cofactors and vitamins')
    taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.5 Amino acid metabolism",)))
    set2taxon2value['AA metabolism'] = taxon_id2nucleo
    my_sets.append('AA metabolism')
    #taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.10 Biosynthesis of other secondary metabolites",)))
    #set2taxon2value['1.10 Biosynthesis of other secondary metabolites'] = taxon_id2nucleo
    #my_sets.append('1.10 Biosynthesis of other secondary metabolites')
    #taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.1 Carbohydrate metabolism",)))
    #set2taxon2value['Carbohydrate metabolism'] = taxon_id2nucleo
    #my_sets.append('Carbohydrate metabolism')
    taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.3 Lipid metabolism",)))
    set2taxon2value['Lipid metabolism'] = taxon_id2nucleo
    my_sets.append('Lipid metabolism')


    sql = 'select taxon_id, count(*) from COG_seqfeature_id2best_COG_hit t1 ' \
          ' inner join COG_cog_names_2014 t2 on t1.hit_cog_id=t2.COG_id ' \
          ' inner join COG_cog_id2cog_category ta on t1.hit_cog_id=ta.COG_id' \
          ' inner join COG_code2category tb on ta.category_id=tb.category_id' \
          ' inner join annotation_seqfeature_id2locus tc on t1.seqfeature_id=tc.seqfeature_id' \
          ' where tb.code="X" group by taxon_id;'

    taxon_id2mobile_elements = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    set2taxon2value['modile_elements'] = taxon_id2mobile_elements
    my_sets.append('modile_elements')




    set2taxon2value_new = {}
    for set in my_sets:
        set2taxon2value_new[set] = {}
        for taxon in set2taxon2value[set]:

            try:
                value = set2taxon2value[set][taxon]
                #if value < 4:
                #    value = 0
                set2taxon2value_new[set][taxon] = value
            except:
                set2taxon2value_new[set][taxon] = 0

    sql = 'select taxon_id, n_no_hits, n_less_100_hits, n_100_hits from blastnr_count_n_blast order by n_no_hits;'

    # blast_hits_taxonomy_overview
    sql2 = 'select * from blastnr_BBH_taxo_hit_number_1;'
    sql3 = 'select * from blastnr_BBH_taxo_hit_number_2;'


    sql4 = 'select t2.taxon_id, sum(n_CDS) from genomes_info t1 ' \
           ' inner join bioentry t2 on t1.ACCESSION=t2.ACCESSION ' \
           ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id where t3.name="%s" group by taxon_id;' % (biodb)


    taxon_id2blast_count = server.adaptor.execute_and_fetchall(sql,)
    taxon_id2BBH_1_taxonomy = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    taxon_id2BBH_2_taxonomy = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))
    taxon_id2proteome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))




    taxon_id2values = {}
    for row in taxon_id2blast_count:
        taxo1 = list(taxon_id2BBH_1_taxonomy[row[0]])
        taxo2 = list(taxon_id2BBH_2_taxonomy[row[0]])
        taxo1.append(int(taxon_id2proteome_size[str(row[0])])-sum(taxo1))
        taxo2.append(int(taxon_id2proteome_size[str(row[0])])-sum(taxo2))
        taxo_filter = [taxo1[0], sum([taxo1[1], taxo1[2], taxo1[3], taxo1[4],taxo1[5]]), taxo1[6]]
        taxo_filter2 = [taxo2[0], sum([taxo2[1], taxo2[2], taxo2[3], taxo2[4],taxo2[5]]), taxo2[6]]
        taxon_id2values[row[0]] = [list(reversed(row[1:])), taxo1, taxo2]
    header_list = ['nr_n_hits', "nr_tax.1","nr_tax.2"]


    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;'

    sql2 = 'select taxon_id, count(*) as n, majority_phylum from blastnr_blastnr_majority_phylum ' \
           ' group by taxon_id,majority_phylum order by taxon_id,n DESC;'

    data = server.adaptor.execute_and_fetchall(sql,)
    data2 = server.adaptor.execute_and_fetchall(sql2,)


    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]
    taxon_match = []
    for row in data2:
        if row[0] not in taxon_match:
            taxon2most_freq_phylum[row[0]] += "/%s" % row[2]
            taxon_match.append(row[0])

    tree1, style1 = phylo_tree_bar.plot_tree_stacked_barplot(tree,
                                                    taxon_id2values,
                                                    header_list,
                                                    taxon2set2value_heatmap=set2taxon2value_new,
                                                    header_list2=my_sets,
                                                    biodb=biodb,
                                                    general_max=False,
                                                    taxon2label=taxon2most_freq_phylum,
                                                    set2taxon2value_list_simple_barplot=set2taxon2value_list_simple_barplot,
                                                    header_list3=header3) # taxon2most_freq_phylum

    # col = '#fc8d59' # col = '#91bfdb' '#99d594'
    path1 = settings.BASE_DIR + '/assets/temp/interpro_tree2.svg'
    asset_path1 = '/temp/interpro_tree2.svg'
    tree1.render(path1, dpi=600, tree_style=style1)
    return render(request, 'chlamdb/blastnr_overview.html', my_locals(locals()))



def blastnr_top_non_phylum(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    blastnr_form_class = make_blastnr_best_non_top_phylum_form(biodb)

    if request.method == 'POST':
        form = blastnr_form_class(request.POST)

        if form.is_valid():
            import pandas as pd

            accessions = form.cleaned_data['accessions']
            selection = form.cleaned_data['selection']

            taxon_filter = ','.join(accessions)

            if selection == 'all':
                sql = 'select t4.locus_tag,t5.product,t5.gene,t1.hit_number,t1.percent_identity,t3.kingdom,t3.class,' \
                      ' t3.order,t3.family,t3.species,t1.subject_accession,t1.subject_title from blastnr_blastnr_best_non_self_phylum t1' \
                      ' inner join blastnr_blastnr_taxonomy t3 on t1.subject_taxon_id=t3.taxon_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t4 on t1.seqfeature_id=t4.seqfeature_id ' \
                      ' inner join orthology_detail t5 on t4.locus_tag=t5.locus_tag  ' \
                      'where t1.superkingdom="Eukaryota"' \
                      ' and query_taxon_id in (%s)' % (taxon_filter)

            if selection == 'specific':
                sql = 'select t4.locus_tag,t5.product,t5.gene,t1.hit_number,t1.percent_identity,t3.kingdom,t3.class,' \
                      't3.order,t3.family,t3.species,t1.subject_accession,t1.subject_title from blastnr_blastnr_best_non_self_phylum t1 ' \
                      ' inner join custom_tables_seqfeature_id2n_species t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                      ' inner join blastnr_blastnr_taxonomy t3 on t1.subject_taxon_id=t3.taxon_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t4 on t1.seqfeature_id=t4.seqfeature_id ' \
                      ' inner join orthology_detail t5 on t4.locus_tag=t5.locus_tag  ' \
                      ' where t1.superkingdom="Eukaryota" and n_species=1  and t1.query_taxon_id in (%s);' % (taxon_filter)

            data = server.adaptor.execute_and_fetchall(sql,)
            envoi = True

    else:  
        form = blastnr_form_class()
    return render(request, 'chlamdb/blastnr_locus_list.html', my_locals(locals()))

def blastnr_barchart(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    blastnr_form_class = make_blastnr_form(biodb)

    if request.method == 'POST':
        form = blastnr_form_class(request.POST)

        if form.is_valid():
            import pandas as pd

            target_accessions = form.cleaned_data['accession']
            rank = form.cleaned_data['rank']
            counttype = form.cleaned_data['type']
            top_n = form.cleaned_data['top_number']

            biodb_id_sql = 'select biodatabase_id from biodatabase where name="%s"'
            biodb_id = server.adaptor.execute_and_fetchall(biodb_id_sql,)[0][0]
            sql_accession = 'select bioentry_id,description from bioentry where biodatabase_id=%s and bioentry_id in (%s)' % (biodb_id,
                                                                                                                              '"'+'","'.join(target_accessions)+'"')
            taxon2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_accession,))

            data_all_accessions = []
            for accession in target_accessions:

                if counttype == 'Majority':
                    sql = 'select seqfeature_id,%s, count(*) as n from blastnr_blastnr A ' \
                          ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id where query_bioentry_id=%s and hit_number<=%s' \
                          ' group by seqfeature_id,%s order by seqfeature_id, n DESC' % (rank,
                                                                                         accession,
                                                                                         top_n,
                                                                                         rank)

                    data = server.adaptor.execute_and_fetchall(sql,)
                    category2count = {}
                    query_locus_list = []
                    for n, i in enumerate(data):
                        # KEEP ONY the first match (highest count ordered with mysql)
                        if i[0] not in query_locus_list:
                            if i[1] not in category2count:
                                category2count[i[1]] = 1
                            else:
                                category2count[i[1]] += 1
                            query_locus_list.append(i[0])
                            try:
                                if data[n+1][0] == data[n][0]:
                                    if data[n+1][2] == data[n][2]:
                                        print("Idem!:",data[n+1], data[n])
                            except:
                                pass

                    data = zip(category2count.keys(), category2count.values())

                elif counttype == 'BBH':
                    sql = ' select %s, count(*) as n from (select * from blastnr_blastnr ' \
                          ' where query_bioentry_id=%s and hit_number=1) A inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id ' \
                          ' group by %s;' % (rank,
                                             accession,
                                             rank)
                    data = server.adaptor.execute_and_fetchall(sql,)
                else:
                    raise 'invalide type'

                for i in data:
                    data_all_accessions.append((accession,) + i)

            taxon_map = 'var taxon2description = {'
            for i in taxon2description:
                taxon_map += '"%s":"%s",' % (i, taxon2description[i])
            taxon_map = taxon_map[0:-1] + '};'

            taxon2category2count = {}
            all_categories = []
            for line in data_all_accessions:
                if line[0] not in taxon2category2count:
                    taxon2category2count[line[0]] = {}
                    taxon2category2count[line[0]][line[1]] = line[2]
                else:
                    taxon2category2count[line[0]][line[1]] = line[2]
                if line[1] not in all_categories:
                    all_categories.append(line[1])
            labels_template = '[\n' \
                              '%s\n' \
                              ']\n'

            serie_template = '[%s\n' \
                             ']\n'

            one_serie_template = '{\n' \
                                 'label: "%s",\n' \
                                 'values: [%s]\n' \
                                 '},\n'

            # count number of hits for each category and order based on the number of hits
            category2count = {}
            for taxon in taxon2category2count:
                for category in all_categories:
                    try:
                        if category not in category2count:
                            category2count[category] = int(taxon2category2count[taxon][category])
                        else:
                            category2count[category] += int(taxon2category2count[taxon][category])
                    except:
                        pass
            for key in category2count:
                print("%s\t%s" % (key, category2count[key]))
            data = pd.DataFrame({'category': list(category2count.keys()),
                    'count': list(category2count.values()) })
            data_sort = data.sort_values("count", ascending=0)

            all_series_templates = []
            for taxon in taxon2category2count:
                one_category_list = []
                for category in data_sort['category']:
                    try:
                        one_category_list.append(taxon2category2count[taxon][category])
                    except:
                        one_category_list.append(0)
                one_category_list = [str(i) for i in one_category_list]
                all_series_templates.append(one_serie_template % (taxon, ','.join(one_category_list)))

            series = serie_template % ''.join(all_series_templates)
            cat_list = [str(i) for i in data_sort['category']]
            labels = labels_template % ('"'+'","'.join(cat_list) + '"')


            circos_url = '?h=' + ('&h=').join(target_accessions) + '&t=%s&n=%s' % (counttype, top_n)


            envoi = True
    else:  
        form = blastnr_form_class()
    return render(request, 'chlamdb/blastnr_best_barplot.html', my_locals(locals()))



def effector_pred(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    from ete3 import Tree
    from chlamdb.phylo_tree_display import ete_motifs

    from chlamdb.phylo_tree_display import phylo_tree_bar
    from chlamdb.plots import hmm_heatmap

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
    t1 = Tree(tree)
    taxon_id_list = [i.name for i in t1.iter_leaves()]

    set2taxon2value, column_names = hmm_heatmap.get_set_data(biodb, score_cutoff=30)

    sql = 'select '

    my_sets = ['T3SS', 'T6SSi', 'T4SS', 'flagellum', 'rinke_et_al_2013', 'dupont_et_al_2012', 'eisen_et_al_2013']

    set2taxon2value_new = {}
    for set in my_sets:
        set2taxon2value_new[set] = {}
        for taxon in set2taxon2value[set]:

            try:
                value = set2taxon2value[set][taxon]
                if value < 4:
                    value = 0
                set2taxon2value_new[set][taxon] = value
            except:
                set2taxon2value_new[set][taxon] = 0

    sql = 'select taxon_id,n_CDS from biodatabase t1 inner join bioentry t2 on ' \
          ' t1.biodatabase_id=t2.biodatabase_id inner join genomes_info t3 ' \
          ' on t3.ACCESSION=t2.accession where t1.name="%s" ' \
          ' and t3.description not like "%%%%plasmid%%%%";' % (biodb)

    taxon_id2genome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # effector prediction data
    # effective T3
    sql = 'select taxon_id, count(*) from effectors_predicted_effectiveT3 where score>0 group by taxon_id;'
    taxon2values_effectiveT3 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # BPBAac
    sql = 'select taxon_id, count(*) from effectors_predicted_BPBAac where SVM_value>0 group by taxon_id;'
    taxon2values_BPBAac = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # T3MM
    sql = 'select taxon_id, count(*) from effectors_predicted_T3MM where value>0 group by taxon_id;'
    taxon2values_T3MM = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # T4SEpre_bpbAac
    sql = 'select taxon_id, count(*) from effectors_predicted_T4SEpre_bpbAac where SVM_value>0 group by taxon_id;'
    taxon2values_T4SEpre_bpbAac = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # T4SEpre_psAac
    sql = 'select taxon_id, count(*) from effectors_predicted_T4SEpre_psAac where SVM_value>0 group by taxon_id;'
    taxon2values_T4SEpre_psAac = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # chapeones
    sql = 'select A.taxon_id, count(*) from (select * from effectors_predicted_chaperones group by seqfeature_id,taxon_id) A group by A.taxon_id;'
    taxon2values_chaperones = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # ELD
    sql = 'select A.taxon_id, count(*) from (select * from effectors_predicted_ELD ' \
          ' where score >9 group by seqfeature_id) A group by A.taxon_id;'

    taxon2values_eld = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # pfam refseq tanonomy
    sql = 'select taxon_id, count(*) as n from (select distinct t5.taxon_id,t1.pfam_id,t5.locus_tag from interpro_interpro_signature2pfam_id t1 ' \
         ' inner join interpro_interpro t2 on t1.signature_id=t2.signature_id ' \
         ' inner join pfam.pfam2superkingdom_frequency_31 t3 on t1.pfam_id=t3.pfam_id  inner join interpro_signature t4 on t1.signature_id=t4.signature_id ' \
         ' inner join annotation_seqfeature_id2locus t5 on t2.seqfeature_id=t5.seqfeature_id where bacteria_freq<=0.02 and eukaryota_count>5) BBB group by taxon_id;'

    taxon2pfam_refseq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # pfam refseq tanonomy: number of unique pfam domain
    sql = 'select taxon_id, count(*) as n from (select distinct t5.taxon_id,t1.pfam_id from interpro_interpro_signature2pfam_id t1 ' \
         ' inner join interpro_interpro t2 on t1.signature_id=t2.signature_id ' \
         ' inner join pfam.pfam2superkingdom_frequency_31 t3 on t1.pfam_id=t3.pfam_id  inner join interpro_signature t4 on t1.signature_id=t4.signature_id ' \
         ' inner join annotation_seqfeature_id2locus t5 on t2.seqfeature_id=t5.seqfeature_id where bacteria_freq<=0.02 and eukaryota_freq>=0.1) BBB group by taxon_id;'

    taxon2pfam_refseq_uniq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select A.taxon_id, count(*) from (select t1.taxon_id, t1.seqfeature_id from effectors_predicted_effectiveT3 t1 ' \
          ' inner join effectors_predicted_BPBAac t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' inner join effectors_predicted_T3MM t3 on t1.seqfeature_id=t3.seqfeature_id ' \
          ' group by t1.taxon_id, t1.seqfeature_id) A group by A.taxon_id;'

    taxon2values_mix = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    # counts
    taxon2values = {}
    taxon2values2 = {}
    for taxon in taxon_id_list:

        if taxon not in taxon2values_mix:
            taxon2values_mix[taxon] = [0, 0]
        else:
            taxon2values_mix[taxon] = [int(taxon2values_mix[taxon]),
                                       round((float(taxon2values_mix[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_effectiveT3:
            taxon2values_effectiveT3[taxon] = [0, 0]
        else:
            taxon2values_effectiveT3[taxon] = [int(taxon2values_effectiveT3[taxon]),
                                               round((float(taxon2values_effectiveT3[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_BPBAac:
            taxon2values_BPBAac[taxon] = [0, 0]
        else:

            taxon2values_BPBAac[taxon] = [int(taxon2values_BPBAac[taxon]),
                                          round((float(taxon2values_BPBAac[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_T3MM:
            taxon2values_T3MM[taxon] = [0, 0]
        else:
            taxon2values_T3MM[taxon] = [int(taxon2values_T3MM[taxon]),
                                        round((float(taxon2values_T3MM[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_T4SEpre_bpbAac:
            taxon2values_T4SEpre_bpbAac[taxon] = [0, 0]
        else:
            taxon2values_T4SEpre_bpbAac[taxon] = [int(taxon2values_T4SEpre_bpbAac[taxon]),
                                                  round((float(taxon2values_T4SEpre_bpbAac[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_T4SEpre_psAac:
            taxon2values_T4SEpre_psAac[taxon] = [0, 0]
        else:
            taxon2values_T4SEpre_psAac[taxon] = [int(taxon2values_T4SEpre_psAac[taxon]),
                                                 round((float(taxon2values_T4SEpre_psAac[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_eld:
            taxon2values_eld[taxon] = [0, 0]
        else:
            taxon2values_eld[taxon] = [int(taxon2values_eld[taxon]),
                                       round((float(taxon2values_eld[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2values_chaperones:
            taxon2values_chaperones[taxon] = [0, 0]
        else:
            taxon2values_chaperones[taxon] = [int(taxon2values_chaperones[taxon]),
                                       round((float(taxon2values_chaperones[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]
        if taxon not in taxon2pfam_refseq:
            taxon2pfam_refseq[taxon] = [0, 0]
        else:
            taxon2pfam_refseq[taxon] = [int(taxon2pfam_refseq[taxon]),
                                       round((float(taxon2pfam_refseq[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]

        if taxon not in taxon2pfam_refseq_uniq:
            taxon2pfam_refseq_uniq[taxon] = [0, 0]
        else:
            taxon2pfam_refseq_uniq[taxon] = [int(taxon2pfam_refseq_uniq[taxon]),
                                       round((float(taxon2pfam_refseq_uniq[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]

        if taxon not in taxon2pfam_refseq:
            taxon2pfam_refseq[taxon] = [0, 0]
        else:
            taxon2pfam_refseq[taxon] = [int(taxon2pfam_refseq[taxon]),
                                       round((float(taxon2pfam_refseq[taxon])/float(taxon_id2genome_size[taxon]))*100,2)]


        taxon2values[taxon] = [taxon2values_effectiveT3[taxon][0],
                               taxon2values_BPBAac[taxon][0],
                               taxon2values_T3MM[taxon][0],
                               taxon2values_T4SEpre_bpbAac[taxon][0],
                               taxon2values_T4SEpre_psAac[taxon][0],
                               taxon2values_chaperones[taxon][0],
                               taxon2values_eld[taxon][0],
                               taxon2values_mix[taxon][0],
                               taxon2pfam_refseq[taxon][0],
                               taxon2pfam_refseq_uniq[taxon][0],
                               ]

        taxon2values2[taxon] = [taxon2values_effectiveT3[taxon][1],
                               taxon2values_BPBAac[taxon][1],
                               taxon2values_T3MM[taxon][1],
                               taxon2values_T4SEpre_bpbAac[taxon][1],
                               taxon2values_T4SEpre_psAac[taxon][1],
                               taxon2values_chaperones[taxon][1],
                               taxon2values_eld[taxon][1],
                               taxon2values_mix[taxon][1],
                               taxon2pfam_refseq[taxon][1],
                               taxon2pfam_refseq_uniq[taxon][1],
                               ]

    header_list2 = ['effectiveT3', 'BPBAac', 'T3MM', 'T4SEpre_bpbAac', 'T4SEpre_psAac', 'chaperones','ELD', 'intesect' , 'refseq_pfam', 'refseq_pfam_uniq']


    sql = 'SELECT orthogroup, count(*) as n FROM (select  orthogroup,taxon_id from orthology_detail ' \
          ' group by orthogroup,taxon_id) A  GROUP BY orthogroup' % biodb
    group2n_organisms = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select C.orthogroup from (select t1.taxon_id, t1.seqfeature_id from effectors_predicted_effectiveT3 t1 ' \
          ' inner join effectors_predicted_BPBAac t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' inner join effectors_predicted_T3MM t3 on t1.seqfeature_id=t3.seqfeature_id ' \
          ' group by t1.taxon_id, t1.seqfeature_id) A inner join custom_tables_locus2seqfeature_id B ' \
          ' on A.seqfeature_id=B.seqfeature_id inner join orthology_detail C on ' \
          ' B.locus_tag=C.locus_tag group by C.orthogroup;'
          

    identity_list = []
    genome_count_list = []
    group_list =[i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    for group in group_list:
        if int(group2n_organisms[group]) > 50:
            print (group2n_organisms[group], group)
        genome_count_list.append(group2n_organisms[group])

    #import pairwiseid_plots
    #pairwiseid_plots.basic_plot(genome_count_list)

    general_max = 0 
    for taxon in taxon2values2:
        m = max([float(i) for i in taxon2values2[taxon]])
        if m > general_max:
            general_max=m
    general_max=False

    tree1, style1 = phylo_tree_bar.plot_tree_barplot(tree,
                                                    taxon2values,
                                                    header_list2,
                                                    taxon2set2value_heatmap=set2taxon2value_new,
                                                    header_list2=my_sets,
                                                    biodb=biodb,
                                                    general_max=general_max)

    path = settings.BASE_DIR + '/assets/temp/interpro_tree.svg'
    asset_path = '/temp/interpro_tree.svg'
    tree1.render(path, dpi=600, tree_style=style1)

    tree2, style2 = phylo_tree_bar.plot_tree_barplot(tree,
                                                    taxon2values2,
                                                    header_list2,
                                                    taxon2set2value_heatmap=set2taxon2value_new,
                                                    header_list2=my_sets,
                                                    biodb=biodb,
                                                    general_max=False)
    path2 = settings.BASE_DIR + '/assets/temp/interpro_tree2.svg'
    asset_path2 = '/temp/interpro_tree2.svg'
    tree2.render(path2, dpi=600, tree_style=style2)

    all=True

    return render(request, 'chlamdb/effector_pred.html', my_locals(locals()))



def interpro_taxonomy(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    interpro_form_class = make_interpro_taxonomy(biodb)

    if request.method == 'POST':
        form = interpro_form_class(request.POST)

        if form.is_valid():
            from ete3 import Tree
            from chlamdb.phylo_tree_display import ete_motifs

            target_taxons = form.cleaned_data['target_taxons']
            kingdom = form.cleaned_data['kingdom']
            percentage_cutoff = form.cleaned_data['percentage_cutoff']

            if target_taxons[0] == 'all':
                from chlamdb.phylo_tree_display import phylo_tree_bar
                from chlamdb.plots import hmm_heatmap

                # counts eukaryotic domains
                sql = 'select taxon_id, eukaryote_count from interpro_taxonomy_summary_50 ;'

                taxon2values = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                sql = 'select taxon_id, eukaryote_count from interpro_taxonomy_summary_90 ;'

                taxon2values_90 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                sql = 'select taxon_id, eukaryote_count from interpro_taxonomy_summary_98 ;'

                taxon2values_98 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb



                tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]







                header_list = ['euk_50', 'euk_90', 'euk_98', 'ELD']

                sql = 'select A.taxon_id, count(*) from (select * from effectors_predicted_ELD ' \
                      ' where score >9 group by seqfeature_id) A group by A.taxon_id;'

                taxon2values_eld = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # secretion systems
                for taxon in taxon2values:
                    taxon2values[taxon] = [taxon2values[taxon],
                                           taxon2values_90[taxon],
                                           taxon2values_98[taxon],
                                           taxon2values_eld[taxon]]
                set2taxon2value, column_names = hmm_heatmap.get_set_data(biodb, score_cutoff=30)

                sql = 'select '

                my_sets = ['T3SS', 'T6SSi', 'T4SS', 'flagellum', 'rinke_et_al_2013', 'dupont_et_al_2012', 'eisen_et_al_2013']

                set2taxon2value_new = {}
                for set in my_sets:
                    set2taxon2value_new[set] = {}
                    for taxon in set2taxon2value[set]:

                        try:
                            value = set2taxon2value[set][taxon]
                            if value < 4:
                                value = 0
                            set2taxon2value_new[set][taxon] = value
                        except:
                            set2taxon2value_new[set][taxon] = 0

                sql = 'select taxon_id,n_CDS from biodatabase t1 inner join bioentry t2 on ' \
                      ' t1.biodatabase_id=t2.biodatabase_id inner join genomes_info t3 ' \
                      ' on t3.ACCESSION=t2.accession where t1.name="%s" ' \
                      ' and t3.description not like "%%%%plasmid%%%%";' % (biodb)

                taxon_id2genome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # effector prediction data
                # effective T3
                sql = 'select taxon_id, count(*) from effectors_predicted_effectiveT3 where score>0 group by taxon_id;'
                taxon2values_effectiveT3 = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # BPBAac
                sql = 'select taxon_id, count(*) from effectors_predicted_BPBAac where SVM_value>0 group by taxon_id;'
                taxon2values_BPBAac = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # T3MM
                sql = 'select taxon_id, count(*) from effectors_predicted_T3MM where value>0 group by taxon_id;'
                taxon2values_T3MM = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # T4SEpre_bpbAac
                sql = 'select taxon_id, count(*) from effectors_predicted_T4SEpre_bpbAac where SVM_value>0 group by taxon_id;'
                taxon2values_T4SEpre_bpbAac = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # T4SEpre_psAac
                sql = 'select taxon_id, count(*) from effectors_predicted_T4SEpre_psAac where SVM_value>0 group by taxon_id;'
                taxon2values_T4SEpre_psAac = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # chapeones
                sql = 'select A.taxon_id, count(*) from (select * from effectors_predicted_chaperones group by seqfeature_id,taxon_id) A group by A.taxon_id;'
                taxon2values_chaperones = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                # mex 3 algo
                sql = 'select A.taxon_id, count(*) from (select t1.taxon_id, t1.seqfeature_id from effectors_predicted_effectiveT3 t1 ' \
                      ' inner join effectors_predicted_BPBAac t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                      ' inner join effectors_predicted_T3MM t3 on t1.seqfeature_id=t3.seqfeature_id ' \
                      ' group by t1.taxon_id, t1.seqfeature_id) A group by A.taxon_id;'
                taxon2values_mix = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                for taxon in taxon2values:
                    if taxon not in taxon2values_mix:
                        taxon2values_mix[taxon] = 0
                    else:
                        taxon2values_mix[taxon] = round((float(taxon2values_mix[taxon])/float(taxon_id2genome_size[taxon]))*100,2)
                    if taxon not in taxon2values_effectiveT3:
                        taxon2values_effectiveT3[taxon] = 0
                    else:
                        taxon2values_effectiveT3[taxon] = round((float(taxon2values_effectiveT3[taxon])/float(taxon_id2genome_size[taxon]))*100,2)
                    if taxon not in taxon2values_BPBAac:
                        taxon2values_BPBAac[taxon] = 0
                    else:

                        taxon2values_BPBAac[taxon] = round((float(taxon2values_BPBAac[taxon])/float(taxon_id2genome_size[taxon]))*100,2)
                    if taxon not in taxon2values_T3MM:
                        taxon2values_T3MM[taxon] = 0
                    else:
                        taxon2values_T3MM[taxon] = round((float(taxon2values_T3MM[taxon])/float(taxon_id2genome_size[taxon]))*100,2)
                    if taxon not in taxon2values_T4SEpre_bpbAac:
                        taxon2values_T4SEpre_bpbAac[taxon] = 0
                    else:
                        taxon2values_T4SEpre_bpbAac[taxon] = round((float(taxon2values_T4SEpre_bpbAac[taxon])/float(taxon_id2genome_size[taxon]))*100,2)
                    if taxon not in taxon2values_T4SEpre_psAac:
                        taxon2values_T4SEpre_psAac[taxon] = 0
                    else:
                        taxon2values_T4SEpre_psAac[taxon] = round((float(taxon2values_T4SEpre_psAac[taxon])/float(taxon_id2genome_size[taxon]))*100,2)

                taxon2values2 = {}
                for taxon in taxon2values:
                    taxon2values2[taxon] = [taxon2values_effectiveT3[taxon],
                                           taxon2values_BPBAac[taxon],
                                           taxon2values_T3MM[taxon],
                                           taxon2values_T4SEpre_bpbAac[taxon],
                                           taxon2values_T4SEpre_psAac[taxon],
                                           taxon2values_chaperones[taxon],
                                           taxon2values_mix[taxon]]

                header_list2 = ['effectiveT3', 'BPBAac', 'T3MM', 'T4SEpre_bpbAac', 'T4SEpre_psAac', 'chaperones', 'intesect']


                sql = 'SELECT orthogroup, count(*) as n FROM (select  orthogroup,taxon_id from orthology_detail ' \
                      ' group by orthogroup,taxon_id) A  GROUP BY orthogroup'
                group2n_organisms = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                sql = 'select C.orthogroup from (select t1.taxon_id, t1.seqfeature_id from effectors_predicted_effectiveT3 t1 ' \
                      ' inner join effectors_predicted_BPBAac t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                      ' inner join effectors_predicted_T3MM t3 on t1.seqfeature_id=t3.seqfeature_id ' \
                      ' group by t1.taxon_id, t1.seqfeature_id) A inner join custom_tables_locus2seqfeature_id B ' \
                      ' on A.seqfeature_id=B.seqfeature_id inner join orthology_detail C on ' \
                      ' B.locus_tag=C.locus_tag group by C.orthogroup;'
                      
                identity_list = []
                genome_count_list = []
                group_list =[i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
                for group in group_list:
                    if int(group2n_organisms[group]) > 50:
                        print (group2n_organisms[group], group)
                    genome_count_list.append(group2n_organisms[group])
                from chlamdb.plots import pairwiseid_plots
                pairwiseid_plots.basic_plot(genome_count_list)


                #general_max = 0
                #for taxon in taxon2values:
                #    m = max([float(i) for i in taxon2values[taxon]])
                #    if m > general_max:
                #        general_max=m

                tree1, style1 = phylo_tree_bar.plot_tree_barplot(tree,
                                                                taxon2values,
                                                                header_list,
                                                                taxon2set2value_heatmap=set2taxon2value_new,
                                                                header_list2=my_sets,
                                                                biodb=biodb,
                                                                 general_max=False)

                path = settings.BASE_DIR + '/assets/temp/interpro_tree.svg'
                asset_path = '/temp/interpro_tree.svg'
                tree1.render(path, dpi=600, tree_style=style1)

                #for taxon in taxon2values2:
                #    m = max([float(i) for i in taxon2values2[taxon]])
                #    if m > general_max:
                #        general_max=m

                tree2, style2 = phylo_tree_bar.plot_tree_barplot(tree,
                                                                taxon2values2,
                                                                header_list2,
                                                                taxon2set2value_heatmap=set2taxon2value_new,
                                                                header_list2=my_sets,
                                                                biodb=biodb,
                                                                general_max=False)
                path2 = settings.BASE_DIR + '/assets/temp/interpro_tree2.svg'
                asset_path2 = '/temp/interpro_tree2.svg'
                tree2.render(path2, dpi=600, tree_style=style2)

                all=True

            else:
                filter = ','.join(target_taxons)


                #p_eukaryote | p_archae | p_virus

                sql = 'select * from (select distinct interpro_accession from interpro where ' \
                      ' interpro_accession!="0" and taxon_id in (%s))A inner join interpro_entry B on A.interpro_accession=B.name ' \
                      ' inner join interpro_interpro_taxonomy_v_60 C on B.interpro_id=C.interpro_id where %s>=%s;' % (filter,
                                                                                                                      kingdom,
                                                                                                                      percentage_cutoff)

                data = server.adaptor.execute_and_fetchall(sql,)
                interpro2description = {}
                for row in data:
                    interpro2description[row[0]] = row[3]

                interpro_accession_list = [i[0] for i in data]
                filter2 = '"' + '","'.join(interpro_accession_list) + '"'
                sql2 = 'select taxon_id, interpro_accession, count(*) from ' \
                       ' (select taxon_id,locus_tag,interpro_accession from interpro ' \
                       ' where interpro_accession in (%s) group by taxon_id,locus_tag,interpro_accession) A ' \
                       ' group by A.taxon_id,A.interpro_accession' % (filter2)

                data_counts = server.adaptor.execute_and_fetchall(sql2,)

                sql3 = 'select B.* from ' \
                       ' (select locus_tag from interpro ' \
                       ' where interpro_accession in (%s) and taxon_id in (%s) group by locus_tag) A ' \
                       ' inner join orthology_detail B on A.locus_tag=B.locus_tag ' % (filter2,
                                                                                       filter)

                data_locus = server.adaptor.execute_and_fetchall(sql3,)

                data_locus = server.adaptor.execute_and_fetchall(sql3,)

                product2count = {}
                for locus in data_locus:
                    if locus[9] not in product2count:
                        product2count[locus[9]] = 1
                    else:
                        product2count[locus[9]] += 1




                locus_list = [i[3] for i in data_locus]
                fasta_url = '?l=' + '&l='.join(locus_list)

                sql_tree = 'select tree from reference_phylogeny as t1 ' \
                           ' inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";'

                tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]


                acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
                filter = '"' + '","'.join(acc_list) + '"'
                sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
                taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]

                t1 = Tree(tree)
                R = t1.get_midpoint_outgroup()
                t1.set_outgroup(R)

                # number of groups with identified signature domains
                sql = 'select name,n from (select AA.interpro_id, count(*) as n from ' \
                      ' (select distinct A.interpro_id, B.orthogroup_id from ' \
                      ' (select distinct seqfeature_id,t2.interpro_id from interpro_interpro t1 ' \
                      ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
                      ' inner join interpro_interpro_taxonomy_v_60 t3 on t2.interpro_id=t3.interpro_id where %s>=%s) A ' \
                      ' inner join orthology_seqfeature_id2orthogroup B on A.seqfeature_id=B.seqfeature_id ' \
                      ' inner join orthology_orthogroup C on B.orthogroup_id=C.orthogroup_id) AA ' \
                      ' group by AA.interpro_id) BB inner join interpro_entry CC on BB.interpro_id=CC.interpro_id;' % (kingdom,
                                                                                                                       percentage_cutoff)

                interpro_accession2n_groups = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                interpro2taxon2count = {}
                interpro_list = ['TOTAL']
                taxon2total = {}
                interpro_id_description = []
                for row in data_counts:
                    try:
                        interpro_des = "%s: %s (%s groups)" % (row[1], interpro2description[row[1]], interpro_accession2n_groups[row[1]])
                    except:
                        interpro_des = "%s: %s (? groups)" % (row[1], interpro2description[row[1]])
                    if interpro_des not in interpro_list:
                        interpro_list.append(interpro_des)
                        interpro_id_description.append(interpro_des.split(':'))
                    if row[0] not in taxon2total:
                        taxon2total[row[0]] = 0
                    if interpro_des not in interpro2taxon2count:
                        interpro2taxon2count[interpro_des] = {}
                        interpro2taxon2count[interpro_des][str(row[0])] = int(row[2])

                    else:
                            interpro2taxon2count[interpro_des][str(row[0])] = int(row[2])
                    taxon2total[row[0]] += int(row[2])
                interpro2taxon2count['TOTAL'] = {}
                for taxon in taxon2total:
                    interpro2taxon2count['TOTAL'][str(taxon)] = taxon2total[taxon]

                tree2, style = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            interpro_list,
                                                            interpro2taxon2count,
                                                            show_labels=True,
                                                            column_scale=True,
                                                            tree=t1,
                                                            rotate=True)

                style.rotation = 90
                path = settings.BASE_DIR + '/assets/temp/interpro_tree.svg'
                asset_path = '/temp/interpro_tree.svg'
                tree2.render(path, dpi=600, tree_style=style)

                #path2 = settings.BASE_DIR + '/assets/temp/COG_tree_%s_complete.svg' % module_name
                #asset_path2 = '/assets/temp/KEGG_tree_%s_complete.svg' % modulextract_interproe_name

                #tree2.render(path2, dpi=800, h=600)
                not_all = True



    else:  
        form = interpro_form_class()
    return render(request, 'chlamdb/interpro_taxonomy.html', my_locals(locals()))

def format_gene(gene):
    if pd.isna(gene):
        return "-"
    else:
        return gene


def blastswissprot(request, locus_tag):
    biodb = settings.BIODB
    print ('blast swissprot %s -- %s' % (biodb, locus_tag))


    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'GET': 

        server, db = manipulate_biosqldb.load_db(biodb)

        columns = 'hit_number,subject_accession,subject_kingdom,subject_scientific_name,subject_taxid,' \
                  ' subject_title,evalue,bit_score,percent_identity,gaps,query_cov,genes,annot_score'
        sql = 'select A.*, B.phylum, B.order, B.family from (select %s from custom_tables_locus2seqfeature_id t1 ' \
              ' inner join blastnr_blast_swissprot t2 on t1.seqfeature_id=t2.seqfeature_id' \
              ' where locus_tag="%s") A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id' % (columns,locus_tag)
        blast_data = server.adaptor.execute_and_fetchall(sql,)

        if len(blast_data) > 0:
            valid_id = True
            #'<a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s">%s<a> ' % (taxon, name)




        return render(request, 'chlamdb/blastswiss.html', my_locals(locals()))


    return render(request, 'chlamdb/blastswiss.html', my_locals(locals()))


def blastnr(request, locus_tag):
    biodb = settings.BIODB

    print ('blastnr hits %s -- %s' % (biodb, locus_tag))

    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'GET': 

        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession, organism from orthology_detail where locus_tag="%s"' % (biodb, locus_tag)
        data = server.adaptor.execute_and_fetchall(sql,)[0]
        accession = data[0]
        organism = data[1]

        server, db = manipulate_biosqldb.load_db(biodb)
        columns = 'hit_number, subject_accession, subject_kingdom, subject_scientific_name, ' \
                  ' subject_taxid, subject_title, evalue, bit_score, percent_identity, gaps, length'
        sql = 'select hit_number,subject_accession,superkingdom,subject_scientific_name,subject_taxid,subject_title,evalue,bit_score,percent_identity,gaps, length,' \
              ' B.phylum, B.order, B.family from (select %s from custom_tables_locus2seqfeature_id t1 inner join blastnr_blastnr t2' \
              ' on t1.seqfeature_id=t2.seqfeature_id where locus_tag="%s" order by hit_number) A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id;' % (columns, locus_tag)

        blast_data = list(server.adaptor.execute_and_fetchall(sql))

        blast_data = list(server.adaptor.execute_and_fetchall(sql))

        if len(blast_data) > 0:
            valid_id = True
            blast_query_locus = blast_data[0][1]
            blast_query_protein_id = blast_data[0][2]


        return render(request, 'chlamdb/blastnr.html', my_locals(locals()))


    return render(request, 'chlamdb/blastnr.html', my_locals(locals()))



def homology(request):
    biodb = settings.BIODB
    from chlamdb.biosqldb import shell_command

    server, db = manipulate_biosqldb.load_db(biodb)

    contact_form_class = make_contact_form(server, biodb)

    if request.method == 'POST': 

        form = contact_form_class(request.POST)
        #form2 = ContactForm(request.POST)
        if form.is_valid():  
            accession = request.POST['accession']

            envoi = True
    else:  
        form = make_contact_form(server, biodb)
    return render(request, 'chlamdb/homology.html', my_locals(locals()))



def orthogroup_identity(request, orthogroup, group=False):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    print ('orthogroup identity -- %s -- %s' % (biodb, orthogroup))
    #if request.method == 'POST':
    import numpy
    import pandas as pd
    sql = 'SELECT locus_a, locus_b, identity FROM orthology_identity where orthogroup="%s";' % (orthogroup)
    
    
    sql2 = 'select locus_a, locus_b from orthology_identity where orthogroup in ("%s")' % (orthogroup)

    #try:
        # numpy.array()
    data = [list(i) for i in server.adaptor.execute_and_fetchall(sql,)]
    locus_list_of_list = server.adaptor.execute_and_fetchall(sql2,)
    # get nr list of locus
    locus_list = list(set([locus for locus_list in locus_list_of_list for locus in locus_list]))

    # create dictionnary
    locus2locus2identity = {}
    for row in data:
        if row[0] not in locus2locus2identity:
            locus2locus2identity[row[0]] = {}
            locus2locus2identity[row[0]][row[1]] = row[2]
        else:
            locus2locus2identity[row[0]][row[1]] = row[2]
        if row[1] not in locus2locus2identity:
            locus2locus2identity[row[1]] = {}
            locus2locus2identity[row[1]][row[0]] = row[2]
        else:
            locus2locus2identity[row[1]][row[0]] = row[2]
    data_matrix = []
    for locus_a in locus_list:
        tmp_lst = [locus_a]
        for locus_b in locus_list:
            tmp_lst.append(locus2locus2identity[locus_a][locus_b])
        data_matrix.append(tmp_lst)
    data = numpy.array(data_matrix)
    homologs = True
    #except:
    #    homologs = False
    #    return render(request, 'chlamdb/orthogroup_identity.html', my_locals(locals()))
    locus_list_filter = '"' + '","'.join(locus_list) + '"'

    sql2 = 'select locus_tag, organism from orthology_detail where locus_tag in (%s)' % (locus_list_filter)

    locus2organism = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    columns = locus_list
    rows = [i + " (%s)" % locus2organism[i] for i in columns]
    frame = pd.DataFrame(data[0:,1:], index=rows, columns=columns)
    frame = frame.astype(float)
    frame = frame/100

    for i in range(0, len(frame)):

        frame.ix[i, i] = None

    path = settings.BASE_DIR + '/assets/temp/%s.json' % orthogroup

    with open(path, 'w') as f:
        f.write(frame.to_json(orient="split"))

    return render(request, 'chlamdb/orthogroup_identity.html', my_locals(locals()))


def ortho_id_plot(request, group):
    return render(request, 'chlamdb/orthogroup_identity_plot.html', my_locals(locals()))





def plot_neighborhood(request, target_locus, region_size=23000):
    biodb = settings.BIODB

    print ("plot region %s -- %s " % (biodb, target_locus))

    server, db = manipulate_biosqldb.load_db(biodb)

    task = plot_neighborhood_task.delay(biodb, target_locus, region_size)

    task_id = task.id
    
    return render(request, 'chlamdb/plot_region_and_profile.html', my_locals(locals()))

def plot_region_generic(biodb, orthogroup, taxon_list, region_size):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)


    filter = '"' + '","'.join(taxon_list) + '"'
    sql3 = 'select locus_tag from orthology_detail where orthogroup = "%s" and taxon_id in (%s)' % (orthogroup, filter)

    locus_tag_target_genomes = [i[0] for i in server.adaptor.execute_and_fetchall(sql3, )]

    home_dir = os.path.dirname(__file__)
    temp_location = os.path.join(home_dir, "../assets")
    temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
    name = os.path.basename(temp_file.name)
    name_png = name.split('.')[0] + '.png'

    locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, 
                                                                                    db, 
                                                                                    biodb, 
                                                                                    locus_tag_target_genomes,
                                                                                    temp_file.name, 
                                                                                    int(region_size),
                                                                                    cache)


    return name, name_png, locus_tags, orthogroup_list

def plot_region_direct(request, orthogroup):
    biodb = settings.BIODB
    target_taxons = [str(i) for i in request.GET.getlist('t')]

    name, name_png, locus_tags, orthogroup_list = plot_region_generic(biodb, 
                                                                      orthogroup, 
                                                                      target_taxons, 
                                                                      18000)

    return render(request, 'chlamdb/plot_region_simple.html', my_locals(locals()))


def plot_region(request):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    plot_region_form_class = make_plot_form(biodb)
    if request.method == 'POST': 

        form = plot_region_form_class(request.POST)

        if form.is_valid():  
            valid_id = True
            import re

            accession = extract_alphanumeric(form.cleaned_data['accession'])

            m = re.match("CT([0-9]+)", accession)
            if m:
                accession = "CT_%s" % (m.group(1))

            server, db = manipulate_biosqldb.load_db(biodb)

            region_size = form.cleaned_data['region_size']

            genomes = form.cleaned_data['genomes']

            all_homologs = form.cleaned_data['all_homologs']
            if all_homologs == 'yes':
                closest_only=False
            else:
                closest_only=True

            sql2 = 'select orthogroup, locus_tag, protein_id, start, stop, strand, organism ' \
                   ' from orthology_detail where locus_tag="%s" or protein_id="%s"' % (accession, accession)

            data = server.adaptor.execute_and_fetchall(sql2, )[0]

            if not data:
                    valid_id = False
            if valid_id:
                locus_tag = data[1]
                if closest_only:
                    sql = 'select B.locus_tag,taxon_2,identity from (select * from custom_tables_locus2seqfeature_id t1 ' \
                          ' inner join comparative_tables_identity_closest_homolog2 t2 on t1.seqfeature_id=t2.locus_1 ' \
                          ' where t1.locus_tag="%s") A ' \
                          ' inner join custom_tables_locus2seqfeature_id B on A.locus_2=B.seqfeature_id ' \
                          ' order by identity DESC;' % (locus_tag)

                    target_locus_data = server.adaptor.execute_and_fetchall(sql,)
                    # keep only the hit with the highest identity
                    genome_hits = []
                    locus_tag_target_genomes=[locus_tag]
                    for one_homolog in target_locus_data:
                        if str(one_homolog[1]) in genomes:
                            if one_homolog[1] not in genome_hits:
                                genome_hits.append(one_homolog[1])
                                locus_tag_target_genomes.append(one_homolog[0])

                else:
                    orthogroup = data[0]

                    select = 'and (taxon_id = %s' % genomes[0]
                    if len(genomes) >1:
                        for i in range(0, len(genomes)-1):
                            select+= ' or taxon_id = %s' % genomes[i]
                        select+= ' or taxon_id = %s)' % genomes[-1]
                    else:
                        select+= ')'
                    sql3 = 'select locus_tag from orthology_detail where orthogroup = "%s" %s' % (orthogroup, select)
                    locus_tag_target_genomes = [i[0] for i in server.adaptor.execute_and_fetchall(sql3, )]

                if plot_region:
                    home_dir = os.path.dirname(__file__)
                    temp_location = os.path.join(home_dir, "../assets")
                    temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
                    name = os.path.basename(temp_file.name)
                    name_png = name.split('.')[0] + '.png'

                    locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, 
                                                                                                    db, 
                                                                                                    biodb, 
                                                                                                    locus_tag_target_genomes,
                                                                                                    temp_file.name, 
                                                                                                    int(region_size),
                                                                                                    cache)

                    columns = 'orthogroup, locus_tag, protein_id, start, stop, ' \
                              'strand, gene, orthogroup_size, n_genomes, TM, SP, product, organism, translation'

                    sql_locus = 'locus_tag="%s"' % locus_tags[0]
                    for locus in range(1, len(locus_tags)):
                        sql_locus += ' or locus_tag="%s"' % locus_tags[locus]

                    sql = 'select %s from orthology_detail where %s' % (columns,
                                                                        sql_locus)

                    raw_data = server.adaptor.execute_and_fetchall(sql,)

                    n = 1
                    search_result = []
                    for one_hit in raw_data:
                        search_result.append((n,) + one_hit)
                        n+=1

            envoi = True

    else:  
        form = plot_region_form_class()

    return render(request, 'chlamdb/plot_region.html', my_locals(locals()))

'''

def plot_region(request, biodb):
    plot_form_class = make_plot_form(biodb)

    print "cache", cache
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST': 


        form = plot_form_class(request.POST)

        if form.is_valid():  






            envoi = True

    else:  
        form = plot_form_class()

    return render(request, 'chlamdb/plot_region.html', my_locals(locals()))
'''



def comparative_extract():
    '''
    get features present in genomes X Y Z and not in A B C

    Se baser sur la table d'orthologie pour identifier les genes

    SQL
    select orthogroup from orthology_chlam where `taxon X` > 1 and `taxon B` = 0

    :return:
    '''



def orthogroups(request):
    if request.method == 'POST':
        form = BiodatabaseForm(request.POST)
        if form.is_valid():
            biodb = form.cleaned_data['biodatabase']
            groups = "orthogroup_size_distrib_%s.svg" % biodb
            envoi = True
    else:
        form = BiodatabaseForm()

    return render(request, 'chlamdb/orthogroups.html', my_locals(locals()))


def sitemap(request):
    from io import StringIO
    map = '''<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
<url><loc>http://chlamdb.ch/</loc></url><url><loc>http://chlamdb.ch/about</loc></url>
</urlset>'''

    strio = StringIO()
    strio.write(map)
    response = HttpResponse(content_type='text/xml')
    #response.content_type = "text/xml"
    #response['Content-Disposition'] = 'attachment; filename="sitemap.xml"'
    response.write(strio.getvalue())
    return response


def get_orthogroup_fasta(request, orthogroup, seqtype):
    biodb = settings.BIODB
    from io import StringIO
    from Bio.Seq import Seq

    server, db = manipulate_biosqldb.load_db(biodb)

    if seqtype == 'aa':
        sql = 'select locus_tag, organism, translation from orthology_detail where orthogroup="%s"' % (orthogroup)

        data = server.adaptor.execute_and_fetchall(sql,)
        fasta = []
        for i in data:
            biorecord = SeqRecord(Seq(i[2]),
                                  id=i[0],
                                  name=i[0],
                                  description=i[1])
            fasta.append(biorecord)
    else:
        sql = 'select accession, locus_tag, start, stop, strand from orthology_detail where orthogroup="%s"' % (orthogroup)

        locus2start_stop = server.adaptor.execute_and_fetchall(sql,)
        fasta = []

        for i in locus2start_stop:
            leng = i[3]-i[2]+1
            strand = int(i[4])
            seq = manipulate_biosqldb.location2sequence(server, i[0], biodb, int(i[2]), leng)
            if strand == -1:

                seq_obj = Seq(seq)
                seq = seq_obj.reverse_complement()
                biorecord = SeqRecord(seq,
                                      id=i[1],
                                      name=i[1],
                                      description=i[0])
                fasta.append(biorecord)

    strio = StringIO()
    SeqIO.write(fasta, strio, 'fasta') #+='>%s %s\n%s\n' % (i[1], i[0], seq)

    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="%s_fasta.fa"' % orthogroup
    response.write(strio.getvalue())
    return response

def get_newick_tree(request, orthogroup, refseq_BBH_phylogeny):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if refseq_BBH_phylogeny == "False":
        sql_tree = 'select phylogeny from phylogenies where orthogroup="%s"' % (orthogroup)
    else:
        sql_tree = 'select phylogeny from phylogenies_BBH where orthogroup="%s";' % (orthogroup)
    try:
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
    except IndexError:
        tree = 'No tree for orthogroup: %s' % orthogroup

    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="%s_tree.nwk"' % orthogroup
    response.write(tree)
    return response


def module2fasta():

    sql = 'select t4.locus_tag, organism, gene,product, translation from ' \
          ' (select module_name,module_id from kegg_module_v1 where module_name="M00009") t1 ' \
          ' inner join module2ko_v1 as t2 on t1.module_id=t2.module_id inner join locus2ko_%s as t3 ' \
          ' on t2.ko_id=t3.ko_id inner join orthology_detail as t4 ' \
          ' on t3.locus_tag=t4.locus_tag where t4.taxon_id in (64);'


def pfam2fasta(request, pfam_id):
    biodb = settings.BIODB
    from io import StringIO
    from Bio.Seq import Seq

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select protein_id,product,translation from interpro_interpro t1 ' \
          ' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
          ' inner join annotation_seqfeature_id2CDS_annotation t3 on t1.seqfeature_id=t3.seqfeature_id ' \
          ' where signature_accession="%s";' % (pfam_id)

    data = server.adaptor.execute_and_fetchall(sql,)
    fasta = []
    for i in data:
        biorecord = SeqRecord(Seq(i[2]),
                              id=i[0],
                              name=i[0],
                              description=i[1])
        fasta.append(biorecord)

    strio = StringIO()
    SeqIO.write(fasta, strio, 'fasta')

    response = HttpResponse(content_type='text/plain')
    name = 'attachment; filename="pfam_%s.fa"' % pfam_id

    response['Content-Disposition'] = name
    response.write(strio.getvalue())
    return response


def ko2fasta(request, ko_id, include_orthologs=False):
    biodb = settings.BIODB
    from io import StringIO
    from Bio.Seq import Seq

    server, db = manipulate_biosqldb.load_db(biodb)
    if not include_orthologs:
        sql = 'select B.locus_tag, organism, B.product, translation from (' \
              ' select t1.ko_id, locus_tag,ko_description from enzyme_locus2ko as t1 inner join enzyme_module2ko_v1 as t2 ' \
              ' on t1.ko_id=t2.ko_id where t1.ko_id="%s") A inner join orthology_detail as B ' \
              ' on A.locus_tag=B.locus_tag;' % (ko_id)
    else:
        sql = 'select distinct orthogroup from (' \
              ' select t1.ko_id, locus_tag, ko_description from enzyme_locus2ko as t1 inner join enzyme_module2ko_v1 as t2 ' \
              ' on t1.ko_id=t2.ko_id where t1.ko_id="%s") A inner join orthology_detail as B ' \
              ' on A.locus_tag=B.locus_tag;' % (ko_id)

        orthogroup_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

        filter = '"' + '","'.join(orthogroup_list) + '"'
        sql = 'select locus_tag, organism, product, translation from orthology_detail where orthogroup in (%s)' % (filter)

    data = server.adaptor.execute_and_fetchall(sql,)
    fasta = []
    for i in data:
        biorecord = SeqRecord(Seq(i[3]),
                              id=i[0],
                              name=i[0],
                              description=i[1])
        fasta.append(biorecord)

        #fasta+='>%s %s (%s)\n%s\n' % (i[0], i[2], i[1], i[3])
    strio = StringIO()
    SeqIO.write(fasta, strio, 'fasta')

    response = HttpResponse(content_type='text/plain')
    if not include_orthologs:
        name = 'attachment; filename="%s.fa"' % ko_id
    else:
        name = 'attachment; filename="%s_ortho.fa"' % ko_id
    response['Content-Disposition'] = name
    response.write(strio.getvalue())
    return response

def fasta(request):
    biodb = settings.BIODB
    '''
    get fasta from locus_list

    :param request:
    :param biodb:
    :return: fasta file
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    locus_list = [str(i) for i in request.GET.getlist('l')]

    filter = '"'+'","'.join(locus_list)+'"'

    sql = 'select locus_tag, organism, translation from orthology_detail where locus_tag in (%s)' % (filter)

    data = server.adaptor.execute_and_fetchall(sql,)
    fasta = ''
    for seq_data in data:
        fasta+='>%s %s\n' % (seq_data[0], seq_data[1])
        for i in range(0, len(seq_data[2]), 60):
            fasta+=(seq_data[2][i:i + 60] + "\n")

    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="fasta.fa"'
    response.write(fasta)
    return response

def get_record(request, accession, seqtype):
    '''
    :param request:
    :param type: fna, faa, gbk, tab, ffn
    :return:
    '''
    from io import StringIO
    import re
    from Bio.Alphabet import IUPAC

    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    record = db.lookup(accession=accession)
    try:
        evalstring = eval(str(record.seq))
        if isinstance(evalstring, bytes):
            record.seq = Seq(str(eval(str(record.seq)),  'utf-8'), IUPAC.extended_dna)
    except NameError:
        pass
    strio = StringIO()
    if seqtype=='fna':
        SeqIO.write(record, strio, 'fasta')
        response = HttpResponse(content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="%s.fna"' % (accession)
    elif seqtype=='ffn':
        fasta_list = []
        for seq_feature in record.features:
            if seq_feature.type == "CDS":
                biorecord = SeqRecord(Seq(str(seq_feature.extract(record.seq)), IUPAC.extended_dna),
                                      id=seq_feature.qualifiers['locus_tag'][0],
                                      name=seq_feature.qualifiers['locus_tag'][0],
                                      description=record.description)
                fasta_list.append(biorecord)
        SeqIO.write(fasta_list, strio, 'fasta')
        response = HttpResponse(content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="%s.ffn"' % (accession)
    elif seqtype=='faa':
        fasta_list = []
        for seq_feature in record.features:
            if seq_feature.type == "CDS":
                biorecord = SeqRecord(Seq(re.sub('\*$', '', str(seq_feature.extract(record.seq).translate())), IUPAC.protein),
                                      id=seq_feature.qualifiers['locus_tag'][0],
                                      name=seq_feature.qualifiers['locus_tag'][0],
                                      description=record.description)
                fasta_list.append(biorecord)
        SeqIO.write(fasta_list, strio, 'fasta')
        response = HttpResponse(content_type='text/plain')
        response['Content-Disposition'] = 'attachment; filename="%s.fna"' % (accession)


    if seqtype == 'gbk':
        SeqIO.write(record, strio, 'genbank')

        response = HttpResponse(content_type='text/plain')

        response['Content-Disposition'] = 'attachment; filename="%s.gbk"' % (accession)
    else:
        pass
    response.write(strio.getvalue())
    return response


def get_fasta(request):
    biodb = settings.BIODB
    from django.http import FileResponse
    from chlamdb.biosqldb import shell_command
    import tarfile
    '''
    get fasta from a corresponding to extract_orthogroup_request

    :param request:
    :param biodb:
    :return: fasta file
    '''

    from chlamdb.biosqldb import biosql_own_sql_tables
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.GET.getlist('ref')[0] == 'False' or request.GET.getlist('ref')[0] == 'F':
        reference = False
    else:
        reference = str(request.GET.getlist('ref')[0])
    include = [str(i) for i in request.GET.getlist('i')]
    exclude = [str(i) for i in request.GET.getlist('e')]

    if exclude[0] == '':
        exclude = []
    if request.GET.getlist('a')[0] == 'F' or request.GET.getlist('a')[0] == 'False':
        accessions = False
    else:
        accessions = True
    freq_missing = float(request.GET.getlist('f')[0])
    if request.GET.getlist('s')[0] == 'F' or request.GET.getlist('s')[0] == 'False':
        single_copy = False
    else:
        single_copy = True

    freq_missing = freq_missing-0.001

    if not accessions:
        # get sub matrix and complete matrix
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                      "orthology",
                                                                      "orthogroup",
                                                                      include,
                                                                      exclude,
                                                                      ratio=freq_missing,
                                                                      single_copy=single_copy,
                                                                      accessions=accessions,
                                                                      cache=cache)
    else:
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                      "orthology",
                                                                      "id",
                                                                      include,
                                                                      exclude,
                                                                      ratio=freq_missing,
                                                                      single_copy=single_copy,
                                                                      accessions=accessions,
                                                                      cache=cache)
    match_groups = mat.index.tolist()

    merged_reference_fasta = ''
    for n, group in enumerate(match_groups):
        if not accessions:
            if reference:
                sql = 'select locus_tag, organism, translation from orthology_detail where taxon_id=%s and orthogroup="%s"' % (reference,
                                                                                                                               group)
            else:
                taxon_filter = '"'+'","'.join(include)+'"'
                sql = 'select taxon_id, organism, translation from orthology_detail where taxon_id in (%s) and orthogroup="%s"' % (taxon_filter,
                                                                                                                                   group)
        else:
            if reference:
                sql = 'select locus_tag, organism, translation from orthology_detail where accession="%s" and orthogroup="%s"' % (reference,
                                                                                                                                  group)
            else:
                taxon_filter = '"'+'","'.join(include)+'"'
                sql = 'select locus_tag, organism, translation from orthology_detail where accession in (%s) and orthogroup="%s"' % (taxon_filter,
                                                                                                                                     group)

        data = server.adaptor.execute_and_fetchall(sql,)
        fasta = ''
        for seq_data in data:
            if not reference:
                fasta+='>%s %s\n' % (seq_data[0], seq_data[1])
                for i in range(0, len(seq_data[2]), 60):
                    fasta+=(seq_data[2][i:i + 60] + "\n")
            else:
                merged_reference_fasta+='>%s %s\n' % (seq_data[0], seq_data[1])
                for i in range(0, len(seq_data[2]), 60):
                    merged_reference_fasta+=(seq_data[2][i:i + 60] + "\n")
        if not reference:          
            with open('/tmp/%s.fasta' % group, 'w') as f:
                f.write(fasta)
    
    response = HttpResponse(content_type='application/x-gzip')
    response['Content-Disposition'] = 'attachment; filename=download.tar.gz'
    tarred = tarfile.open(fileobj=response, mode='w:gz')
    
    if reference:
        with open('/tmp/reference.fasta', 'w') as f:
            f.write(merged_reference_fasta)
        tarred.add('/tmp/reference.fasta', arcname='reference.fasta')
    else:
        for group in match_groups:
            tarred.add('/tmp/%s.fasta' % group, arcname='%s.fasta' % group)
    tarred.close()           

    return response

    #return FileResponse(open('/tmp/groups_fasta.tar.gz', 'rb')) #HttpResponse(request, fasta, content_type='text/plain; charset=utf8')


def download_COG(request, accession=False): 

    import os
    import pandas
    import gzip
    from io import BytesIO
    from io import StringIO
    
    biodb = settings.BIODB
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn

    sql = f'select t3.accession,t3.description,locus_tag,start,stop,strand,gene,protein_id,product,t6.COG_name,t6.description ' \
          f' from annotation_seqfeature_id2locus t1 ' \
          f' inner join annotation_seqfeature_id2CDS_annotation t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          f' inner join bioentry t3 on t1.bioentry_id=t3.bioentry_id ' \
          f' inner join biodatabase t4 on t3.biodatabase_id=t4.biodatabase_id ' \
          f' inner join COG_seqfeature_id2best_COG_hit t5 on t1.seqfeature_id=t5.seqfeature_id ' \
          f' inner join COG_cog_names_2014 t6 on t5.hit_cog_id=t6.cog_id where t4.name="{biodb}" order by t3.accession'
    
    df = pandas.read_sql(sql, conn)
    zbuf = BytesIO()
    
    with gzip.GzipFile(fileobj=zbuf, compresslevel=6, mode='wb') as f:
        f.write(df.to_csv(sep="\t").encode('utf-8'))
    
    zbuf.seek(0)

    response = StreamingHttpResponse(zbuf, content_type='application/x-gzip')
    response['Content-Disposition'] = 'attachment; filename="COG_annotation.tsv.gz"'

    return response


def download_all_COG(request):
    import os
    import pandas
    import gzip
    from io import BytesIO
    from io import StringIO

    biodb = settings.BIODB
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn

    sql = f'select t3.accession,t3.description,locus_tag,start,stop,strand,gene,protein_id,product,t6.COG_name,t6.description ' \
          f' from annotation.seqfeature_id2locus_{biodb} t1 ' \
          f' inner join annotation.seqfeature_id2CDS_annotation_{biodb} t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          f' inner join biosqldb.bioentry t3 on t1.bioentry_id=t3.bioentry_id ' \
          f' inner join biosqldb.biodatabase t4 on t3.biodatabase_id=t4.biodatabase_id ' \
          f' inner join COG.seqfeature_id2best_COG_hit_{biodb} t5 on t1.seqfeature_id=t5.seqfeature_id ' \
          f' inner join COG.cog_names_2014 t6 on t5.hit_cog_id=t6.cog_id where t4.name="{biodb}" order by t3.accession'
    
    df = pandas.read_sql(sql, conn)
    zbuf = BytesIO()
    
    with gzip.GzipFile(fileobj=zbuf, compresslevel=6, mode='wb') as f:
        f.write(df.to_csv(sep="\t").encode('utf-8'))
    
    zbuf.seek(0)

    response = StreamingHttpResponse(zbuf, content_type='application/x-gzip')
    response['Content-Disposition'] = 'attachment; filename="COG_annotation.tsv.gz"'

    return response






def get_fasta_all(request):
    biodb = settings.BIODB

    '''
    get fasta from a corresponding to extract_orthogroup_request

    :param request:
    :param biodb:
    :return: fasta file
    '''

    from chlamdb.biosqldb import biosql_own_sql_tables
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.GET.getlist('ref')[0] == 'False' or request.GET.getlist('ref')[0] == 'F':
        reference = False
    else:
        reference = str(request.GET.getlist('ref')[0])
    include = [str(i) for i in request.GET.getlist('i')]
    exclude = [str(i) for i in request.GET.getlist('e')]

    if exclude[0] == '':
        exclude = []
    if request.GET.getlist('a')[0] == 'F' or request.GET.getlist('a')[0] == 'False':
        accessions = False
    else:
        accessions = True
    freq_missing = float(request.GET.getlist('f')[0])
    if request.GET.getlist('s')[0] == 'F' or request.GET.getlist('s')[0] == 'False':
        single_copy = False
    else:
        single_copy = True

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

    filter = '"'+'","'.join(match_groups)+'"'
    if not accessions:
        if reference:
            sql = 'select locus_tag, organism, translation from orthology_detail where taxon_id=%s and orthogroup in (%s)' % (reference,
                                                                                                                              filter)
        else:
            taxon_filter = '"'+'","'.join(include)+'"'
            sql = 'select locus_tag, organism, translation from orthology_detail where taxon_id in (%s) and orthogroup in (%s)' % (taxon_filter,
                                                                                                                                   filter)
    else:
        if reference:
            sql = 'select locus_tag, organism, translation from orthology_detail where accession="%s" and orthogroup in (%s)' % (reference,
                                                                                                                                 filter)
        else:
            taxon_filter = '"'+'","'.join(include)+'"'
            sql = 'select locus_tag, organism, translation from orthology_detail where accession in (%s) and orthogroup in (%s)' % (taxon_filter,
                                                                                                                                    filter)

    data = server.adaptor.execute_and_fetchall(sql,)
    fasta = ''
    for seq_data in data:
        fasta+='>%s %s\n' % (seq_data[0], seq_data[1])
        for i in range(0, len(seq_data[2]), 60):
            fasta+=(seq_data[2][i:i + 60] + "\n")
    response = HttpResponse(content_type='text/plain')
    response['Content-Disposition'] = 'attachment; filename="fasta.fa"'
    response.write(fasta)
    return response #HttpResponse(request, fasta, content_type='text/plain; charset=utf8')


def get_pfam_hit_list(request,
                      pfam_domain,
                      superkingdom=False,
                      phylum=False,
                      order=False,
                      family=False,
                      genus=False):

    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if superkingdom == '-':
        superkingdom = ''
    if phylum == '-':
        phylum = ''
    if order == '-':
        order = ''
    if family == '-':
        family = ''
    if genus == '-':
        genus = ''

    if superkingdom is False:
        sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,t2.assembly_accession,t2.organism_name,t1.protein_id,t1.evalue_full,t1.score_full from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
              ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
              ' where t3.hmm_accession like "%s%%%%" ' % (pfam_domain)
    elif phylum is False:
        sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,t2.assembly_accession,t2.organism_name,t1.protein_id,t1.evalue_full,t1.score_full from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
              ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
              ' where t3.hmm_accession like "%s%%%%" ' \
              ' and t4.superkingdom="%s" '  % (pfam_domain, superkingdom)
    elif order is False:
        sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,t2.assembly_accession,t2.organism_name,t1.protein_id,t1.evalue_full,t1.score_full from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
              ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
              ' where t3.hmm_accession like "%s%%%%" ' \
              ' and t4.superkingdom="%s" ' \
              ' and t4.phylum="%s" ' % (pfam_domain, superkingdom, phylum)
    elif family is False:
        sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,t2.assembly_accession,t2.organism_name,t1.protein_id,t1.evalue_full,t1.score_full from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
              ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
              ' where t3.hmm_accession like "%s%%%%" ' \
              ' and t4.superkingdom="%s" ' \
              ' and t4.phylum="%s" ' \
              ' and t4.order="%s" '  % (pfam_domain, superkingdom, phylum, order)
    elif genus is False:
        sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,t2.assembly_accession,t2.organism_name,t1.protein_id,t1.evalue_full,t1.score_full from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
              ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
              ' where t3.hmm_accession like "%s%%%%" ' \
              ' and t4.superkingdom="%s" ' \
              ' and t4.phylum="%s" ' \
              ' and t4.order="%s" ' \
              ' and t4.family="%s" ' % (pfam_domain, superkingdom, phylum, order, family)
    else:
        sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,t2.assembly_accession,t2.organism_name,t1.protein_id,t1.evalue_full,t1.score_full from pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
              ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
              ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
              ' where t3.hmm_accession like "%s%%%%" ' \
              ' and t4.superkingdom="%s" ' \
              ' and t4.phylum="%s" ' \
              ' and t4.order="%s" ' \
              ' and t4.family="%s" ' \
              ' and t4.genus="%s";' % (pfam_domain, superkingdom, phylum, order, family, genus)

    taxon_data = server.adaptor.execute_and_fetchall(sql,)

    return render(request, 'chlamdb/pfam_taxon_detail.html', my_locals(locals()))

def get_pfam_taxon_table(request, pfam_domain):

    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select t4.superkingdom,t4.phylum,t4.order,t4.family,t4.genus,count(*) as n from ' \
          ' pfam.refseq_ref_repres_genomes_domains_pfam_31 t1 ' \
          ' inner join pfam.refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
          ' inner join pfam.pfam_summary_version_31 t3 on t1.pfam_id=t3.hmm_id ' \
          ' inner join blastnr_blastnr_taxonomy t4 on t2.species_taxid=t4.taxon_id  ' \
          ' where t3.hmm_accession like "%s%%%%" ' \
          ' group by t4.superkingdom,t4.phylum,t4.family,t4.genus order by superkingdom,n DESC' % (pfam_domain)

    raw_data = server.adaptor.execute_and_fetchall(sql,)
    taxon_data = []
    for row in raw_data:
        row = list(row)
        if row[0] == '':
            row[0] = '-'
        if row[1] == '':
            row[1] = '-'
        if row[2] == '':
            row[2] = '-'
        if row[3] == '':
            row[3] = '-'
        if row[4] == '':
            row[4] = '-'
        taxon_data.append(row)

    return render(request, 'chlamdb/pfam_taxon_table.html', my_locals(locals()))


def pfam_profile(request, pfam_domain, rank):
    biodb = settings.BIODB
    from ete3 import Tree,TreeStyle
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.phylo_tree_display import pfam_phylogenetic_profile
    server, db = manipulate_biosqldb.load_db(biodb)
    path = settings.BASE_DIR + '/assets/temp/tree.svg'
    asset_path = '/temp/tree.svg'
    sql = 'select hmm_accession from pfam.pfam_summary_version_31 where hmm_accession like "%s%%%%"; ' % pfam_domain
    pf_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

    tree, style = pfam_phylogenetic_profile.plot_phylum_counts(pf_id,rank)

    tree.render(path, tree_style=style, dpi=800)
    return render(request, 'chlamdb/pfam_profile.html', my_locals(locals()))

def eggnog_profile(request, eggnog_id, rank):
    biodb = settings.BIODB

    from ete3 import Tree,TreeStyle
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import eggnog_data
    server, db = manipulate_biosqldb.load_db(biodb)
    path = settings.BASE_DIR + '/assets/temp/tree.svg'
    asset_path = '/temp/tree.svg'


    tree, style = eggnog_data.plot_phylum_counts(eggnog_id, rank=rank,colapse_low_species_counts=0)

    tree.render(path, tree_style=style) # dpi=800,
    return render(request, 'chlamdb/eggnog_profile.html', my_locals(locals()))

def annotation_overview(request):
    biodb = settings.BIODB
    from ete3 import Tree,TreeStyle
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.phylo_tree_display import phylo_tree_bar
    from chlamdb.plots import pairwiseid_plots

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select * from genomes_info'

    genomes_data = server.adaptor.execute_and_fetchall(sql,)

    sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
               ' where t2.name="%s";' % biodb
    server, db = manipulate_biosqldb.load_db(biodb)
    tree = Tree(server.adaptor.execute_and_fetchall(sql_tree,)[0][0])
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()

    sql_COG = 'select A.taxon_id,count(*) from (select t2.taxon_id,t1.locus_tag,COG_id ' \
              ' from COG_locus_tag2gi_hit t1 ' \
              ' inner join bioentry t2 on t1.accession=t2.accession ' \
              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' where t3.name="%s" group by t2.taxon_id,t1.locus_tag,COG_id) A group by A.taxon_id;' % (biodb)

    taxon_id2CDS_with_COG = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_COG,))

    sql_PRIAM = 'select A.taxon_id,count(*) from (select t2.taxon_id,t1.locus_tag,ec_id ' \
              ' from enzyme_locus2ec t1 ' \
              ' inner join bioentry t2 on t1.accession=t2.accession ' \
              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
              ' where t3.name="%s" group by t2.taxon_id,t1.locus_tag,ec_id) A group by A.taxon_id;' % (biodb)

    taxon_id2CDS_with_EC = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_PRIAM,))

    sql_KO = 'select taxon_id, count(*) from (select distinct taxon_id,locus_tag ' \
             ' from enzyme_locus2ko)A group by taxon_id;'

    taxon_id2CDS_with_KO = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_KO,))

    sql = 'select taxon_id,sum(n_CDS) from genomes_info t1 ' \
          ' inner join bioentry t2 on t1.accession=t2.accession ' \
          ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
          ' where t3.name="%s" group by taxon_id;' % (biodb)
    taxon_id2n_CDS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql_swiss = 'select query_taxon_id, count(*) from (select * from blastnr_blast_swissprot ' \
                ' group by query_taxon_id,seqfeature_id) A group by A.query_taxon_id;'
                
    taxon_id2n_swiss_hits = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_swiss,))

    sql_Pfam = 'select taxon_id, count(*) from (select distinct taxon_id,locus_tag from interpro ' \
               ' where analysis="Pfam") A group by taxon_id;' 

    taxon_id2n_Pfam = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_Pfam,))

    sql_string_mapping = 'select taxon_id, count(*) from custom_tables_uniprot_id2seqfeature_id t0 ' \
                         ' inner join custom_tables_uniprot_db_xref t1 on t0.uniprot_id=t1.uniprot_id ' \
                         ' inner join custom_tables_db_xref t2 on t1.db_xref_id=t2.db_xref_id ' \
                         ' inner join custom_tables_locus2seqfeature_id t3 on t0.seqfeature_id=t3.seqfeature_id ' \
                         ' where db_xref_name="string" group by taxon_id;'

    sql_uniprot_mapping = 'select taxon_id, count(*) from custom_tables_uniprot_id2seqfeature_id t0 ' \
                         ' inner join custom_tables_locus2seqfeature_id t3 on t0.seqfeature_id=t3.seqfeature_id ' \
                         ' group by taxon_id;'


    taxon_id2n_string = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_string_mapping,))
    taxon_id2n_uniprot_map = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_uniprot_mapping,))

    taxon_id2barplot_data = {}

    for taxon in taxon_id2n_CDS:
        if taxon not in taxon_id2barplot_data:
            taxon_id2barplot_data[taxon] = []

        taxon_id2barplot_data[taxon].append(int(taxon_id2n_CDS[taxon]))
        taxon_id2barplot_data[taxon].append(round((float(taxon_id2CDS_with_COG[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        taxon_id2barplot_data[taxon].append(round((float(taxon_id2CDS_with_KO[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        taxon_id2barplot_data[taxon].append(round((float(taxon_id2CDS_with_EC[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        taxon_id2barplot_data[taxon].append(round((float(taxon_id2n_swiss_hits[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        taxon_id2barplot_data[taxon].append(round((float(taxon_id2n_Pfam[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        try:
            taxon_id2barplot_data[taxon].append(round((float(taxon_id2n_uniprot_map[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        except:
            taxon_id2barplot_data[taxon].append(0)
        try:
            taxon_id2barplot_data[taxon].append(round((float(taxon_id2n_string[taxon])/float(taxon_id2n_CDS[taxon]))*100,2))
        except:
            taxon_id2barplot_data[taxon].append(0)
    taxon2set2value_heatmap = {}

    sql = 'select taxon_id,operon_id from custom_tables_DOOR2_operons t1 ' \
          ' right join custom_tables_locus2seqfeature_id t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' where taxon_id is not NULL group by taxon_id;'

    sql_string_taxons = 'select taxon_id from string_interactions where taxon_id is not NULL group by taxon_id;'
    string_taxons = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql_string_taxons,)]
    taxon_id2string = {}
    taxon_id2door = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    sql_string_pmid = 'select taxon_id,count(*) from(select * from string_seqfeature_id2string_pmid ' \
                      ' group by seqfeature_id) A where taxon_id is not NULL group by taxon_id;'
    taxon_id2count_pmid = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_string_pmid,))

    for taxon in taxon_id2door:
        if taxon not in taxon_id2count_pmid:
            taxon_id2count_pmid[taxon] = 0
        if taxon in string_taxons:
            taxon_id2string[taxon] = 1
        else:
            taxon_id2string[taxon] = 0

        if taxon_id2door[taxon] is None:
            taxon_id2door[taxon] = 0
        else:
            taxon_id2door[taxon] = 1

    taxon2set2value_heatmap['Doors'] = taxon_id2door
    taxon2set2value_heatmap['STRING'] = taxon_id2string
    taxon2set2value_heatmap['PMID'] = taxon_id2count_pmid

    tree1, style1 = phylo_tree_bar.plot_tree_barplot(tree,
                      taxon_id2barplot_data,
                      ['nCDS','COG', 'KO', 'PRIAM', 'SwissProt', 'Pfam', 'UniProt-mapping','STRING-mapping'],
                      taxon2set2value_heatmap=taxon2set2value_heatmap,
                      header_list2=['Doors', 'STRING', 'PMID'],
                      presence_only=True,
                      biodb=biodb,
                      column_scale=True,
                      general_max=False,
                      barplot2percentage=[False, True, True, True, True, True, True, True])

    #t.render("test2.svg", tree_style=ts)
    path = settings.BASE_DIR + '/assets/temp/tree.svg'
    asset_path = '/temp/tree.svg'

    proportions = []

    for taxon in taxon_id2n_CDS:

        proportions.append(["KO",float(taxon_id2CDS_with_KO[taxon])/float(taxon_id2n_CDS[taxon]), int(taxon_id2n_CDS[taxon])])
        proportions.append(["Hits SwissProt",float(taxon_id2n_swiss_hits[taxon])/float(taxon_id2n_CDS[taxon]), int(taxon_id2n_CDS[taxon])])
        proportions.append(["COG",float(taxon_id2CDS_with_COG[taxon])/float(taxon_id2n_CDS[taxon]), int(taxon_id2n_CDS[taxon])])

    pairwiseid_plots.plot_multiseries_points(proportions,output_path="/home/tpillone/ko2size.svg")

    tree1.render(path, dpi=800, tree_style=style1)
    return render(request, 'chlamdb/species_specific.html', my_locals(locals()))

def orthogroup_KO_COG(request):
    biodb = settings.BIODB
    from chlamdb.plots import pairwiseid_plots

    sql = 'select orthogroup,count(*) as n from (select orthogroup,COG_id from ' \
          ' COG_locus_tag2gi_hit t1 ' \
          ' inner join orthology_detail t2 on t1.locus_tag=t2.locus_tag ' \
          ' group by orthogroup,COG_id) A group by A.orthogroup'

    server, db = manipulate_biosqldb.load_db(biodb)

    group2n_COG = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    pairwiseid_plots.basic_plot(group2n_COG.values())


    sql = ' select orthogroup, count(*) as n from (select orthogroup,ko_id from enzyme_locus2ko ' \
          ' group by orthogroup,ko_id) A group by orthogroup'

    server, db = manipulate_biosqldb.load_db(biodb)

    group2n_KO = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    pairwiseid_plots.basic_plot(group2n_KO.values(),output_path="/home/tpillone/test2.svg")

    return render(request, 'chlamdb/species_specific.html', my_locals(locals()))

def paralogs(request):
    biodb = settings.BIODB

    from ete3 import Tree,TreeStyle
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.phylo_tree_display import phylo_tree_bar

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select * from genomes_info' % biodb

    genomes_data = server.adaptor.execute_and_fetchall(sql,)

    sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
               ' where t2.name="%s";' % biodb
    server, db = manipulate_biosqldb.load_db(biodb)
    tree = Tree(server.adaptor.execute_and_fetchall(sql_tree,)[0][0])
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()
    sql = 'select taxon_id, count(*) as n from comparative_tables_silix_35 group by taxon_id, silix_id'
    data = server.adaptor.execute_and_fetchall(sql,)
    sql = 'select taxon_id, count(*) as n from comparative_tables_silix_50 group by taxon_id, silix_id'
    data2 = server.adaptor.execute_and_fetchall(sql,)
    sql = 'select taxon_id, count(*) as n from comparative_tables_silix_80 group by taxon_id, silix_id'
    data3 = server.adaptor.execute_and_fetchall(sql,)
    sql = 'select taxon_id, count(*) as n from comparative_tables_silix_98 group by taxon_id, silix_id'
    data4 = server.adaptor.execute_and_fetchall(sql,)
    taxon_id2n_paralogs = {}

    for row in data:
        if row[0] not in taxon_id2n_paralogs:
            taxon_id2n_paralogs[row[0]] = [0,0,0,0]
        if int(row[1])>1:
            taxon_id2n_paralogs[row[0]][0] += int(row[1])-1
    for row in data2:
        if int(row[1])>1:
            taxon_id2n_paralogs[row[0]][1] += int(row[1])-1
    for row in data3:
        if int(row[1])>1:
            taxon_id2n_paralogs[row[0]][2] += int(row[1])-1
    for row in data4:
        if int(row[1])>1:
            taxon_id2n_paralogs[row[0]][3] += int(row[1])-1

    sql = 'select taxon_id,sum(n_CDS) from genomes_info t1 ' \
          ' inner join bioentry t2 on t1.accession=t2.accession ' \
          ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
          ' where t3.name="%s" group by taxon_id;' % (biodb)

    data = server.adaptor.execute_and_fetchall(sql,)

    for row in data:
        taxon_id2n_paralogs[row[0]].append(int(row[1]))
        taxon_id2n_paralogs[row[0]].append(round((taxon_id2n_paralogs[row[0]][0]/float(row[1]))*100,2))
        taxon_id2n_paralogs[row[0]].append(round((taxon_id2n_paralogs[row[0]][1]/float(row[1]))*100,2))
        #taxon_id2n_paralogs[row[0]].append(round((taxon_id2n_paralogs[row[0]][2]/float(row[1]))*100,2))

    tree1, style1 = phylo_tree_bar.plot_tree_barplot(tree,
                      taxon_id2n_paralogs,
                      ['silix_35','silix_50','silix_80', 'silix_98','TOTAL', 'perc_35', 'perc_50'],
                      taxon2set2value_heatmap=False,
                      header_list2=False,
                      presence_only=True,
                      biodb=biodb,
                      column_scale=True,
                      general_max=False)

    #t.render("test2.svg", tree_style=ts)
    path = settings.BASE_DIR + '/assets/temp/tree.svg'
    asset_path = '/temp/tree.svg'

    tree1.render(path, dpi=800, tree_style=style1)
    return render(request, 'chlamdb/species_specific.html', my_locals(locals()))


def species_specific_groups(request):
    
    print("species_specific_groups -------------")
    
    biodb = settings.BIODB

    from ete3 import Tree,TreeStyle
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    from chlamdb.phylo_tree_display import phylo_tree_bar
    from chlamdb.phylo_tree_display import species_tree

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select * from genomes_info'

    genomes_data = server.adaptor.execute_and_fetchall(sql,)

    sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
               ' where t2.name="%s";' % biodb
    server, db = manipulate_biosqldb.load_db(biodb)
    tree = Tree(server.adaptor.execute_and_fetchall(sql_tree,)[0][0])
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)

    sql = 'select taxon_id, family_description from custom_tables_taxonomy'
    sql2 = 'select family_description, taxon_id from custom_tables_taxonomy'
    taxon_id2species_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    species_id2taxon_id = {}
    data = server.adaptor.execute_and_fetchall(sql2,)
    for row in data:
        if row[0] not in species_id2taxon_id:
            species_id2taxon_id[row[0]] = [str(row[1])]
        else:
            species_id2taxon_id[row[0]].append(str(row[1]))

    all_taxons = [str(i.name) for i in tree.iter_leaves()]
    # changing taxon id to species id
    for leaf in tree.iter_leaves():
        #print '%s --> %s' % (leaf.name, str(taxon_id2species_id[str(leaf.name)]))
        leaf.name = "%s" % str(taxon_id2species_id[str(leaf.name)])

    # attributing unique id to each node
    # if all node descendant have the same name, use that name as node name
    n = 0
    for node in tree.traverse():
        if node.name=='':
            desc_list = list(set([i.name for i in node.iter_descendants()]))
            try:
                desc_list.remove('')
            except ValueError:
                pass
            if len(desc_list) != 1:
                node.name = '%sbb' % n
            else:
                node.name = desc_list[0]
            n+=1

    node2labels = tree.get_cached_content(store_attr="name")

    def collapsed_leaf(node):
        if len(node2labels[node]) == 1:
            return True
        else:
            return False
    species_id2count_unique = {}
    
    for species in species_id2taxon_id:
        #if 'Akkermansia' in species:
        #    species_id2count_unique[species] = [0]
        #    continue
        species_taxons = species_id2taxon_id[species]
        other_taxons = set(all_taxons) - set(species_taxons)
        print("extract subtable")
        '''
        mat, mat_all = biosql_own_sql_tables.get_comparative_subtable(biodb,
                                                                     "orthology",
                                                                     "orthogroup",
                                                                     species_taxons,
                                                                     other_taxons,
                                                                     ratio=1/float(len(species_taxons)),
                                                                     single_copy=False,
                                                                     accessions=False,
                                                                     cache=cache)

        '''
        species_id2count_unique[species] = [1]
   
    # plot tree
    tree1, style1 = phylo_tree_bar.plot_tree_barplot(tree_species,
                    species_id2count_unique,
                    ['unique'],
                    taxon2set2value_heatmap=False,
                    header_list2=False,
                    presence_only=True,
                    biodb=biodb,
                    column_scale=True,
                    general_max=False)

    path = settings.BASE_DIR + '/assets/temp/stree.svg'
    asset_path = '/temp/stree.svg'

    tree1.render(path, dpi=800, tree_style=style1)
    return render(request, 'chlamdb/species_specific.html', my_locals(locals()))


def api_vs_16S_identity(request):
    biodb = settings.BIODB
    from chlamdb.plots import pairwiseid_plots

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select t1.taxon_1, t1.taxon_2, median_identity,identity from ' \
          ' comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 ' \
          ' on t1.taxon_1=t2.taxon_1 and t1.taxon_2=t2.taxon_2;'
    data = server.adaptor.execute_and_fetchall(sql)

    api = []
    rrna_identity = []
    for row in data:
        api.append(float(row[2]))
        rrna_identity.append(float(row[3]))

    pairwiseid_plots.basic_plot(api, rrna_identity, output_path="~/tata2.svg")

    # intera order comparisons only

    taxon_complete_genomes = [46,
87925,
1279815,
1279822,
52,
49,
804807,
886707,
314,
1280098,
48,
55,
1280079,
293,
60,
64,
66
]
    taxon_complete_genomes = [1279774, # exclde chlamydiaceae
1280030,
1280079,
1280091,
290,
1280035,
1280044,
293,
1117985,
1035343,
1280085,
60,
288,
315,
1280034,
1069693,
1172027,
1172028,
1279839,
1279497,
64,
66,
67,
1137444
]


    taxon_complete_genomes = [str(i) for i in taxon_complete_genomes]
    filter_complete_genomes = ','.join(taxon_complete_genomes)
    identity_column = 'max_density'#'median_identity'#

    sql = 'select B.order_description,A.%s,identity from' \
          ' (select t1.taxon_1, t1.taxon_2, %s,identity ' \
          ' from comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 on t1.taxon_1=t2.taxon_1 ' \
          ' and t1.taxon_2=t2.taxon_2) A ' \
          ' inner join custom_tables_taxonomy B ' \
          ' on A.taxon_1=B.taxon_id ' \
          ' inner join custom_tables_taxonomy C on A.taxon_2=C.taxon_id ' \
          ' where B.order_description=C.order_description;' % (identity_column,
                                                               identity_column)

    data = server.adaptor.execute_and_fetchall(sql)
    api = []
    rrna_identity = []
    for row in data:
        api.append(float(row[1]))
        rrna_identity.append(float(row[2]))


    pairwiseid_plots.basic_plot(api, rrna_identity, output_path="~/tata3.svg")

    data_list = [list(i) for i in data]

    sql = 'select B.order_description,A.%s,identity from' \
          ' (select t1.taxon_1, t1.taxon_2, %s,identity ' \
          ' from comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 on t1.taxon_1=t2.taxon_1 ' \
          ' and t1.taxon_2=t2.taxon_2) A ' \
          ' inner join custom_tables_taxonomy B ' \
          ' on A.taxon_1=B.taxon_id ' \
          ' inner join custom_tables_taxonomy C on A.taxon_2=C.taxon_id ' \
          ' where B.order_description=C.order_description;' % (identity_column,
                                                               identity_column)
    data_chlam = [list (i) for i in server.adaptor.execute_and_fetchall(sql)]

    data_list+=data_chlam

    pairwiseid_plots.plot_multiseries_points(data_list, "~/tata4")

    # inter family same order
    sql = 'select B.order_description,A.%s,identity from' \
          ' (select t1.taxon_1, t1.taxon_2, %s,identity ' \
          ' from comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 on t1.taxon_1=t2.taxon_1 ' \
          ' and t1.taxon_2=t2.taxon_2) A ' \
          ' inner join custom_tables_taxonomy B ' \
          ' on A.taxon_1=B.taxon_id ' \
          ' inner join custom_tables_taxonomy C on A.taxon_2=C.taxon_id ' \
          ' where B.order_description=C.order_description ' \
          ' and B.family_description!=C.family_description;' % (identity_column,
                                                                identity_column)


    data = server.adaptor.execute_and_fetchall(sql)
    api = []
    rrna_identity = []
    for row in data:
        api.append(float(row[1]))
        rrna_identity.append(float(row[2]))

    data_list = [list(i) for i in data]

    sql0 = 'select B.order_description,A.%s,identity from' \
          ' (select t1.taxon_1, t1.taxon_2, %s,identity ' \
          ' from comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 on t1.taxon_1=t2.taxon_1 ' \
          ' and t1.taxon_2=t2.taxon_2) A ' \
          ' inner join custom_tables_taxonomy B ' \
          ' on A.taxon_1=B.taxon_id ' \
          ' inner join custom_tables_taxonomy C on A.taxon_2=C.taxon_id ' \
          ' where B.order_description=C.order_description ' \
          ' and B.family_description!=C.family_description;' % (identity_column,
                                                               identity_column)

    sql = 'select B.order_description,A.%s,identity from' \
          ' (select t1.taxon_1, t1.taxon_2, %s,identity ' \
          ' from comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 on t1.taxon_1=t2.taxon_1 ' \
          ' and t1.taxon_2=t2.taxon_2 where t1.taxon_1 in (%s) and t1.taxon_2 in (%s)) A ' \
          ' inner join custom_tables_taxonomy B ' \
          ' on A.taxon_1=B.taxon_id ' \
          ' inner join custom_tables_taxonomy C on A.taxon_2=C.taxon_id ' \
          ' where B.order_description=C.order_description ' \
          ' and B.family_description!=C.family_description;' % (identity_column,
                                                                identity_column,
                                                                filter_complete_genomes,
                                                                filter_complete_genomes)
          
    data_chlam = [list (i) for i in server.adaptor.execute_and_fetchall(sql)]
    data_chlam0 = [list (i) for i in server.adaptor.execute_and_fetchall(sql0)]

    data_list+=data_chlam

    pairwiseid_plots.plot_multiseries_points(data_list, "~/tata5")

    sql = 'select B.phylum_description,A.%s,identity from' \
          ' (select t1.taxon_1, t1.taxon_2, %s,identity ' \
          ' from comparative_tables_reciprocal_BBH_average_identity t1 ' \
          ' inner join custom_tables_identity_16S t2 on t1.taxon_1=t2.taxon_1 ' \
          ' and t1.taxon_2=t2.taxon_2) A ' \
          ' inner join custom_tables_taxonomy B ' \
          ' on A.taxon_1=B.taxon_id ' \
          ' inner join custom_tables_taxonomy C on A.taxon_2=C.taxon_id ' \
          ' where B.phylum_description=C.phylum_description and B.order_description!=C.order_description ' \
          ' and B.order_description in ("Rhizobiales","Rickettsiales")' \
          ' and C.order_description in ("Rhizobiales","Rickettsiales");' % (identity_column,
                                                                            identity_column)
    data = server.adaptor.execute_and_fetchall(sql)
    api = []
    rrna_identity = []
    for row in data:
        api.append(float(row[1]))
        rrna_identity.append(float(row[2]))

    data_list = [list(i) for i in data]

    pairwiseid_plots.plot_multiseries_points(data_list, "~/inter_orders")

    return render(request, 'chlamdb/species_specific.html', my_locals(locals()))


def circos_main(request):
    biodb = settings.BIODB
    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos

    if request.method == 'POST':
        
        reference_taxon = request.POST["reference_taxon"]
        target_taxons = eval(request.POST["target_list"])
        highlight = eval(request.POST["highlight"])
        
        if len(target_taxons) == 0:
            try:
                sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa where taxon_1=%s order by identity desc;' % (reference_taxon)
                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
            except:
                sql_order = 'select taxon_2 from comparative_tables_shared_orthogroups where taxon_1=%s order by n_shared_orthogroups DESC;' % (reference_taxon)

                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]   
                         
        task = run_circos_main.delay(reference_taxon, target_taxons, highlight)
        print("task", task)
        task_id = task.id
        envoi_circos = True
        envoi_region = True
            
    if request.method == 'GET':
        server, db = manipulate_biosqldb.load_db(biodb)

        reference_taxon = int(request.GET.getlist('ref')[0])

        if request.GET.getlist('t')[0] == '':
            # if no target list given, get the 10 closest genomes
            try:
                sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa where taxon_1=%s order by identity desc;' % (reference_taxon)
                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
            except:
                sql_order = 'select taxon_2 from comparative_tables_shared_orthogroups where taxon_1=%s order by n_shared_orthogroups DESC;' % (reference_taxon)

                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
        else:
            target_taxons = [int(i) for i in request.GET.getlist('t')]
        highlight = request.GET.getlist('h')

        #sql = 'select locus_tag,traduction from orthology_detail_k_cosson_05_16 where orthogroup in (%s) and accession="NC_016845"' % ('"'+'","'.join(highlight)+'"')

        task = run_circos_main.delay(reference_taxon, target_taxons, highlight)
        print("task", task)
        task_id = task.id
        envoi_circos = True
        envoi_region = True

        #return HttpResponse(json.dumps({'task_id': task.id}), content_type='application/json')

    return render(request, 'chlamdb/circos_main.html', my_locals(locals()))


def circos_blastnr(request):
    biodb = settings.BIODB
    from chlamdb.plots import gbk2circos
    circos_form_class = make_genome_selection_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST':

        form = circos_form_class(request.POST)

        if form.is_valid():

            reference_taxon = form.cleaned_data['genome']

            description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

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


            ref_name = ''
            for i in reference_accessions:
                ref_name += i
            circos_file = "circos/%s.svg" % ref_name
            from chlamdb.biosqldb import circos
            from chlamdb.biosqldb import shell_command
            import ete3



            outtput_dir = settings.BASE_DIR + "/assets/circos/"

            myplot = circos.CircosAccession2blastnr_plot(server,
                                                         biodb,
                                                         record_list,
                                                         outtput_dir,
                                                         taxon_list=[],
                                                         highlight_BBH=True)

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

            circos_new_file = '/assets/circos/circos_clic_%s.html' % ref_name

            with open(settings.BASE_DIR + circos_new_file, "w") as f:
                f.write(circos_html)

            #target_map_file = settings.BASE_DIR + "/templates/circos/%s.html" % ref_name
            circos_svg_file = "circos/%s.svg" % ref_name
            circos_png_file = "circos/%s.png" % ref_name
            original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
            #target_map_file_svg = settings.BASE_DIR + "/templates/circos/%s.svg" % ref_name
            map_file = "circos/%s.html" % ref_name
            svg_file = "circos/%s.svg" % ref_name
            #a, b, c = shell_command.shell_command("mv %s %s" % (original_map_file, target_map_file))
            #a, b, c = shell_command.shell_command("cp %s %s" % (original_map_file_svg, target_map_file_svg))
            map_name = ref_name
            envoi_circos = True
            envoi_region = True
    else:
        form = circos_form_class()
    return render(request, 'chlamdb/circos_blastnr.html', my_locals(locals()))


def circos(request):
    biodb = settings.BIODB

    circos_form_class = make_circos_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST':

        form = circos_form_class(request.POST)

        if form.is_valid():
            reference_taxon = form.cleaned_data['circos_reference']
            target_taxons = form.cleaned_data['targets']
            task = run_circos.delay(reference_taxon, target_taxons)
            print("task", task)
            return HttpResponse(json.dumps({'task_id': task.id}), content_type='application/json')
        else:
            return HttpResponse(json.dumps({'task_id': None}), content_type='application/json')
    else:
        form = circos_form_class()
    
    local_vars = my_locals(locals())
    local_vars["form"] = form
    
    return render(request, 'chlamdb/circos.html', local_vars)


def alignment(request, input_fasta):

    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        pass
    return render(request, 'chlamdb/alignment.html', my_locals(locals()))






def format_seqfeature_values(server, seqfeature_id):
    biodb = settings.BIODB
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



def search_taxonomy(request):
    biodb = settings.BIODB
    from collections import Counter
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST': 

        if len(request.POST['Phylum']) == 0:
            genome_accession = request.POST['Genome']
            superkingdom = request.POST['Superkingdom'].split('_')[-1]
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
                          ' from blastnr_blastnr_hits_%s as t1 ' \
                          ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s"' % (genome_accession,
                                                        genome_accession,
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
                          ' from blastnr_blastnr_hits_%s as t1 ' \
                          ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s"' % (genome_accession,
                                                                                genome_accession,
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
                          ' from blastnr_blastnr_hits_%s as t1 ' \
                          ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s"' % (genome_accession,
                                                                                genome_accession,
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
                          ' from blastnr_blastnr_hits_%s as t1 ' \
                          ' inner join blastnr_blastnr_hits_taxonomy_filtered_%s as t2 on t1.nr_hit_id = t2.nr_hit_id ' \
                          ' inner join blastnr_blastnr_taxonomy as t3 on t2.subject_taxon_id = t3.taxon_id' \
                          ' inner join blastnr_blastnr_hsps_%s as t4 on t1.nr_hit_id=t4.nr_hit_id' \
                          ' where t1.hit_number=%s and t3.superkingdom="%s" and t3.phylum="%s";' % (genome_accession,
                                                                                                    genome_accession,
                                                                                                    genome_accession,
                                                                                                    1,
                                                                                                    superkingdom,
                                                                                                    phylum)

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

        '''
        envoi = True


    return render(request, 'chlamdb/search_taxonomy.html', my_locals(locals()))



def interpro(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    interproform = make_interpro_from(biodb)

    if request.method == 'POST': 

        form = interproform(request.POST)
        #form2 = ContactForm(request.POST)
        if form.is_valid():  
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

                    sql = 'select %s from interpro where %s and (interpro_description REGEXP "%s" or signature_description REGEXP "%s")' % (columns, 
                                                                                                                                            taxon_limit, 
                                                                                                                                            search_term, 
                                                                                                                                            search_term)

                if search_type == "GO":
                    sql = 'select %s from interpro where %s and (GO_terms REGEXP "%s")' % (columns, 
                                                                                           taxon_limit, 
                                                                                           search_term)


                if search_type == "EC":
                    sql = 'select %s from interpro where %s and (pathways REGEXP "%s")' % (columns, 
                                                                                           taxon_limit, 
                                                                                           search_term)

                if search_type == "interpro_accession":
                    sql = 'select %s from interpro where %s and (interpro_accession REGEXP "%s")' % (columns, 
                                                                                                     taxon_limit, 
                                                                                                     search_term)


                try:
                    raw_data = server.adaptor.execute_and_fetchall(sql, )
                except IndexError:
                    invalid_id = True


            envoi = True

    else:  
        form = interproform()

    return render(request, 'chlamdb/interpro.html', my_locals(locals()))


def search(request):

    import re
    from chlamdb.biosqldb import search

    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST':  # S'il s'agit d'une requête POST
        display_from = 'yes'
        form = SearchForm(request.POST)
        #form2 = ContactForm(request.POST)
        if form.is_valid():  
            invalid_id = False
            # Ici nous pouvons traiter les données du formulaire
            search_type = form.cleaned_data['search_type']
            search_term = form.cleaned_data['search_term']
            #biodb = form.cleaned_data['biodatabase']

            search_result = search.perform_search(request, search_term)

            if search_result is None:
                search_echec = True 
            else:
                locus_list, 
                raw_data_gene_raw_data_product, 
                raw_data_EC,
                raw_data_ko,
                raw_data_cog,
                raw_data_interpro,
                raw_data_pathway,
                raw_data_pmid,
                pmid2n_homologs,
                pmid_data,
                search_echec = search_result
            envoi = True

    else:  
        search_term = request.GET.get('accession').strip()
        if search_term[-1] == '+':
            search_term = search_term[0:-1]
        elif "%" in search_term:
            import re
            search_term = re.sub("%", " ", search_term)
        elif "@" in search_term:
            import re
            search_term = re.sub("@", " ", search_term)

        if search_term:

            search_result = search.perform_search(request, biodb, search_term)
            print(search_result)
            if search_result is None:
                search_echec = True 
            elif isinstance(search_result, HttpResponse):
                return search_result

            elif len(search_result) == 1:
                return locusx(request, search_result[0])
            else:
                locus_list, raw_data_EC, raw_data_ko, raw_data_cog, raw_data_interpro, raw_data_pathway, raw_data_pmid, pmid2n_homologs, pmid_data, synonymous = search_result
            envoi = True

            envoi = True
            display_form = "no"
            #print "display_form",  display_form
            #search_result = perform_search(locus, False)
            #if isinstance(search_result, HttpResponse):
            #    return search_result
            #else:
            #    envoi = True
        else:
            display_from = 'yes'
            form = SearchForm()

    search_term_edit = re.sub("\-|\+", "", search_term)
    print("SEARCH:", search_term_edit, search_term)

    return render(request, 'chlamdb/search.html', my_locals(locals()))




def primer_search(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    if request.method == 'POST': 

        form = PCRForm(request.POST)

        if form.is_valid():  
            from Bio.Blast.Applications import NcbiblastpCommandline
            #from StringIO import StringIO
            from tempfile import NamedTemporaryFile


            from Bio.Alphabet import IUPAC
            import os
            from chlamdb.biosqldb import shell_command
            import re
            def ExtractAlphanumeric(InputString):
                from string import ascii_letters, digits
                return "".join([ch for ch in InputString if ch in (ascii_letters + digits)])

            input_sequence = form.cleaned_data['blast_input']
            input_sequence = ExtractAlphanumeric(input_sequence)

            input_sequence = input_sequence.rstrip(os.linesep)

            my_record = SeqRecord(Seq(input_sequence, IUPAC.protein), id="INPUT", description="INPUT")

            query_file = NamedTemporaryFile()
            SeqIO.write(my_record, query_file, "fasta")
            query_file.flush()

            blastdb = settings.BASE_DIR + '/assets/blast_db/%s.faa' % biodb


            blastp_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=0.001, outfmt=0)
            stdout, stderr = blastp_cline()

            blast_file = NamedTemporaryFile()
            blast_file.write(stdout)
            mview_cmd = 'mview -in blast -ruler on -html data -css on -coloring identity %s' % blast_file.name
            stdout, stderr, code = shell_command.shell_command(mview_cmd)
            blast_result = stdout


            envoi = True

    else:  
        form = PCRForm()

    return render(request, 'chlamdb/pcr.html', my_locals(locals()))



def motif_search(request):
    biodb = settings.BIODB
    from chlamdb.biosqldb import shell_command
    import os
    server, db = manipulate_biosqldb.load_db(biodb)

    motif_form_class = make_motif_form(biodb)

    if request.method == 'POST':

        form = motif_form_class(request.POST)

        if form.is_valid():

            from Bio.Emboss.Applications import FuzznucCommandline
            from tempfile import NamedTemporaryFile

            input_pattern = form.cleaned_data['motif_input']
            n_missmatch = form.cleaned_data['n_missmatch']
            target_taxon_id = form.cleaned_data['search_in']

            if target_taxon_id != "all":
                accessions = manipulate_biosqldb.taxon_id2accessions(server, target_taxon_id, biodb)

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

            cmd = 'fuzzpro -sequence %s -pmismatch %s -rformat seqtable -auto -stdout -pattern "%s"' % (db_path ,n_missmatch ,input_pattern)

            std_out, std_err, code = shell_command.shell_command(cmd)
            #fuzznuc_cline = FuzznucCommandline(sequence=genome_db, mismatch=n_missmatch, pattern=input_pattern, stdout=True)#, rformat="srspair")
            #stdout,stderr = fuzznuc_cline()


            envoi = True

    else:  
        form = motif_form_class()  # empty form

    return render(request, 'chlamdb/motifs.html', my_locals(locals()))


def blast_profile(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST': 

        form = BlastProfileForm(request.POST, request.FILES)

        if form.is_valid():  
            from tempfile import NamedTemporaryFile
            from io import StringIO
            from chlamdb.phylo_tree_display import biosqldb_plot_blast_hits_phylo
            from chlamdb.biosqldb import biosql_own_sql_tables
            from ete3 import Tree, TreeStyle

            fasta_file = request.FILES['fasta_file']
            fasta_string = StringIO(request.FILES['fasta_file'].read().decode("UTF-8"))
            fasta_rec = [i for i in SeqIO.parse(fasta_string, 'fasta')]

            try:
                ordered_labels_all = [i.id.split('|')[1] for i in fasta_rec]
            except:
                ordered_labels_all = [i.id for i in fasta_rec]

            blast_type = form.cleaned_data['blast']

            #my_record = SeqIO.read(request.FILES['fasta_file'].open(), 'fasta')


            tree, style1, tree2, style2, tree3, style3, locus2taxon2locus_closest = biosqldb_plot_blast_hits_phylo.plot_BBH_phylo(fasta_rec, 
                                                                                                                                  biodb,
                                                                                                                                  settings.BASE_DIR + '/assets/')


            sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

            tt = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

            t1 = Tree(tt)

            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            t1.ladderize()
            ordered_taxons_all = [str(i.name) for i in t1.iter_leaves()]

            sql = 'select taxon_id, t2.description from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                  ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % (biodb)

            taxon2genome = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            taxon_list = taxon2genome.keys()


            ordered_taxons_keep=[]
            for taxon in ordered_taxons_all:
                if taxon in taxon_list:
                    ordered_taxons_keep.append(taxon)

            ordered_labels_keep=[]
            for label in ordered_labels_all:
                if label in locus2taxon2locus_closest:
                    ordered_labels_keep.append(label)

            path = settings.BASE_DIR + '/assets/temp/profile_tree.svg'
            path2 = settings.BASE_DIR + '/assets/temp/profile_tree2.svg'
            path3 = settings.BASE_DIR + '/assets/temp/profile_tree3.svg'
            asset_path = '/temp/profile_tree.svg'
            asset_path2 = '/temp/profile_tree2.svg'
            asset_path3 = '/temp/profile_tree3.svg'

            tree3.render(path3, dpi=500, tree_style=style3)
            tree.render(path2, dpi=500, tree_style=style1)
            tree2.render(path, dpi=500, tree_style=style2)


            envoi = True

    else:  
        form = BlastProfileForm()

    return render(request, 'chlamdb/blast_profile.html', my_locals(locals()))








def blast(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    #server = manipulate_biosqldb.load_db(db)
    blast_form_class = make_blast_form(db) 

    if request.method == 'POST': 

        form = blast_form_class(request.POST)

        if form.is_valid():  
            input_sequence = form.cleaned_data['blast_input']
            customized_evalue = form.cleaned_data['evalue']
            number_blast_hits = form.cleaned_data['max_number_of_hits']
            target_accession = form.cleaned_data['target'] #form.cleaned_data returns a dictionary of validated form input fields and their values, where string primary keys are returned as objects 
            blast_type = form.cleaned_data['blast']

          
            unknown_format = False

            if '>' in input_sequence:
                my_record = [i for i in SeqIO.parse(StringIO(input_sequence), 'fasta')]
                seq = my_record[0]
                                          
            else:
                input_sequence == input_sequence.rstrip(os.linesep)
                seq = Seq(input_sequence)
                                
            dna = set("ATGCNRYKMSWBDHV") #added ambiguous nucleotides from CoGePEDIA
            prot = set('ACDEFGHIKLMNPQRSTVWYXZJOU') #added ambiguous aa 
            check_seq_DNA = set(seq) - dna
            check_seq_prot = set(seq) - prot
            
            if  not check_seq_DNA:
                try: my_record[0].description = "DNA"
                except: my_record = [SeqRecord(seq, id="INPUT", description="DNA")]
                seq_type= my_record[0].description
                    
            elif  not check_seq_prot:
                try: my_record[0].description = "Protein"
                except: my_record = [SeqRecord(seq, id="INPUT", description="Protein")]
                seq_type= my_record[0].description
                print('my_record_changed', my_record)
            else:
                unknown_format = True
                                                   
                
            if not unknown_format:
                if seq_type == 'DNA' and blast_type in ["blastp", "tblastn"]:
                    wrong_format = True
                elif seq_type =='Protein' and blast_type in ["blastn_ffn", "blastn_fna", "blastx"]:
                    wrong_format = True
                else:
                    query_file = NamedTemporaryFile(mode='w')
                    SeqIO.write(my_record, query_file, "fasta")
                    query_file.flush()
                    
                    
                    if target_accession =='all':
                        key_dict = 'merged'
                    else:
                        dictionary_acc_names=db.get_taxon_id_to_filenames() #here db replaces 'db_utils.DB' that is used to create the link because it has already been done
                        key_dict=dictionary_acc_names[int(target_accession)]               
                        
                    if number_blast_hits=='all':
                       number_blast_hits = 100000                 
                    if blast_type=='blastn_ffn':
                        blastType = 'locus'
                        blastdb = settings.FOLDER_PATH + "blast_DB/ffn/%s" % (key_dict)
                        blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=customized_evalue, max_target_seqs=number_blast_hits, outfmt=0 )
                        print('BLASTDB', blastdb)
                    if blast_type=='blastn_fna':
                        blastType = 'genome'
                        blastdb = settings.FOLDER_PATH  + "blast_DB/fna/%s" % (key_dict)
                        print(blastdb)
                        blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=customized_evalue, max_target_seqs=number_blast_hits,  outfmt=0)
                    if blast_type=='blastp':
                        blastType = 'locus'
                        blastdb = settings.FOLDER_PATH  + "blast_DB/faa/%s" % (key_dict)
                        blast_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=customized_evalue, max_target_seqs=number_blast_hits, outfmt=0)
                    if blast_type=='tblastn':
                        blastType = 'genome'
                        blastdb = settings.FOLDER_PATH  + "blast_DB/fna/%s" % (key_dict)
                        blast_cline = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=customized_evalue, max_target_seqs=number_blast_hits,  outfmt=0)
                        blast_cline2 = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=customized_evalue, max_target_seqs=number_blast_hits, outfmt=5)
                    if blast_type=='blastx':
                        blastType = 'locus'
                        blastdb = settings.FOLDER_PATH  + "blast_DB/faa/%s" % (key_dict)
                        blast_cline = NcbiblastxCommandline(query=query_file.name, db=blastdb, evalue=customized_evalue, max_target_seqs=number_blast_hits,  outfmt=0) #max_target_seqs=number_blast_hits
                    
                    
                    blast_stdout, blast_stderr = blast_cline()
                                    

                    if blast_type=='tblastn':
                        from Bio.SeqUtils import six_frame_translations

                        blast_stdout2, blast_stderr2 = blast_cline2()
                        print('blast_stdout2', blast_stdout2)

                        blast_records = NCBIXML.parse(StringIO(blast_stdout2))
                        print('blast_records', blast_records)
                        all_data = []
                        best_hit_list = []
                        for record in blast_records:
                            print(record)
                            for n, alignment in enumerate(record.alignments): #n is a increasing number given from 0 going on to each allignment considering at n=0 the best hit
                                accession = alignment.title.split(' ')[0] #before it was 1 but it went to the beginning of the name of the sample, with 0 we take the contig name
                                description=alignment.title.replace(accession, '') #not sure if this description is what it wanted before
                                
                                for n2, hsp in enumerate(alignment.hsps): #all n2 are 0
                                    if n == 0 and n2 == 0:  #select the best hit
                                        best_hit_list.append([alignment.title.split(' ')[0], hsp.sbjct_start, hsp.sbjct_end])
                                    start = hsp.sbjct_start
                                    end = hsp.sbjct_end
                                    if start > end:
                                        start = hsp.sbjct_end
                                        end = hsp.sbjct_start
                                    length = end-start
                                    seq_A = db.location2sequence(accession, start, end) #replaced biodb wtih db
                                    anti = reverse_complement(seq_A)
                                    comp = anti[::-1]
                                    length = len(seq_A)
                                    frames = {}
                                    for i in range(0, 3):
                                        fragment_length = 3 * ((length-i) // 3)
                                        tem1 = translate(seq_A[i:i+fragment_length])
                                        frames[i+1] = '<span style="color: #181407;">%s</span><span style="color: #bb60d5;">%s</span><span style="color: #181407;">%s</span>' % (tem1[0:100], tem1[100:len(tem1)-99], tem1[len(tem1)-99:])
                                        tmp2 = translate(anti[i:i+fragment_length])[::-1]
                                        frames[-(i+1)] = tmp2

                                    all_data.append([accession, start, end, length, frames[1], frames[2], frames[3], frames[-1], frames[-2], frames[-3], description, seq_A])   #all_data is required in the hyml file to display the second graph
                                    print('all_data', all_data)

                        if len(best_hit_list) > 0:
                            print('best_hit_list', best_hit_list)
                            fig_list = []
                            for best_hit in best_hit_list:
                                accession = best_hit[0]
                                best_hit_start = best_hit[1]
                                best_hit_end = best_hit[2]
                                temp_location = os.path.join(settings.BASE_DIR, "assets/temp/")
                                temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
                                name = 'temp/' + os.path.basename(temp_file.name)
                                fig_list.append([accession, name])
                                orthogroup_list = db.location2plot(accession, #it wants the contig name, but it should take the locus tag
                                                                    temp_file.name,
                                                                    best_hit_start-15000,
                                                                    best_hit_end+15000,
                                                                    cache,
                                                                    color_locus_list = [],
                                                                    region_highlight=[best_hit_start, best_hit_end])
                                print('orthogroup_list', orthogroup_list)
                             


                    no_match = re.compile('.* No hits found .*', re.DOTALL) #I still do not know how it knows that there are no hits bit it works properly
                    

                    if no_match.match(blast_stdout):
                        print ("no blast hit")
                        blast_no_hits = blast_stdout
                    elif len(blast_stderr) != 0:
                        print ("blast error")
                        blast_err = blast_stderr #linked to the html file
                    else:

                        rand_id = id_generator(6) #it generates a random id of 6 character
                        blast_file_l = settings.BASE_DIR + '/assets/temp/%s.xml' % rand_id #this file changes every time you run it and it contians the blast output
                        f = open(blast_file_l, 'w')
                        f.write(blast_stdout)
                        #print('stout', blast_stdout)
                        #print('f', f)
                        f.close()
                        
                        asset_blast_path = '/assets/temp/%s.xml' % rand_id
                        js_out = True #it could be for java script

            envoi= True #here the data are passed when it is done correctly

    else:  
        form = blast_form_class()

    return render(request, 'chlamdb/blast.html', my_locals(locals()))


def get_record_from_memory(biodb, cache_obj, record_key, accession):

        biorecord = cache_obj.get(record_key)
        if not biorecord:
            new_record = biodb.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
            cache_obj.set(record_key, biorecord)
        return biorecord






def mummer(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    mummer_form_class = make_mummer_form(biodb)

    if request.method == 'POST': 
        plot = True
        form = mummer_form_class(request.POST)
        if form.is_valid():  
            server, db = manipulate_biosqldb.load_db(biodb)
            reference_taxon = form.cleaned_data['reference_genome']
            query_taxon = form.cleaned_data['query_genome']

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)
            query_accessions = manipulate_biosqldb.taxon_id2accessions(server, query_taxon, biodb)

            ref_accession = manipulate_biosqldb.taxon_id2chromosome_accession(server, biodb, reference_taxon)
            query_accession = manipulate_biosqldb.taxon_id2chromosome_accession(server, biodb, query_taxon)

            reference_path = settings.BASE_DIR + '/assets/%s/fna/%s.fna' % (biodb, ref_accession)
            query_path = settings.BASE_DIR + '/assets/%s/fna/%s.fna' % (biodb, query_accession)

            rand = id_generator(5)

            out_delta = settings.BASE_DIR + '/assets/temp/promer_%s' % rand
            out_plot = settings.BASE_DIR + '/assets/temp/promer_%s' % rand




            cmd1 = 'promer -l 2 -p %s %s %s' % (out_delta, reference_path, query_path)
            cmd2 = 'mummerplot -layout -small -png -p %s %s.delta' % (out_plot, out_delta)

            from chlamdb.biosqldb.shell_command import shell_command

            out, err, log = shell_command(cmd1)
            out, err, log = shell_command(cmd2)

            plot_path = 'temp/promer_%s.png' % rand

            if not os.path.exists(settings.BASE_DIR + '/assets/' + plot_path):
                plot = False

            envoi = True

    else:  
        form = mummer_form_class()

    return render(request, 'chlamdb/mummer.html', my_locals(locals()))






def circos2genomes(request):
    biodb = settings.BIODB
    from chlamdb.biosqldb import circos
    from chlamdb.biosqldb import shell_command
    server, db = manipulate_biosqldb.load_db(biodb)
    circos2genomes_form_class = make_circos2genomes_form(biodb)



    if request.method == 'POST': 
        #make_circos2genomes_form

        form = circos2genomes_form_class(request.POST)
        if form.is_valid():  
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

            orthogroup_list = []
            if len(protein_locus_list[0]) > 0:
                for protein in protein_locus_list:

                    sql = 'select orthogroup from orthology_detail where protein_id="%s" or locus_tag="%s"' % (protein, 
                                                                                                               protein)

                    try:
                        protein_group = server.adaptor.execute_and_fetchall(sql,)[0][0]
                        orthogroup_list.append(protein_group)
                    except IndexError:
                        valid_id = False

            if valid_id:

                #accession2description = manipulate_biosqldb.accession2description_dict(server, biodb)
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

                from chlamdb.biosqldb import circos

                path = settings.BASE_DIR + "/assets/circos"



                biplot = circos.CircosAccession2biplot(server, db, biodb, reference_records, query_records,
                                                       orthogroup_list, path)

                reference_file = "circos/%s" % biplot.reference_circos
                query_file = "circos/%s" % biplot.query_circos

            envoi = True

    else:  
        form = circos2genomes_form_class()

    return render(request, 'chlamdb/circos2genomes.html', my_locals(locals()))


def circos2genomes_main(request):
    return render(request, 'chlamdb/circos2genomes_main.html', my_locals(locals()))


def update_db(server):
    biodb_list = manipulate_biosqldb.get_biodatabase_list(server)
    for biodb in biodb_list:
        update_genomes_db(server, biodb)

def update_genomes_db(server):
    from models import Genome
    from models import Database
    from models import GenDB
    biodb = settings.BIODB

    database = Database.objects.get_or_create(db_name=biodb)[0]
    genome_list = manipulate_biosqldb.get_genome_description_list(server, biodb)
    for genome_description in genome_list:
        genome = Genome.objects.get_or_create(genome_name=genome_description, database=database)[0]
        GenDB.objects.get_or_create(database=database, ref_genome=genome, query_genome=genome, genome_name=genome.genome_name, database_name=database.db_name)



def crossplot(request):

    cache.clear()
    bioentry_in_memory = cache.get('biodb')

    if not bioentry_in_memory:
        cache.set("biodb", {})
    bioentry_in_memory = cache.get("biodb")
    server, db = manipulate_biosqldb.load_db(biodb)

    update_db(server)

    crossplot_form_class = make_crossplot_form("Chlamydia_11_14")
    if request.method == 'POST': 
        form = DBForm(request.POST) #crossplot_form_class(request.POST)

        if form.is_valid():  
            server, db = manipulate_biosqldb.load_db("Chlamydia_11_14")

            update_genomes_db(server, "Chlamydia_11_14")

            reference_genome = str(form.cleaned_data['ref_genome'])
            query_genome = str(form.cleaned_data['query_genome'])
            protein_locus = form.cleaned_data['accession']
            region_size = form.cleaned_data['region_size']

            #description2accession = manipulate_biosqldb.description2accession(server, "saureus1")
            #reference_accession = description2accession[reference_genome]
            #query_accession = description2accession[query_genome]

            orthogroup = manipulate_biosqldb.locus_tag2orthogroup_id(server, protein_locus, "Chlamydia_11_14")
            ortho_detail = list(manipulate_biosqldb.orthogroup_id2locus_tag_list(server, orthogroup, "Chlamydia_11_14"))
            locus_tag_list = []
            for i in range(0, len(ortho_detail)):
                if reference_genome in ortho_detail[i] or query_genome in ortho_detail[i]:
                    locus_tag_list.append(ortho_detail[i][2])
            if "Chlamydia_11_14" not in bioentry_in_memory.keys():
                bioentry_in_memory["Chlamydia_11_14"] = {}

            home_dir = os.path.dirname(__file__)
            temp_location = os.path.join(home_dir, "../assets")
            temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
            name = os.path.basename(temp_file.name)
            bioentry_dict = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, "Chlamydia_11_14", locus_tag_list,
                                                                                  temp_file.name, int(region_size),
                                                                                  bioentry_in_memory["Chlamydia_11_14"])




            envoi = True

    else:  
        form = DBForm() #crossplot_form_class()
        form2 = DBForm()
    return render(request, 'chlamdb/crossplot.html', my_locals(locals()))




def string_page(request, cog_id, genome_accession):
    biodb = settings.BIODB
    from chlamdb.biosqldb import manipulate_biosqldb
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

        cogs_in_chlamdb = []
        cogs_in_reference = []
        cog2description = {}
        for cog in all_cogs:

            sql1 = 'select function, name from COG_cog_names_2014 where COG_id="%s"' % cog
            try:
                data = list(server.adaptor.execute_and_fetchall(sql1,)[0])
                cog2description[cog] = "%s (%s)" % (data[1], data[0])
            except:
                cog2description[cog] = "-"


            try:
                sql = 'select * from COG_locus_tag2gi_hit where COG_id="%s" limit 1;' % (cog)
                sql2 = 'select * from COG_locus_tag2gi_hit where COG_id="%s" and accession="%s" limit 1;' % (cog, 
                                                                                                             genome_accession)

                data = server.adaptor.execute_and_fetchall(sql)
                data2 = server.adaptor.execute_and_fetchall(sql2)
                if len(data)>0:
                    cogs_in_chlamdb.append(cog)
                if len(data2)>0:
                    cogs_in_reference.append(cog)
            except:
                print ('%s not present in %s' % (cog, biodb))
        for i, data in enumerate(string_interactions):
            string_interactions[i] = data + [cog2description[data[0]], cog2description[data[1]]]

        cog_url = '?'
        for i in cogs_in_chlamdb:
            cog_url+= 'cog_list=%s&' % i
        cog_url = cog_url[0:-1]
    except urllib2.URLError:
        connect = False

    return render(request, 'chlamdb/string.html', my_locals(locals()))


def multiple_COGs_heatmap(request):
    biodb = settings.BIODB
    from ete3 import Tree, TextFace
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    cog_list = request.GET.getlist('cog_list')

    cog_filter = '"' + '","'.join(cog_list) + '"'

    cog_annotation_sql = 'select * from COG_cog_names_2014 where COG_id in (%s)' % cog_filter

    cog_annotation = list(server.adaptor.execute_and_fetchall(cog_annotation_sql,))

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    t1 = Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


    sql = 'show columns from comparative_tables_COG' % biodb
    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    ortho_sql = '"' + '","'.join(cog_list) + '"'

    sql = 'select * from comparative_tables_COG where id in (%s)' % (ortho_sql)

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

                    #'first row, print gene names'
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
        asset_path = '/temp/cog_tree.png'
        t1.render(path, dpi=1200)
    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/cog_tree.svg'
        asset_path = '/temp/cog_tree.svg'

        t1.render(path, dpi=800)

    return render(request, 'chlamdb/cog_tree.html', my_locals(locals()))

def pfam_tree(request, orthogroup):
    biodb = settings.BIODB

    task = pfam_tree_task.delay(biodb, 
                                orthogroup)
    print("task", task)
    task_id = task.id

    return render(request, 'chlamdb/pfam_tree.html', my_locals(locals()))

def TM_tree(request, orthogroup):
    biodb = settings.BIODB

    task = TM_tree_task.delay(biodb, 
                              orthogroup)
    print("task", task)
    task_id = task.id

    return render(request, 'chlamdb/TM_tree.html', my_locals(locals()))

def phylogeny(request, orthogroup):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_groups = 'select count(*) from orthology_seqfeature_id2orthogroup t1 ' \
                    ' inner join orthology_orthogroup t2 on t1.orthogroup_id=t2.orthogroup_id where t2.orthogroup_name="%s";' % (orthogroup)

    homologues = server.adaptor.execute_and_fetchall(sql_groups, )[0][0]



    if optional2status["interpro_data"]:
        sql_TM_SP = 'select count(*) from orthology_seqfeature_id2orthogroup t1 ' \
                ' inner join orthology_orthogroup t2 on t1.orthogroup_id=t2.orthogroup_id ' \
                ' inner join interpro_interpro t3 on t1.seqfeature_id=t3.seqfeature_id' \
                ' inner join interpro_signature t4 on t3.signature_id=t4.signature_id ' \
                ' where signature_accession in ("TRANSMEMBRANE", "SIGNAL_PEPTIDE_C_REGION", "SIGNAL_PEPTIDE", "SIGNAL_PEPTIDE_N_REGION") ' \
                ' and t2.orthogroup_name="%s" ; ' % (orthogroup)
                
        tm_count = server.adaptor.execute_and_fetchall(sql_TM_SP, )[0][0]
        
        if tm_count > 0:
            show_tm_tree = True


        task = pfam_tree_task.delay(biodb, 
                                    orthogroup)
    else:
        show_tm_tree = False
        task = basic_tree_task.delay(biodb, 
                                    orthogroup)

    task_id = task.id

    return render(request, 'chlamdb/phylogeny.html', my_locals(locals()))


def refseq_swissprot_tree(request, orthogroup):
    biodb = settings.BIODB

    task = phylogeny_task.delay(biodb, 
                                orthogroup)
    print("refseq_swissprot_tree task", task)
    task_id = task.id

    return render(request, 'chlamdb/best_refseq_swissprot_tree.html', my_locals(locals()))

def multiple_orthogroup_heatmap(request, reference_orthogroup, max_distance=2.2):
    biodb = settings.BIODB
    '''

    multi group heatmap for profiles
    color as a function of profile distance

    :param request:
    :param biodb:
    :param reference_orthogroup:
    :param max_distance:
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    import pandas
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl
    from ete3 import Tree, TextFace

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2' \
               ' on t1.biodatabase_id=t2.biodatabase_id where name="%s";'

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    t1 = Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


    sql = 'select * from interactions_phylo_profiles_eucl_dist' \
          ' where (group_1="%s" or group_2="%s") and euclidian_dist <=%s limit 40;' % (reference_orthogroup,
                                                                                       reference_orthogroup,
                                                                                       max_distance)
    data = list(server.adaptor.execute_and_fetchall(sql,))

    data_frame = pandas.DataFrame(data)
    sorted_data_frame = data_frame.sort(2)

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

    sql = 'show columns from comparative_tables_orthology' % biodb
    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    ortho_sql = '"' + '","'.join(orthogroup2distance.keys()) + '"' + ',"%s"' % reference_orthogroup

    sql = 'select * from comparative_tables_orthology where orthogroup in (%s)' % (biodb, ortho_sql)

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
            if head:

                    # first row, print gene names
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
        asset_path = '/temp/profile_tree_%s.png' % reference_orthogroup
        t1.render(path, dpi=1200)
    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/profile_tree_%s.svg' % reference_orthogroup
        asset_path = '/temp/profile_tree_%s.svg' % reference_orthogroup
        t1.render(path, dpi=800)

    # get data about orthogroups

    match_groups_data, raw_data = biosql_own_sql_tables.orthogroup_list2detailed_annotation(ordered_orthogroups, biodb)


    return render(request, 'chlamdb/profile_tree.html', my_locals(locals()))

def locus2locus(request):
    biodb = settings.BIODB
    NetForm = make_locus2network_form(biodb)

    if request.method == 'POST':

        form = NetForm(request.POST)

        if form.is_valid():

            from chlamdb.biosqldb import manipulate_biosqldb
            from chlamdb.network_d3 import string_networks

            server, db = manipulate_biosqldb.load_db(biodb)


            taxon_id = form.cleaned_data['genome']
            target_list = [i.rstrip() for i in form.cleaned_data['locus_list'].split('\n')]

            filter1 = "'"+"','".join(target_list)+"'"

            sql = 'select locus_tag, orthogroup from orthology_detail where locus_tag in (%s)' % (biodb, filter1)

            data = server.adaptor.execute_and_fetchall(sql,)
            orthogroup2locus_list = {}
            for row in data:
                if row[1] not in orthogroup2locus_list:
                    orthogroup2locus_list[row[1]] = [row[0]]
                else:
                    orthogroup2locus_list[row[1]].append(row[0])

            orthogroup_list = orthogroup2locus_list.keys()

            filter2 = "'"+"','".join(orthogroup_list)+"'"

            sql = 'select locus_tag, orthogroup from orthology_detail where orthogroup in (%s) and taxon_id=%s' % (filter2,
                                                                                                                   taxon_id)
            data = server.adaptor.execute_and_fetchall(sql,)
            orthogroup2locus_list_corresp = {}
            for row in data:
                if row[1] not in orthogroup2locus_list_corresp:
                    orthogroup2locus_list_corresp[row[1]] = [row[0]]
                else:
                    orthogroup2locus_list_corresp[row[1]].append(row[0])
            for group in orthogroup_list:
                if group not in orthogroup2locus_list_corresp:
                    orthogroup2locus_list_corresp[group] = ['-']
    else:
        form = NetForm()
    return render(request, 'chlamdb/locus2locus.html', my_locals(locals()))

def interactions_genome(request):
    biodb = settings.BIODB
    NetForm = make_locus2network_form(biodb)

    if request.method == 'POST':
        form = NetForm(request.POST)

        if form.is_valid():

            from chlamdb.biosqldb import manipulate_biosqldb
            from chlamdb.network_d3 import string_networks

            server, db = manipulate_biosqldb.load_db(biodb)


            taxon_id = form.cleaned_data['genome']
            target_list = [i.rstrip() for i in form.cleaned_data['locus_list'].split('\n')]

            sql = 'select seqfeature_id from annotation_seqfeature_id2locus where taxon_id="%s"' % (taxon_id)

            seqfeature_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            #all_groups_neig = string_networks.find_links_recusrsive(biodb, locus_list, 0.8, n_comp_cutoff=2)

            sql = 'select seqfeature_id, gene, product from orthology_detail where taxon_id="%s"' % (taxon_id)
            locus2gene_product = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select seqfeature_id,interpro_accession,interpro_description from interpro ' \
                  ' where taxon_id=%s and interpro_accession != "0" group by locus_tag,interpro_accession' % (taxon_id)

            locus2interpro = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            if len (target_list) > 1:
                locus_filter = string_networks.get_subgraph(biodb, seqfeature_id_list, 0.8, target_list)
            else:
                locus_filter = seqfeature_id_list
            script = string_networks.generate_network(biodb,
                                                      locus_filter,
                                                      target_list,
                                                      0.8,
                                                      scale_link=True,
                                                      width=200,
                                                      height=200,
                                                      interpro=locus2interpro,
                                                      annot=locus2gene_product)


    else:
        form = NetForm()
    return render(request, 'chlamdb/interactions_genome.html', my_locals(locals()))




def interactions_genome_string(request):
    biodb = settings.BIODB
    NetForm = make_locus2network_form(biodb)

    if request.method == 'POST':
        form = NetForm(request.POST)

        if form.is_valid():

            from chlamdb.biosqldb import manipulate_biosqldb
            from chlamdb.network_d3 import string_networks

            server, db = manipulate_biosqldb.load_db(biodb)


            taxon_id = form.cleaned_data['genome']
            target_list = [i.rstrip() for i in form.cleaned_data['locus_list'].split('\n')]

            sql = 'select locus_tag from orthology_detail where taxon_id="%s"' % (taxon_id)

            locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            #all_groups_neig = string_networks.find_links_recusrsive(biodb, locus_list, 0.8, n_comp_cutoff=2)

            sql = 'select locus_tag, gene, product from orthology_detail where taxon_id="%s"' % (taxon_id)
            locus2gene_product = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            sql = 'select locus_tag,interpro_accession,interpro_description from interpro ' \
                  ' where taxon_id=%s and interpro_accession != "0" group by locus_tag,interpro_accession' % (taxon_id)

            locus2interpro = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            if len (target_list) > 1:
                locus_filter = string_networks.get_subgraph(biodb, locus_list, 0.8, target_list)
            else:
                locus_filter = locus_list
            script = string_networks.generate_network_string(biodb,
                                                      locus_filter,
                                                      target_list,
                                                      ratio_limit=0.7,
                                                      taxon_id=taxon_id,
                                                      width=200,
                                                      height=200,
                                                      interpro=locus2interpro,
                                                      annot=locus2gene_product)


    else:
        form = NetForm()
    return render(request, 'chlamdb/interactions_genome_string.html', my_locals(locals()))




def interactions(request, locus_tag):
    biodb = settings.BIODB
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.network_d3 import string_networks

    server, db = manipulate_biosqldb.load_db(biodb)

    print ('get interactors -- %s -- %s' % (biodb, locus_tag))

    sql = 'select orthogroup from orthology_detail where locus_tag="%s"' % (locus_tag)

    orthogroup = server.adaptor.execute_and_fetchall(sql,)[0][0]

    all_groups_profile_jac, cutoff_jac, too_much_hits_jac = string_networks.successive_cutof_search(biodb,
                                                                                                    "jac",
                                                                                                    orthogroup,
                                                                                                    0.15,
                                                                                                    0.1,
                                                                                                    0.05,
                                                                                                    0)

    all_groups_profile_eucl, cutoff_eucl, too_much_hits_eucl = string_networks.successive_cutof_search(biodb,
                                                                                                        "eucl",
                                                                                                        orthogroup,
                                                                                                        2.2,
                                                                                                        2,
                                                                                                        1,
                                                                                                        0)

    if len(all_groups_profile_eucl) > 1:
        profile_match_eucl = True
        profile_match = True

    if len(all_groups_profile_jac) > 1:
        profile_match_jac = True
        profile_match = True

    if biodb != 'chlamydia_04_16':
        sql = 'select seqfeature_id from custom_tables_locus2seqfeature_id where locus_tag="%s"' % (locus_tag)

        seqfeature_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

        all_groups_neig = string_networks.find_links_recusrsive(biodb, [seqfeature_id], 0.8, n_comp_cutoff=0)
    else:
        all_groups_neig = string_networks.find_links_recusrsive(biodb, [locus_tag], 0.8, n_comp_cutoff=0)
    if len(all_groups_neig) == 0:
        neig_match = False
    else:
        neig_match = True

    return render(request, 'chlamdb/interactions.html', my_locals(locals()))


def plot_heatmap(request, type):
    import seaborn as sns
    from datetime import datetime

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)
    form_class = make_venn_from(db)

    if request.method != "POST":
        form_venn = form_class()
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    form_venn = form_class(request.POST)
    if not form_venn.is_valid():
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    taxon_ids = form_venn.get_taxids()
    if type=="COG":
        mat = db.get_cog_hits(taxon_ids, indexing="taxid", search_on="taxid")
    elif type=="orthology":
        mat = db.get_og_count(taxon_ids)
    elif type == "ko":
        mat = db.get_ko_hits(taxon_ids)
    else:
        form_venn = form_class()
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    target2description = db.get_genomes_description().description.to_dict()
    mat.columns = (target2description[i] for i in mat.columns.values)
    cur_time = datetime.now().strftime("%H%M%S")

    filename = f"heatmap_{cur_time}.png"
    path = settings.BASE_DIR + '/assets/temp/' + filename
    asset_path = '/temp/' + filename
    cm = sns.clustermap(mat)
    cm.savefig(path)
    envoi_heatmap = True
    return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))


def profile_interactions(request, orthogroup, distance):
    biodb = settings.BIODB
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.network_d3 import string_networks
    from chlamdb.biosqldb import biosql_own_sql_tables



    server, db = manipulate_biosqldb.load_db(biodb)

    if distance == "eucl":
        all_groups_profile, cutoff, too_much_hits = string_networks.successive_cutof_search(biodb,
                                                                                            "eucl",
                                                                                            orthogroup,
                                                                                            2.2,
                                                                                            2,
                                                                                            1,
                                                                                            0)
    if distance == "jac":
        all_groups_profile, cutoff, too_much_hits = string_networks.successive_cutof_search(biodb,
                                                                                            "jac",
                                                                                            orthogroup,
                                                                                            0.15,
                                                                                            0.1,
                                                                                            0.05,
                                                                                            0)                                                                   

    if len(all_groups_profile) <= 1:
        profile_match = False
    else:
        profile_match = True
        from chlamdb.phylo_tree_display import ete_motifs
        match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(all_groups_profile, biodb)
        match = True

        script = string_networks.generate_network_profile(biodb, 
                                                          all_groups_profile, 
                                                          [orthogroup], 
                                                          euclidian_distance_limit=cutoff, 
                                                          scale_link=True,
                                                          interpro=False,
                                                          annot=False,
                                                          distance=distance)

        taxon2orthogroup2count_all = ete_motifs.get_taxon2orthogroup2count(biodb, all_groups_profile)

        labels = all_groups_profile
        orthogroup_n = all_groups_profile.index(orthogroup)
        tree, style = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2orthogroup2count_all,reference_column=orthogroup_n)

        path = settings.BASE_DIR + '/assets/temp/ortho_tree.svg'
        asset_path = '/temp/ortho_tree.svg'
        tree.render(path, dpi=500, tree_style=style)


    return render(request, 'chlamdb/profile_interactions.html', my_locals(locals()))

def neig_interactions(request, locus_tag):
    biodb = settings.BIODB
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.network_d3 import string_networks
    from chlamdb.biosqldb import biosql_own_sql_tables




    server, db = manipulate_biosqldb.load_db(biodb)

    if biodb != 'chlamydia_04_16':
        sql = 'select seqfeature_id from custom_tables_locus2seqfeature_id where locus_tag="%s"' % (locus_tag)

        seqfeature_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
        locus_tag_list = string_networks.find_links_recusrsive(biodb, [seqfeature_id], 0.8, n_comp_cutoff=1)
    else:
        locus_tag_list = string_networks.find_links_recusrsive(biodb, [locus_tag], 0.8, n_comp_cutoff=1)
    if len(locus_tag_list) == 0:
        match = False
    else:
        from chlamdb.phylo_tree_display import ete_motifs

        if biodb != 'chlamydia_04_16':
            locus_tag_list = [str(i) for i in locus_tag_list]
            filter_seqfeatures = ','.join(locus_tag_list)
            sql = 'select locus_tag from custom_tables_locus2seqfeature_id where seqfeature_id in (%s)' % (filter_seqfeatures)
            locus_tag_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

        #match_groups_data, extract_result = biosql_own_sql_tables.orthogroup_list2detailed_annotation(all_groups, biodb)
        locus2annot, \
        locus_tag2cog_catego, \
        locus_tag2cog_name, \
        locus_tag2ko, \
        pathway2category, \
        module2category, \
        ko2ko_pathways, \
        ko2ko_modules,\
        locus2interpro = get_locus_annotations(biodb, locus_tag_list)

        locus2start = {}
        # get reference orthogroup
        for i in locus2annot:
            locus2start[i[2]] = int(i[4])
            if i[2] == locus_tag:
                orthogroup = i[1]

        middle_position = sorted(locus2start.values())[int(len(locus2start.values())/2)]
        middle_locus_tag = list(locus2start.keys())[list(locus2start.values()).index(middle_position)]

        # get complete orthogroup list
        all_groups = list(set([i[1] for i in locus2annot]))

        # plot phylogenetic profile
        taxon2orthogroup2count_all = ete_motifs.get_taxon2orthogroup2count(biodb, all_groups)
        labels = all_groups
        orthogroup_n = all_groups.index(orthogroup)

        tree, style = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2orthogroup2count_all, reference_column=orthogroup_n)
        path = settings.BASE_DIR + '/assets/temp/ortho_tree.svg'
        asset_path = '/temp/ortho_tree.svg'
        tree.render(path, dpi=500, tree_style=style)

        sql = 'select taxon_id from orthology_detail where locus_tag ="%s" group by taxon_id' % (locus_tag)

        taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]

        plot_url = "?t=%s" % taxon_list[0] +('&t=').join((taxon_list[1:]))

        match = True

        home_dir = os.path.dirname(__file__)
        temp_location = os.path.join(home_dir, "../assets")
        temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
        name = os.path.basename(temp_file.name)
        name_png = name.split('.')[0] + '.png'
        locus_tags, orthogroup_list = mysqldb_plot_genomic_feature.proteins_id2cossplot(server, db, biodb, [middle_locus_tag],
                                                                          temp_file.name, int(29000),
                                                                          cache, color_locus_list=locus_tag_list)

        sql = 'select seqfeature_id from annotation_seqfeature_id2locus where locus_tag in ("%s")' % ('","'.join(locus_tag_list))
        seqfeature_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    script = string_networks.generate_network(biodb, seqfeature_id_list, [locus_tag], 0.7, scale_link=True)

    return render(request, 'chlamdb/neig_interactions.html', my_locals(locals()))


def similarity_network(request, orthogroup, annotation):
    biodb = settings.BIODB

    server, db = manipulate_biosqldb.load_db(biodb)

    if annotation == 'ko':
        sql = 'select locus_tag,t5.ko_accession from annotation_seqfeature_id2locus t1 ' \
              ' inner join orthology_seqfeature_id2orthogroup t2 on t1.seqfeature_id=t2.seqfeature_id ' \
              ' inner join orthology_orthogroup t3 on t2.orthogroup_id=t3.orthogroup_id ' \
              ' inner join enzyme_seqfeature_id2ko t4 on t1.seqfeature_id=t4.seqfeature_id ' \
              ' inner join enzyme_ko_annotation t5 on t4.ko_id=t5.ko_id where t3.orthogroup_name="%s";' % (orthogroup)

    locus2annotation = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    unique_annotations = list(set(locus2annotation.values()))
    annotation2group = {}
    unique_annotations.append(None)
    for n, value in enumerate(unique_annotations):
        annotation2group[value] = n
    
    sql = 'select orthogroup, locus_a, locus_b, identity from orthology_identity where identity!=0 and orthogroup="%s";' % (orthogroup)
    
    data = [i for i in server.adaptor.execute_and_fetchall(sql,)]

    edge_list = []
    node2id = {}
    id = 0
    node_list = []
    for row in data:
        node_1 = row[1]
        node_2 = row[2]
        identity = row[3]
        
        if node_1 == node_2:
            continue

        if node_1 not in locus2annotation:
            locus2annotation[node_1] = None
        if node_2 not in locus2annotation:
            locus2annotation[node_2] = None
        
        if node_1 not in node2id:
            node2id[node_1] = id
            id+=1
            node_list.append({"id": node2id[node_1], "label": node_1, "title": 'Locus_tag: ' + node_1 + '<br> %s: %s ' % (annotation, locus2annotation[node_1]) , "group": annotation2group[locus2annotation[node_1]]})
        if node_2 not in node2id:
            node2id[node_2] = id
            id+=1
            node_list.append({"id": node2id[node_2], "label": node_2, "title": 'Locus_tag: ' + node_2 + '<br> %s: %s ' % (annotation, locus2annotation[node_2]), "group": annotation2group[locus2annotation[node_2]]})
        
        
        if float(identity) >= 20:
            edge_list.append({"from": node2id[node_1], "to": node2id[node_2], "length": 100-float(identity), "label": identity})        
    
    '''
    [{from: 1, to: 15},
    {from: 1, to: 97} ]
    
      {id: 735, label: 'Yuya Osako', title: 'Country: ' + 'Japan' + '<br>' + 'Team: ' + '1860 München', value: 22, group: 27, x: 806.69904, y: 633.54565},
      {id: 736, label: 'Zvjezdan Misimovic', title: 'Country: ' + 'Bosnia and Herzegovina' + '<br>' + 'Team: ' + 'Guizhou Renhe', value: 22, group: 20, x: 1277.4697, y: -479.12265}
  
    '''
    envoi = True

    return render(request, 'chlamdb/similarity_network.html', my_locals(locals()))



def orthogroup_conservation_tree_legacy(request, orthogroup_or_locus):
    biodb = settings.BIODB
    '''

    produit un profile presence/absence pour un orthogroup donne
    si locus_tag, ajouter une colonne avec l'identité du groupe le plus proche

    :param request:
    :param biodb:
    :param orthogroup:
    :return:
    '''


    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.phylo_tree_display import ete_heatmap_conservation
    from chlamdb.biosqldb import shell_command
    from chlamdb.network_d3 import string_networks

    server, db = manipulate_biosqldb.load_db(biodb)

    if orthogroup_or_locus.startswith("group_"):
        input_type = 'orthogroup'
        orthogroup = orthogroup_or_locus
        taxon2identity_closest = False
        taxon2locus_tag_closest = False
        taxon_id = False
    else:
        input_type = 'locus_tag'
        sql = 'select orthogroup, taxon_id from orthology_detail where locus_tag="%s"' % (orthogroup_or_locus)
        print(sql)
        data = server.adaptor.execute_and_fetchall(sql, )[0]
        orthogroup = data[0]
        taxon_id = data[1]

        sql2 = 'select taxon_2,B.locus_tag,identity from (select * from custom_tables_locus2seqfeature_id t1 ' \
               ' inner join comparative_tables_identity_closest_homolog2 t2 on t1.seqfeature_id=t2.locus_1 ' \
               ' where locus_tag="%s") A inner join custom_tables_locus2seqfeature_id B on A.locus_2=B.seqfeature_id;' % (orthogroup_or_locus)

        identity_data = server.adaptor.execute_and_fetchall(sql2, )
        taxon2identity_closest = {}
        taxon2locus_tag_closest = {}

        for row in identity_data:
            taxon2identity_closest[str(row[0])] = row[2]
            taxon2locus_tag_closest[str(row[0])] = row[1]
        taxon2locus_tag_closest[str(taxon_id)] = orthogroup_or_locus
        taxon2identity_closest[str(taxon_id)] = 100

    try: # TODO deal with missing interactions tables
        all_groups_profile_eucl, cutoff_eucl, too_much_hits_eucl = string_networks.successive_cutof_search(biodb,
                                                                                                            "eucl",
                                                                                                            orthogroup,
                                                                                                            2.2,
                                                                                                            2,
                                                                                                            1,
                                                                                                            0)
    except:
        all_groups_profile_eucl, cutoff_eucl, too_much_hits_eucl = [], None, None
    try: # TODO deal with missing interactions tables
        all_groups_profile_jac, cutoff_jac, too_much_hits_jac = string_networks.successive_cutof_search(biodb,
                                                                                                        "jac",
                                                                                                        orthogroup,
                                                                                                        0.15,
                                                                                                        0.1,
                                                                                                        0.05,
                                                                                                        0)
    except:
        all_groups_profile_jac, cutoff_jac, too_much_hits_jac = [], None, None
        
    profile_match_jac = False
    profile_match_eucl = False
    if len(all_groups_profile_jac) > 1:
        profile_match_jac = True
    if len(all_groups_profile_eucl) > 1:
        profile_match_eucl = True
    
    if input_type != 'orthogroup':
        locus_list = list(taxon2locus_tag_closest.values())

    asset_path = '/temp/phylo.svg'
    path = settings.BASE_DIR + '/assets/' + asset_path

    a,b,c = shell_command.shell_command("rm %s" % path)


    sql_grp = 'select taxon_id,count(*) from  orthology_detail where orthogroup="%s" ' \
              ' group by taxon_id;' % (orthogroup)

    taxid2n = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_grp,))
    tree_sql = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 ' \
               ' on t1.biodatabase_id=t2.biodatabase_id where t2.name="%s"' % biodb
    tree = server.adaptor.execute_and_fetchall(tree_sql,)[0][0]

    taxon_profile = []
    for i in taxid2n:
        if taxid2n[i]>0:
            taxon_profile.append(i)

    url_pattern = ''
    for i in taxon_profile:
        url_pattern += 'taxons_profile=%s&' % i
    url_pattern=url_pattern[0:-1]

    #t1, leaf_number = ete_heatmap_conservation.plot_heat_tree(biodb, taxid2n, tree)
    t1, leaf_number, tree_style = ete_heatmap_conservation.plot_heatmap_tree_locus(biodb,
                                                                                    tree,
                                                                                    taxid2n,
                                                                                    taxid2identity= taxon2identity_closest,
                                                                                    taxid2locus = taxon2locus_tag_closest,
                                                                                    reference_taxon=taxon_id,
                                                                                    n_paralogs_barplot=True)
    shell_command.shell_command('rm %s' % path)

    t1.render(path, 
              tree_style=tree_style,
              dpi=300, 
              w=800)

    return render(request, 'chlamdb/orthogroup_conservation.html', my_locals(locals()))


def priam_kegg(request):
    biodb = settings.BIODB
    priam_form_class = make_priam_form(biodb)

    if request.method == 'POST': 
        form = priam_form_class(request.POST)
        if form.is_valid():

            server, db = manipulate_biosqldb.load_db(biodb)

            genome = form.cleaned_data['genome']

            sql = 'select taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
                  ' where t2.name="%s" and t1.accession="%s"' % (biodb, 
                                                                 genome)
            taxon_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

            sql = 'select pathway_name,pathway_category,description from enzyme_locus2ko t1 ' \
                  ' inner join enzyme_pathway2ko_v1 t2 on t1.ko_id=t2.ko_id  ' \
                  ' inner join enzyme_kegg_pathway t3 on t2.pathway_id=t3.pathway_id ' \
                  ' where taxon_id=%s and pathway_category ' \
                  ' not in ("6.9 Infectious diseases: Viral", "6.8 Infectious diseases: Bacterial", ' \
                  ' "6.7 Endocrine and metabolic diseases", "6.5 Substance dependence", ' \
                  ' "6.4 Neurodegenerative diseases", "6.3 Immune diseases", "6.2 Cancers: Specific types", ' \
                  ' "6.12 Drug resistance: Antineoplastic", "6.1 Cancers: Overview", "5.9 Aging", ' \
                  ' "5.6 Nervous system", "5.7 Sensory system", "5.8 Development", "6.10 Infectious diseases: Parasitic", ' \
                  ' "5.1 Immune system", "5.3 Circulatory system", "5.2 Endocrine system", ' \
                  ' "5.5 Excretory system", "5.10 Environmental adaptation","4.4 Cellular community - eukaryotes", ' \
                  ' "5.4 Digestive system", "3.3 Signaling molecules and interaction", "1.0 Global and overview maps") ' \
                  ' group by pathway_name order by pathway_category;' % (taxon_id)
            data = server.adaptor.execute_and_fetchall(sql,)


            envoi = True

    else:  
        form = priam_form_class()

    return render(request, 'chlamdb/priam_kegg.html', my_locals(locals()))



def locus_list2circos(request, target_taxon):
    biodb = settings.BIODB
    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    locus_list = [i for i in request.GET.getlist('l')]

    # 1. check if locus_list is part of the target genome
    # 2. if not part of the target genome, get orthogroup id

    locus_filter = '"'+'","'.join(locus_list)+'"'
    sql = 'select locus_tag, product from orthology_detail where taxon_id=%s and locus_tag in (%s)' % (target_taxon,
                                                                                                       locus_filter)
    data = server.adaptor.execute_and_fetchall(sql,)
    locus_target_genome = [i[0] for i in data]

    locus2label = {}
    for row in data:
        locus2label[row[0]] = row[0] #+ "-%s" % re.sub(" ","-",row[1])

    locus_other_genomes = []
    for locus in locus_list:
        if locus not in locus_target_genome:
            locus_other_genomes.append(locus)
    locus_other_genomes_filter = '"' + '","'.join(locus_other_genomes) + '"'
    sql2 = 'select B.locus_tag, A.orthogroup, B.product from (select orthogroup from orthology_detail where locus_tag in (%s) group by orthogroup) A' \
           ' inner join orthology_detail B on A.orthogroup=B.orthogroup where taxon_id=%s' % (locus_other_genomes_filter,
                                                                                              target_taxon)

    data = server.adaptor.execute_and_fetchall(sql2,)

    for row in data:
        locus2label[row[0]] = row[1] #+ "-%s" % re.sub(" ","-",row[2])

    reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, target_taxon, biodb) # ["NC_009648"] NZ_CP009208 NC_016845

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

    draft_data = []
    for biorecord in record_list:
        draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

    home_dir = os.path.dirname(__file__)

    temp_location = os.path.join(home_dir, "../assets/circos/")

    myplot = circos.CircosAccession2multiplot(server,
                                              db,
                                              biodb,
                                              record_list,
                                              [],
                                              locus_highlight=locus2label.keys(),
                                              out_directory=temp_location,
                                              draft_fasta=draft_data,
                                              href="/chlamdb/locusx/",
                                              ordered_taxons=[],
                                              locus2label=locus2label,
                                              show_homologs=False,
                                              radius=0.3)

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

    original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
    map_file = "circos/%s.html" % ref_name
    svg_file = "circos/%s.svg" % ref_name
    map_name = ref_name

    envoi_circos = True


    return render(request, 'chlamdb/locus2circos.html', my_locals(locals()))


def hmm2circos(request):
    from chlamdb.phylo_tree_display import ete_motifs
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    hmm_form = hmm_sets_form_circos(biodb)

    if request.method == 'POST': 
        form = hmm_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import hmm_heatmap
            from chlamdb.plots import gbk2circos
            from chlamdb.biosqldb import circos

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            hmm_set = form.cleaned_data['hmm_set']
            reference_taxon = form.cleaned_data['genome']
            score_cutoff = form.cleaned_data['score_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']


            sql = 'select locus_tag from hmm.hmm_sets t1 ' \
                  ' inner join hmm.hmm_sets_entry t2 on t1.set_id=t2.set_id ' \
                  ' inner join hmm_hmm_hits_annotated_genome t3 on t2.hmm_id=t3.hmm_id' \
                  ' inner join custom_tables_locus2seqfeature_id t4 on t3.seqfeature_id=t4.seqfeature_id ' \
                  ' where t1.name="%s" and t3.taxon_id=%s and bitscore>=%s ' \
                  ' and query_coverage>=%s order by bitscore;' % (hmm_set,
                                                                  reference_taxon,
                                                                  score_cutoff,
                                                                  query_coverage_cutoff)

            target_locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


            sql2 = 'select locus_tag,t5.name  from hmm.hmm_sets t1 ' \
                  ' inner join hmm.hmm_sets_entry t2 on t1.set_id=t2.set_id ' \
                  ' inner join hmm_hmm_hits_annotated_genome t3 on t2.hmm_id=t3.hmm_id' \
                  ' inner join custom_tables_locus2seqfeature_id t4 on t3.seqfeature_id=t4.seqfeature_id ' \
                  ' inner join hmm.hmm_profiles t5 on t2.hmm_id=t5.hmm_id' \
                  ' where t1.name="%s" and t3.taxon_id=%s and bitscore>=%s ' \
                   ' and query_coverage>=%s order by bitscore;' % (hmm_set, 
                                                                    reference_taxon, 
                                                                    score_cutoff, 
                                                                    query_coverage_cutoff)

            locus2label = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

            description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

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

            myplot = circos.CircosAccession2multiplot(server,
                                                      db,
                                                      biodb,
                                                      record_list,
                                                      [],
                                                      locus_highlight=target_locus_list,
                                                      out_directory=temp_location,
                                                      draft_fasta=draft_data,
                                                      href="/chlamdb/locusx/",
                                                      ordered_taxons = ordered_taxons,
                                                      locus2label=locus2label,
                                                      show_homologs=False)

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

            original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
            map_file = "circos/%s.html" % ref_name
            svg_file = "circos/%s.svg" % ref_name
            map_name = ref_name

            envoi_circos = True

    else:  
        form = hmm_form()

    return render(request, 'chlamdb/hmm2circos.html', my_locals(locals()))




def transporters_list(request):
    from chlamdb.phylo_tree_display import ete_motifs
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    transporters_form = transporters_superfam_form(biodb, True)

    if request.method == 'POST': 
        form = transporters_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import blast_heatmap
            #if request.method == 'POST': 
            genome = form.cleaned_data['genome']
            transporter_superfamily = form.cleaned_data['transporter_superfamily'][0]
            score_cutoff = form.cleaned_data['score_cutoff']
            evalue_cutoff = form.cleaned_data['evalue_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']
            hit_coverage_cutoff = form.cleaned_data['hit_coverage_cutoff']

            if transporter_superfamily == 'all':
                sql = 'select t8.locus_tag,t3.description,t1.n_hsps, t1.evalue, t1.bitscore_first_hsp, ' \
                      ' t1.identity, t1.query_TMS, t1.hit_TMS, t1.query_cov, t1.hit_cov,t7.uniprot_accession, ' \
                      ' t7.substrate, t7.uniprot_description, t5.description,t6.description, t3.tc_name from transporters_transporters t1 ' \
                      ' inner join transporters.transporter_table t2 on t1.transporter_id=t2.transporter_id ' \
                      ' inner join transporters.tc_table t3 on t2.family=t3.tc_id ' \
                      ' inner join transporters.tc_table t4 on t2.superfamily=t4.tc_id ' \
                      ' inner join transporters.tc_table t5 on t2.subfamily=t5.tc_id ' \
                      ' inner join transporters.tc_table t6 on t2.transporter_id=t6.tc_id ' \
                      ' inner join transporters.uniprot_table t7 on t1.hit_uniprot_id=t7.uniprot_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t8 on t1.seqfeature_id=t8.seqfeature_id ' \
                      ' where query_cov>=%s and hit_cov>=%s and evalue<=%s and bitscore_first_hsp>=%s ' \
                      ' and t1.taxon_id=%s;' % (query_coverage_cutoff,
                                                hit_coverage_cutoff,
                                                evalue_cutoff,
                                                score_cutoff,
                                                genome)

                data = list(server.adaptor.execute_and_fetchall(sql,))
                for n, row in enumerate(data):
                    data[n] = list(data[n])
                    data[n][-4] = ','.join(set([i.rstrip().lstrip() for i in data[n][-4].split(',')]))
                    for i in range(0,len(row)):
                        data[n][i] = str(data[n][i])#.decode("latin-1")
                envoi = True

            else:
                sql = 'select t8.locus_tag,t3.description,t1.n_hsps, t1.evalue, t1.bitscore_first_hsp, ' \
                      ' t1.identity, t1.query_TMS, t1.hit_TMS, t1.query_cov, t1.hit_cov,t7.uniprot_accession, ' \
                      ' t7.substrate, t7.uniprot_description, t5.description,t6.description, t3.tc_name from transporters_transporters t1 ' \
                      ' inner join transporters.transporter_table t2 on t1.transporter_id=t2.transporter_id ' \
                      ' inner join transporters.tc_table t3 on t2.family=t3.tc_id ' \
                      ' inner join transporters.tc_table t4 on t2.superfamily=t4.tc_id ' \
                      ' inner join transporters.tc_table t5 on t2.subfamily=t5.tc_id ' \
                      ' inner join transporters.tc_table t6 on t2.transporter_id=t6.tc_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t6 on t1.seqfeature_id=t6.seqfeature_id ' \
                      ' inner join transporters.uniprot_table t7 on t1.hit_uniprot_id=t7.uniprot_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t8 on t1.seqfeature_id=t8.seqfeature_id ' \
                      ' where query_cov>=%s and hit_cov>=%s and evalue<=%s and bitscore_first_hsp>=%s ' \
                      ' and t4.description="%s" ' \
                      ' and t1.taxon_id=%s;' % (query_coverage_cutoff,
                                                hit_coverage_cutoff,
                                                evalue_cutoff,
                                                score_cutoff,
                                                transporter_superfamily,
                                                genome)

                data = list(server.adaptor.execute_and_fetchall(sql,))
                for n, row in enumerate(data):
                    data[n] = list(data[n])
                    for i in range(0,len(row)):
                        data[n][i] = str(data[n][i]) #.decode("latin-1")

                envoi = True


    else:  
        form = transporters_form()

    return render(request, 'chlamdb/transporters_table.html', my_locals(locals()))


def transporters_family(request, family):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs
    from chlamdb.plots import transporters_heatmap

    server, db = manipulate_biosqldb.load_db(biodb)

    from chlamdb.plots import blast_heatmap
    #if request.method == 'POST': 


    from ete3 import Tree
    sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
               ' where t2.name="%s";' % biodb
    server, db = manipulate_biosqldb.load_db(biodb)
    tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]

    acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
    filter = '"' + '","'.join(acc_list) + '"'
    sql = 'select taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
          ' where t2.name="%s" and accession in (%s)' % (biodb, filter)
    
    taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    try:
        t2 = t1.prune(taxon_list)
    except:
        pass
    t1.ladderize()

    evalue_cutoff = 0.05
    score_cutoff = 10
    query_coverage_cutoff = 0.5
    hit_coverage_cutoff = 0.5

    transporter_family_list = [ "2.A.15",
                                "2.A.111",
                                "2.A.12",
                                "2.A.17",
                                "2.A.19",
                                "2.A.2",
                                "2.A.20",
                                "2.A.21",
                                "2.A.22",
                                "2.A.23",
                                "2.A.25",
                                "2.A.26",
                                "2.A.27",
                                "2.A.28",
                                "2.A.3",
                                "2.A.30",
                                "2.A.33",
                                "2.A.35",
                                "2.A.36",
                                "2.A.37",
                                "2.A.39",
                                "2.A.4",
                                "2.A.40",
                                "2.A.41",
                                "2.A.47",
                                "2.A.55",
                                "2.A.58",
                                "2.A.6",
                                "2.A.62",
                                "2.A.63",
                                "2.A.66",
                                "2.A.7",
                                "2.A.83"
                                ]


    superfam_list, taxon2code2count = transporters_heatmap.transporter_family_heatmap(biodb,
                                                                                      [family],
                                                                                      evalue_cutoff,
                                                                                      score_cutoff,
                                                                                      query_coverage_cutoff,
                                                                                      hit_coverage_cutoff)

    #family_filter = '"','","'.join([family]),'"'
    family_filter = '"%s"' % family
    sql = 'select locus_tag,orthogroup_name, A.taxon_id from (' \
          ' select t1.taxon_id,t3.description,t1.seqfeature_id from transporters_transporters t1  ' \
          ' inner join transporters.transporter_table t2 on t1.transporter_id=t2.transporter_id ' \
          ' inner join transporters.tc_table t3 on t2.family=t3.tc_id ' \
          ' inner join transporters.tc_table t4 on t2.family=t4.tc_id ' \
          ' where query_cov>=%s and hit_cov>=%s and evalue<=%s and bitscore_first_hsp>=%s and t4.tc_name in (%s))A ' \
          ' inner join annotation_seqfeature_id2locus B on A.seqfeature_id=B.seqfeature_id' \
          ' inner join orthology_seqfeature_id2orthogroup C on A.seqfeature_id=C.seqfeature_id' \
          ' inner join orthology_orthogroup D on C.orthogroup_id=D.orthogroup_id;' % (query_coverage_cutoff,
                                                                                      hit_coverage_cutoff,
                                                                                      score_cutoff,
                                                                                      evalue_cutoff,
                                                                                      family_filter)

    taxon2orthogroup2transporter_family = {}
    data = server.adaptor.execute_and_fetchall(sql,)
    locus_list = []
    orthogroup_list = []
    for i in data:
        locus = i[0]
        orthogroup = i[1]
        taxon_id = int(i[2])
        locus_list.append(locus)
        if orthogroup not in orthogroup_list:
            orthogroup_list.append(orthogroup)
        if taxon_id not in taxon2orthogroup2transporter_family:
            taxon2orthogroup2transporter_family[taxon_id] = {}
            taxon2orthogroup2transporter_family[taxon_id][orthogroup] = superfam_list
        else:
            if orthogroup not in taxon2orthogroup2transporter_family[taxon_id]:
                taxon2orthogroup2transporter_family[taxon_id][orthogroup] = superfam_list
            else:
                pass

    taxon2orthogroup2count = ete_motifs.get_taxon2orthogroup2count(biodb, 
                                                                   orthogroup_list)

    merged_dico = taxon2orthogroup2count
    for taxon in taxon2code2count:
        merged_dico[str(taxon)] = taxon2code2count[taxon]

    labels = superfam_list + orthogroup_list

    tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                labels,
                                                merged_dico,
                                                taxon2group2value=taxon2orthogroup2transporter_family,
                                                highlight_first_column=True)


    filter = '"'+'","'.join(locus_list)+'"'
    if db_driver == 'mysql':
        sql = 'select locus_tag, accession, start, stop, gene, product, n_genomes, orthogroup, ' \
            ' CHAR_LENGTH(translation) from orthology_detail ' \
            ' where locus_tag in (%s)' % (filter)
    if db_driver == 'sqlite':
        sql = 'select locus_tag, accession, start, stop, gene, product, n_genomes, orthogroup, ' \
           ' LENGTH(translation) from orthology_detail ' \
           ' where locus_tag in (%s)' % (filter)       
    locus_annot = [list(i) for i in server.adaptor.execute_and_fetchall(sql,)]

    sql = 'select locus_tag,t3.COG_name,t5.code,t3.description ' \
          ' from COG_seqfeature_id2best_COG_hit t1 ' \
          ' inner join annotation_seqfeature_id2locus t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' inner join COG_cog_names_2014 t3 on t1.hit_cog_id=t3.cog_id ' \
          ' inner join COG_cog_id2cog_category t4 on t3.cog_id=t4.cog_id ' \
          ' inner join COG_code2category t5 on t4.category_id=t5.category_id where locus_tag in (%s);' % (filter)

    locus2COG_data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    for one_locus_annot in locus_annot:
        try:
            one_locus_annot.append(locus2COG_data[one_locus_annot[0]][0])
            one_locus_annot.append(locus2COG_data[one_locus_annot[0]][1])
            one_locus_annot.append(locus2COG_data[one_locus_annot[0]][2])
        except:
            one_locus_annot.append('-')
            one_locus_annot.append('-')
            one_locus_annot.append('-')

    locus2annot, \
    locus_tag2cog_catego, \
    locus_tag2cog_name, \
    locus_tag2ko, \
    pathway2category, \
    module2category, \
    ko2ko_pathways, \
    ko2ko_modules,\
    locus2interpro = get_locus_annotations(biodb, locus_list)

    '''
    tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                superfam_list,
                                                taxon2code2count,
                                                tree=t1,
                                                identity_scale=False,
                                                show_labels=True,
                                                column_scale=True,
                                                as_float=False,
                                                rotate=True)
    '''
    #style1.rotation = 90
    path1 = settings.BASE_DIR + '/assets/temp/ortho_tree2.svg'
    asset_path1 = '/temp/ortho_tree2.svg'
    tree1.render(path1, dpi=800, tree_style=style1)
    envoi = True

    return render(request, 'chlamdb/transporters_families.html', my_locals(locals()))


def transporters(request):
    from chlamdb.phylo_tree_display import ete_motifs
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    transporters_form = transporters_superfam_form(biodb)

    if request.method == 'POST': 
        form = transporters_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import blast_heatmap
            #if request.method == 'POST': 

            transporter_superfamily = form.cleaned_data['transporter_superfamily'][0]
            score_cutoff = form.cleaned_data['score_cutoff']
            evalue_cutoff = form.cleaned_data['evalue_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']
            hit_coverage_cutoff = form.cleaned_data['hit_coverage_cutoff']

            from ete3 import Tree
            sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                       ' where t2.name="%s";' % biodb
                       
            server, db = manipulate_biosqldb.load_db(biodb)
            tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
            t1 = Tree(tree)

            acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
            filter = '"' + '","'.join(acc_list) + '"'
            sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
            taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            t1 = Tree(tree)
            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            try:
                t2 = t1.prune(taxon_list)
            except:
                pass
            t1.ladderize()



            if transporter_superfamily == 'all':

                from chlamdb.plots import transporters_heatmap

                superfam_list, taxon2code2count = transporters_heatmap.transporter_all_superfamily_heatmap(biodb,
                                                                                                            evalue_cutoff,
                                                                                                            score_cutoff,
                                                                                                            query_coverage_cutoff,
                                                                                                            hit_coverage_cutoff)

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                    superfam_list,
                                                                    taxon2code2count,
                                                                    tree=t1,
                                                                    identity_scale=False,
                                                                    show_labels=True,
                                                                    column_scale=True,
                                                                    as_float=False,
                                                                    rotate=True)
                style1.rotation = 90
                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)
                envoi = True

            else:
                from chlamdb.plots import transporters_heatmap

                superfam_list, taxon2code2count = transporters_heatmap.transporter_superfamily_heatmap(biodb,
                                                                                                        transporter_superfamily,
                                                                                                        evalue_cutoff,
                                                                                                        score_cutoff,
                                                                                                        query_coverage_cutoff,
                                                                                                        hit_coverage_cutoff)

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                    superfam_list,
                                                                    taxon2code2count,
                                                                    tree=t1,
                                                                    identity_scale=False,
                                                                    show_labels=True,
                                                                    column_scale=True,
                                                                    as_float=False,
                                                                    rotate=True)
                style1.rotation = 90
                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)
                envoi = True



    else:  
        form = transporters_form()

    return render(request, 'chlamdb/transporters_superfam.html', my_locals(locals()))






def blast_sets(request):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)
    sets_form = blast_sets_form(biodb)

    if request.method == 'POST': 
        form = sets_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import blast_heatmap
            #if request.method == 'POST': 

            hmm_sets = form.cleaned_data['blast_set']
            score_cutoff = form.cleaned_data['score_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']
            hit_coverage_cutoff = form.cleaned_data['hit_coverage_cutoff']

            counts = request.POST['counts']

            if counts == 'detailed':
                blast_set_filter = '"'+'","'.join(hmm_sets)+'"'

                gene2taxon2score, gene_list = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'bitscore',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)
                if hmm_sets == 'flagellum':
                    gene_list = ['Flg_fliE','Flg_flgC', 'Flg_flgB','Flg_sctJ_FLG','Flg_sctN_FLG','Flg_sctQ_FLG','Flg_sctR_FLG',
                                 'Flg_sctS_FLG','Flg_sctT_FLG','Flg_sctU_FLG','Flg_sctV_FLG']

                sql = 'select accession, t1.name from blast.blast_sets t1 inner join blast.blast_sets_entry t2 on t1.set_id=t2.set_id ' \
                      ' inner join blast.blast_db t3 on t2.seq_id=t3.seq_id where t1.name in (%s) order by t1.name;' % (blast_set_filter)

                ordered_gene_list = []
                for row in server.adaptor.execute_and_fetchall(sql):
                    ordered_gene_list.append('%s (%s)' % (row[0], row[1]))

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                     ordered_gene_list,
                                                                     gene2taxon2score,
                                                                     identity_scale=True,
                                                                     show_labels=True,
                                                                     column_scale=True)

                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)


                gene2taxon2score2, gene_list2 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'query_coverage',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list2,
                                                            gene2taxon2score2,
                                                            identity_scale=True,
                                                            show_labels=True,
                                                            column_scale=True)

                path2 = settings.BASE_DIR + '/assets/temp/ortho_tree2.svg'
                asset_path2 = '/temp/ortho_tree2.svg'
                tree2.render(path2, dpi=800, tree_style=style2)

                gene2taxon2score3, gene_list3 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'locus_tag',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree3, style3 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list3,
                                                            gene2taxon2score3,
                                                            identity_scale=False,
                                                            show_labels=True,
                                                            column_scale=False)

                path3 = settings.BASE_DIR + '/assets/temp/ortho_tree3.svg'
                asset_path3 = '/temp/ortho_tree3.svg'
                tree3.render(path3, dpi=800, tree_style=style3)

                gene2taxon2score4, gene_list4 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'identity',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree4, style4 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list4,
                                                            gene2taxon2score4,
                                                            identity_scale=True,
                                                            show_labels=True,
                                                            column_scale=True)

                path4 = settings.BASE_DIR + '/assets/temp/ortho_tree4.svg'
                asset_path4 = '/temp/ortho_tree4.svg'
                tree4.render(path4, dpi=800, tree_style=style4)


                gene2taxon2score5, gene_list5 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'evalue',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree5, style5 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list5,
                                                            gene2taxon2score5,
                                                            identity_scale=False,
                                                            show_labels=True,
                                                            column_scale=True,
                                                                     as_float=True)

                path5 = settings.BASE_DIR + '/assets/temp/ortho_tree5.svg'
                asset_path5 = '/temp/ortho_tree5.svg'
                tree5.render(path5, dpi=800, tree_style=style5)

                envoi = True
            else:
                from ete3 import Tree
                sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                           ' where t2.name="%s";' % biodb
                           
                server, db = manipulate_biosqldb.load_db(biodb)
                tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]

                #acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
                #filter = '"' + '","'.join(acc_list) + '"'
                #sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
                #taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
                t1 = Tree(tree)
                R = t1.get_midpoint_outgroup()
                t1.set_outgroup(R)
                #try:
                #    t2 = t1.prune(taxon_list)
                #except:
                #    pass
                t1.ladderize()

                set2taxon2count = blast_heatmap.get_multiple_set_counts(biodb,
                                                                        hmm_sets,
                                                                        score_cutoff,
                                                                        query_coverage_cutoff,
                                                                        hit_coverage_cutoff)

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                     hmm_sets,
                                                                     set2taxon2count,
                                                                     tree=t1,
                                                                     identity_scale=False,
                                                                     show_labels=True,
                                                                     column_scale=True,
                                                                     as_float=False,
                                                                     rotate=True)
                style1.rotation = 90
                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)
                envoi = True





    else:  
        form = sets_form()

    return render(request, 'chlamdb/blast_sets_profiles.html', my_locals(locals()))



def hmm(request):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)
    hmm_form = hmm_sets_form(biodb)

    if request.method == 'POST': 
        form = hmm_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import hmm_heatmap
            #if request.method == 'POST': 

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            hmm_set = form.cleaned_data['hmm_set']
            score_cutoff = form.cleaned_data['score_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']



            if hmm_set == 'all':
                 sql = 'select t1.*,t2.orthogroup from custom_tables_annot_table as t1 inner ' \
                         ' join orthology_detail as t2 on t1.locus_tag=t2.locus_tag;'
            else:

                gene2taxon2score, gene_list = hmm_heatmap.get_single_set_data(biodb,
                                                                              hmm_set,
                                                                              'bitscore',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff)
                if hmm_set == 'flagellum':
                    gene_list = ['Flg_fliE','Flg_flgC', 'Flg_flgB','Flg_sctJ_FLG','Flg_sctN_FLG','Flg_sctQ_FLG','Flg_sctR_FLG',
                                 'Flg_sctS_FLG','Flg_sctT_FLG','Flg_sctU_FLG','Flg_sctV_FLG']

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                     gene_list,
                                                                     gene2taxon2score,
                                                                     identity_scale=True,
                                                                     show_labels=True,
                                                                     column_scale=True)

                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)




                gene2taxon2score2, gene_list2 = hmm_heatmap.get_single_set_data(biodb,
                                                                               hmm_set,
                                                                               'query_coverage',
                                                                               score_cutoff,
                                                                               query_coverage_cutoff)


                tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                     gene_list2,
                                                                     gene2taxon2score2,
                                                                     identity_scale=True,
                                                                     show_labels=True,
                                                                     column_scale=True)

                path2 = settings.BASE_DIR + '/assets/temp/ortho_tree2.svg'
                asset_path2 = '/temp/ortho_tree2.svg'
                tree2.render(path2, dpi=800, tree_style=style2)

                gene2taxon2score3, gene_list3 = hmm_heatmap.get_single_set_data(biodb,
                                                                                hmm_set,
                                                                                'locus_tag',
                                                                                score_cutoff,
                                                                                query_coverage_cutoff)


                tree3, style3 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list3,
                                                            gene2taxon2score3,
                                                            identity_scale=False,
                                                            show_labels=True,
                                                            column_scale=False)

                path3 = settings.BASE_DIR + '/assets/temp/ortho_tree3.svg'
                asset_path3 = '/temp/ortho_tree3.svg'
                tree3.render(path3, dpi=800, tree_style=style3)


                envoi = True

    else:  
        form = hmm_form()

    return render(request, 'chlamdb/hmm_profiles.html', my_locals(locals()))


def locus_int(request):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)
    module_int_form = locus_int_form(biodb)

    if request.method == 'POST': 
        form = module_int_form(request.POST)
        if form.is_valid():
            #if request.method == 'POST': 

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' 

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            category = form.cleaned_data['category']

            show_identity = request.POST['show_id']
            if show_identity=="noshow":
                identity_heatmap =False
            else:
                identity_heatmap=True


            if category == 'all':
                 sql = 'select t1.*,t2.orthogroup from custom_tables_annot_table as t1 inner ' \
                         ' join orthology_detail as t2 on t1.locus_tag=t2.locus_tag;' 

            else:

                 sql = 'select t1.*,t2.orthogroup from custom_tables_annot_table as t1 inner ' \
                         ' join orthology_detail as t2 on t1.locus_tag=t2.locus_tag ' \
                         ' where category="%s";' % (category) # where pathway_category!="1.0 Global and overview maps"


            data = server.adaptor.execute_and_fetchall(sql,)

            orthogroups = set([i[-1] for i in data])
            locus_tag_list = [i[2] for i in data]

            locus2gene = {}
            for i in data:
                locus2gene[i[2]] = i[1]

            filter = '"'+'","'.join(locus_tag_list)+'"'
            sql = 'select locus_tag, taxon_id from orthology_detail where locus_tag in (%s)' % (filter)

            locus2taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


            taxon2locus2identity_closest = ete_motifs.get_locus2taxon2identity(biodb, locus_tag_list)
            # format
            # 'PC_RS08835': {'1279787': 35.32, '1279808': 35.49, '1069694'}
            locus2taxon2identity_closest = {}
            labels2 = []
            locus_edit2taxon = {}
            for locus_tag in taxon2locus2identity_closest:
                new_label = '%s / %s' % (locus_tag, locus2gene[locus_tag])
                locus_edit2taxon[new_label] = locus2taxon[locus_tag]
                locus2taxon2identity_closest[new_label] = taxon2locus2identity_closest[locus_tag]
            for locus_tag in locus_tag_list:
                new_label = '%s / %s' % (locus_tag, locus2gene[locus_tag])
                labels2.append(new_label)
            #print 'PC_RS08840 / SdhC' in labels2
            #print locus2taxon2identity_closest['PC_RS08840 / SdhC']
            taxon2orthogroup2count_all = ete_motifs.get_taxon2orthogroup2count(biodb, orthogroups)

            labels = orthogroups

            tree1, style = ete_motifs.multiple_profiles_heatmap(biodb, labels, taxon2orthogroup2count_all)
            #labels = locus_tag_list
            tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                        labels2,
                                                        locus2taxon2identity_closest,
                                                        identity_scale=True,
                                                        show_labels=identity_heatmap,
                                                        reference_taxon=locus_edit2taxon,
                                                        rotate=False)
            scale_path = "/scales/scale_identity_red.png"


            #except:
            #    tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, merged_dico)

            if len(orthogroups) > 1000:
                big = True
                path = settings.BASE_DIR + '/assets/temp/ortho_tree.png'
                path2 = settings.BASE_DIR + '/assets/temp/ortho_tree2.png'
                asset_path = '/temp/ortho_tree.png'
                asset_path2 = '/temp/ortho_tree2.png'
                tree1.render(path, dpi=1200, tree_style=style)
                tree2.render(path2, dpi=1200, tree_style=style2)
            else:
                big = False
                path = settings.BASE_DIR + '/assets/temp/ortho_tree.svg'
                path2 = settings.BASE_DIR + '/assets/temp/ortho_tree2.svg'
                asset_path = '/temp/ortho_tree.svg'
                asset_path2 = '/temp/ortho_tree2.svg'
                tree1.render(path, dpi=800, tree_style=style)
                tree2.render(path2, dpi=800, tree_style=style2)
            envoi = True

    else:  
        form = module_int_form()

    return render(request, 'chlamdb/inter_tree.html', my_locals(locals()))




def kegg_pathway_heatmap(request):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs
    server, db = manipulate_biosqldb.load_db(biodb)
    pathway_form = make_pathway_overview_form(biodb)#get_locus_annotations_form(biodb)

    if request.method == 'POST': 
        form = pathway_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import pathway_heatmap
            from ete3 import Tree
            #if request.method == 'POST': 

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            category = form.cleaned_data['category']
            sql = 'select pathway_name from enzyme_kegg_pathway where pathway_category="%s";' % (category)
            modules = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                       ' where t2.name="%s";' % biodb
                       
            server, db = manipulate_biosqldb.load_db(biodb)
            tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]

            #acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
            #filter = '"' + '","'.join(acc_list) + '"'
            #sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
            #taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            t1 = Tree(tree)
            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            #try:
            #    t2 = t1.prune(taxon_list)
            #except:
            #    pass
            t1.ladderize()

            tree, style = pathway_heatmap.plot_pathway_heatmap(biodb,
                                                               t1,
                                                               modules,
                                                               taxon_id_list = [],
                                                               rotate=True)
            style.rotation = 90




            path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
            asset_path = '/temp/metabo_tree.svg'
            tree.render(path, dpi=800, tree_style=style)


            envoi = True


    else:  
        form = pathway_form()

    return render(request, 'chlamdb/pathway_cat.html', my_locals(locals()))


def kegg_module_subcat(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    module_overview_form = make_module_overview_form(db, True)
    if request.method != "POST":
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        # TODO: add error message
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    category        = form.cleaned_data["category"]
    leaf_to_name    = db.get_genomes_description()
    ko_count_subcat = db.get_ko_count_cat(category=category)

    grouped_count = ko_count_subcat.groupby(["taxon_id", "module_id"]).sum()

    unique_module_ids = ko_count_subcat.index.get_level_values("module_id").unique().tolist()
    if len(unique_module_ids) == 0:
        # add error message : no module found
        envoi = True
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    module_infos = db.get_modules_info(unique_module_ids)
    expression_tree = {}
    for module_id, descr, definition, *other in module_infos:
        parser = ModuleParser(definition)
        expression_tree[module_id] = parser.parse()

    labels = [format_module(val[0]) for val in module_infos]

    grouped_count = grouped_count.unstack(level=1, fill_value=0)
    grouped_count.columns = [col for col in grouped_count["count"].columns]
    ref_tree = db.get_reference_phylogeny()
    e_tree = EteTree.default_tree(ref_tree)
    for module, counts in grouped_count.iteritems():
        values = counts.to_dict()
        header = format_module(module)
        expr_tree = expression_tree[module]
        n_missing = {}
        for bioentry, count in counts.iteritems():
            # hack for now
            s = ko_count_subcat.loc[bioentry].index.get_level_values("KO").tolist()
            n_missing[bioentry] = expr_tree.get_n_missing({ko: 1 for ko in s})
        new_col = KOAndCompleteness(values, n_missing, header)
        e_tree.add_column(new_col)
    e_tree.rename_leaves(leaf_to_name.description.to_dict())
    path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
    asset_path = '/temp/metabo_tree.svg'
    e_tree.render(path, dpi=500, w=800)
    envoi = True
    return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))


def kegg_module(request):
    from chlamdb.phylo_tree_display import ete_motifs
    from chlamdb.plots import module_heatmap

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    module_overview_form = make_module_overview_form(db)

    if request != "POST":
        form = module_overview_form()
        return render(request, 'chlamdb/module_overview.html', my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        # TODO: add error message
        form = module_overview_form()
        return render(request, 'chlamdb/module_overview.html', my_locals(locals()))

    envoie = True
    return render(request, 'chlamdb/module_overview.html', my_locals(locals()))


def module2heatmap(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_kegg_form(biodb)

    if request.method == 'POST': 
        form = comp_metabo_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import pathway_heatmap
            from chlamdb.biosqldb import biosql_own_sql_tables
            from chlamdb.phylo_tree_display import ete_motifs
            from ete3 import Tree

            pathway_category = form.cleaned_data['pathway_choice']
            pathway_filter = '"' + '","'.join(pathway_category) + '"'
            sql = 'select pathway_name from enzyme_kegg_pathway where description in (%s);' % (pathway_filter)

            pathways = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            module_category = form.cleaned_data['module_choice']
            module_filter = '"' + '","'.join(module_category) + '"'
            sql = 'select module_name from enzyme_kegg_module_v1 where description in (%s);' % (module_filter)

            modules = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

            sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                       ' where t2.name="%s";' % biodb
            server, db = manipulate_biosqldb.load_db(biodb)
            tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]

            acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
            filter = '"' + '","'.join(acc_list) + '"'
            sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
            taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            t1 = Tree(tree)
            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            try:
                t2 = t1.prune(taxon_list)
            except:
                pass
            t1.ladderize()

            modules = [
"M00051",
"M00052",
"M00053",
"M00048",
"M00049",
"M00050",
"M00115",
"M00116",
"M00117",
"M00119",
"M00120",
"M00121",
"M00123",
"M00125",
"M00126",
"M00127",
"M00122",
"M00133",
"M00134",
            ]

            pathways = ["map00230",
"map00240",
"map00220",
"map00250",
"map00260",
"map00270",
"map00280",
"map00290",
"map00300",
"map00310",
"map00330",
"map00340",
"map00350",
"map00360",
"map00380",
"map00400",

]

            tree, style = pathway_heatmap.plot_module_heatmap(biodb,
                                                                                     t1,
                                                                                     pathways,
                                                                                     modules,
                                                                                     [],
                                                                                     rotate=True)

            big = False
            path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
            asset_path = '/temp/metabo_tree.svg'

            style.show_leaf_name = False
            style.rotation = 90

            tree.render(path, dpi=200, w=900, tree_style=style)

            envoi = True

    else:  
        form = comp_metabo_form()

    return render(request, 'chlamdb/module2heatmap.html', my_locals(locals()))


def module_comparison(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    comp_metabo_form = make_metabo_from(db)

    if request != "POST":
        form = comp_metabo_form()
        return render(request, 'chlamdb/module_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        # TODO: add error message
        form = comp_metabo_form()
        return render(request, 'chlamdb/module_comp.html', my_locals(locals()))

    envoie = True
    return render(request, 'chlamdb/module_comp.html', my_locals(locals()))


def module_comparison_legacy(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST': 
        form = comp_metabo_form(request.POST)
        if form.is_valid():
            from chlamdb.biosqldb import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


            sql_pathway_count = 'select BB.module_name,count_all,count_db,count_db/count_all from (select module_id, count(*) ' \
                        ' as count_db from (select distinct ko_id from enzyme_locus2ko) as t1' \
                        ' inner join enzyme_module2ko_v1 as t2 on t1.ko_id=t2.ko_id group by module_id) AA ' \
                        ' right join (select t1.module_id,module_name, count_all from (select module_id, count(*) as count_all ' \
                        'from enzyme_module2ko_v1 group by module_id) t1 inner join enzyme_kegg_module_v1 as t2 ' \
                        'on t1.module_id=t2.module_id)BB on AA.module_id=BB.module_id;'

            map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))

            category2maps = {}

            sql_category2maps = 'select module_sub_cat,module_name,description from enzyme_kegg_module_v1;'

            data = server.adaptor.execute_and_fetchall(sql_category2maps,)

            for one_map in data:
                if one_map[0] not in category2maps:
                    category2maps[one_map[0]] = [[one_map[1], one_map[2]]]
                else:
                    category2maps[one_map[0]].append([one_map[1], one_map[2]])

            sql = 'select distinct module_name,module_sub_sub_cat from enzyme_kegg_module_v1;'
            module2sub_category = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            taxon_maps = []
            for taxon in taxon_list:
                database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

                sql = 'select module_name, n from (select B.module_id,count(*) as n from ' \
                      ' (select * from enzyme_locus2ko where taxon_id=%s) A ' \
                      ' left join enzyme_module2ko_v1 as B on A.ko_id=B.ko_id group by module_id) AA ' \
                      ' right join enzyme_kegg_module_v1 as BB on AA.module_id=BB.module_id;' % (taxon)

                map2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                taxon_maps.append(map2count_taxon)


            envoi_comp = True

    else:  
        form = comp_metabo_form()

    return render(request, 'chlamdb/module_comp.html', my_locals(locals()))




def metabo_overview(request):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

    sql_pathway_count = 'select PATH2.pathway_name,PATH2.n,PATH1.n,PATH1.n/PATH2.n from (select pathway_name,count(*) as n ' \
                        ' from (select distinct pathway_name,t4.ec_id from enzyme_locus2ec as t3 ' \
                        ' inner join enzyme_kegg2ec as t4 on t3.ec_id=t4.ec_id inner join enzyme_kegg_pathway as t5 ' \
                        ' on t4.pathway_id=t5.pathway_id where pathway_category!="1.0 Global and overview maps") A ' \
                        ' group by pathway_name) ' \
                        ' PATH1 right join ' \
                        ' (select pathway_name,count(*) as n from enzyme_kegg2ec as t1 ' \
                        ' inner join enzyme_kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                        '  group by pathway_name) ' \
                        ' PATH2 on PATH2.pathway_name=PATH1.pathway_name;'  # where pathway_category!="1.0 Global and overview maps"

    map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))

    sql = 'select C.pathway_category,taxon_id, A.pathway_name,A.n_enzymes, C.description from ' \
          '( select distinct taxon_id,pathway_name, count(*) as n_enzymes from (' \
          'select distinct taxon_id, ec_id  from enzyme_locus2ec as b1 ' \
          'inner join bioentry as b2 on b1.accession=b2.accession where biodatabase_id=%s) ' \
          't1 inner join enzyme_kegg2ec as t2  on t1.ec_id=t2.ec_id ' \
          'inner join enzyme_kegg_pathway as t3 on t2.pathway_id=t3.pathway_id ' \
          'group by taxon_id,pathway_name) A ' \
          'inner join enzyme_kegg_pathway as C  on A.pathway_name=C.pathway_name;' % (database_id)

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
    '''
    tree = ete_motifs.pathways_heatmap(biodb,
                                       category2maps,
                                       pathway_category2taxon2map)


    #except:
    #    tree = ete_motifs.multiple_profiles_heatmap(biodb, labels, merged_dico)

    if len(all_maps) > 1000:
        big = True
        path = settings.BASE_DIR + '/assets/temp/metabo_tree.png'
        asset_path = '/temp/metabo_tree.png'
        tree.render(path, dpi=1200, h=600)
    else:
        big = False
        path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
        asset_path = '/temp/metabo_tree.svg'

        tree.render(path, dpi=800, h=600)
    '''
    envoi = True

    #else:  
    #    pass

    return render(request, 'chlamdb/metabo_overview.html', my_locals(locals()))



def metabo_comparison(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST': 
        form = comp_metabo_form(request.POST)
        if form.is_valid():
            from chlamdb.biosqldb import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


            sql_pathway_count = 'select PATH2.pathway_name,PATH2.n,PATH1.n,PATH1.n/PATH2.n from (select pathway_name,count(*) as n ' \
                                ' from (select distinct pathway_name,t4.ec_id from enzyme_locus2ec as t3 ' \
                                ' inner join enzyme_kegg2ec as t4 on t3.ec_id=t4.ec_id inner join enzyme_kegg_pathway as t5 ' \
                                ' on t4.pathway_id=t5.pathway_id where pathway_category!="1.0 Global and overview maps") A ' \
                                ' group by pathway_name) ' \
                                ' PATH1 right join ' \
                                ' (select pathway_name,count(*) as n from enzyme_kegg2ec as t1 ' \
                                ' inner join enzyme_kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                                ' where pathway_category!="1.0 Global and overview maps" group by pathway_name) ' \
                                ' PATH2 on PATH2.pathway_name=PATH1.pathway_name;'

            map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))

            category2maps = {}

            sql_category2maps = 'select pathway_category,pathway_name,description from  enzyme_kegg2ec as t1 ' \
                                ' inner join enzyme_kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                                ' where pathway_category!="1.0 Global and overview maps" group by pathway_name;'

            data = server.adaptor.execute_and_fetchall(sql_category2maps,)

            for one_map in data:
                if one_map[0] not in category2maps:
                    category2maps[one_map[0]] = [[one_map[1], one_map[2]]]
                else:
                    category2maps[one_map[0]].append([one_map[1], one_map[2]])

            taxon_maps = []
            for taxon in taxon_list:

                sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

                biodatabase_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

                sql = 'select PATH2.pathway_name,PATH1.n from (select pathway_name,count(*) as n from ' \
                      ' (select distinct pathway_name,t4.ec_id from enzyme_locus2ec as t3 ' \
                      ' inner join enzyme_kegg2ec as t4 on t3.ec_id=t4.ec_id inner join enzyme_kegg_pathway as t5 ' \
                      ' on t4.pathway_id=t5.pathway_id inner join bioentry as t6 on t3.accession=t6.accession ' \
                      ' where biodatabase_id=%s and pathway_category!="1.0 Global and overview maps" and t6.taxon_id=%s) A ' \
                      ' group by pathway_name) PATH1 right join (select pathway_name,count(*) as n ' \
                      ' from enzyme_kegg2ec as t1 inner join enzyme_kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                      ' where pathway_category!="1.0 Global and overview maps" ' \
                      ' group by pathway_name) PATH2 on PATH2.pathway_name=PATH1.pathway_name order by n;' % (biodatabase_id, taxon)
                map2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                taxon_maps.append(map2count_taxon)


            envoi_comp = True

    else:  
        form = comp_metabo_form()

    return render(request, 'chlamdb/metabo_comp.html', my_locals(locals()))


def metabo_comparison_ko(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST': 
        form = comp_metabo_form(request.POST)
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            from chlamdb.biosqldb import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


            sql_pathway_count = 'select PATH2.pathway_name,PATH2.n,PATH1.n,PATH1.n/PATH2.n from (select pathway_name,count(*) ' \
                                ' as n from (select ko_id from enzyme_locus2ko as t1 ' \
                                ' group by ko_id) A inner join enzyme_pathway2ko_v1 as B on A.ko_id=B.ko_id inner join ' \
                                ' enzyme_kegg_pathway as C on B.pathway_id=C.pathway_id ' \
                                ' where pathway_category!="1.0 Global and overview maps" group by pathway_name) ' \
                                ' PATH1 right join ' \
                                ' (select pathway_name,count(*) as n from enzyme_pathway2ko_v1 as t1 inner join enzyme_kegg_pathway as ' \
                                ' t2 on t1.pathway_id=t2.pathway_id where pathway_category!="1.0 Global and overview maps" ' \
                                ' group by pathway_name)' \
                                ' PATH2 on PATH2.pathway_name=PATH1.pathway_name;'

            map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))

            category2maps = {}

            sql_category2maps = 'select pathway_category,pathway_name,description from  enzyme_pathway2ko_v1 as t1 ' \
                                ' inner join enzyme_kegg_pathway as t2 on t1.pathway_id=t2.pathway_id ' \
                                ' where pathway_category!="1.0 Global and overview maps" group by pathway_name;'

            data = server.adaptor.execute_and_fetchall(sql_category2maps,)

            for one_map in data:
                if one_map[0] not in category2maps:
                    category2maps[one_map[0]] = [[one_map[1], one_map[2]]]
                else:
                    category2maps[one_map[0]].append([one_map[1], one_map[2]])

            taxon_maps = []
            for taxon in taxon_list:

                sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

                biodatabase_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

                sql = 'select PATH2.pathway_name,PATH1.n from (select pathway_name,count(*) as n from (select ko_id ' \
                      ' from enzyme_locus2ko as t1 where taxon_id=%s group by ko_id) A inner join enzyme_pathway2ko_v1 ' \
                      ' as B on A.ko_id=B.ko_id inner join enzyme_kegg_pathway as C on B.pathway_id=C.pathway_id ' \
                      ' where pathway_category!="1.0 Global and overview maps" group by pathway_name) PATH1 ' \
                      ' right join (select pathway_name,count(*) from enzyme_pathway2ko_v1 as t1 inner join enzyme_kegg_pathway as ' \
                      ' t2 on t1.pathway_id=t2.pathway_id where pathway_category!="1.0 Global and overview maps" ' \
                      ' group by pathway_name) PATH2 on PATH2.pathway_name=PATH1.pathway_name order by n;' % (taxon)
                map2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
                taxon_maps.append(map2count_taxon)

            envoi_comp = True

    else:  
        form = comp_metabo_form()

    return render(request, 'chlamdb/metabo_comp_ko.html', my_locals(locals()))


def pfam_comparison(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    comp_metabo_form = make_metabo_from(biodb)

    if request.method == 'POST': 
        form = comp_metabo_form(request.POST)
        #form2 = ContactForm(request.POST)
        if form.is_valid():
            from chlamdb.biosqldb import biosql_own_sql_tables
            taxon_list = form.cleaned_data['targets']

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

            columns = '`' + '`,`'.join(taxon_list) + '`'
            filter = '(`' + '`>0 or`'.join(taxon_list) + '`>0)'


            sql = 'select * from (select id from comparative_tables_Pfam where %s group by id) A' \
                  ' inner join (select distinct signature_accession,signature_description,count(*) as n ' \
                  ' from interpro where analysis="Pfam" group by signature_accession,signature_description) B on A.id = B.signature_accession' % (filter)

            sql_pathway_count = 'select distinct signature_accession,signature_description,count(*) as n ' \
                                ' from interpro where analysis="Pfam" group by signature_accession,signature_description;'

            pfam_data_raw = server.adaptor.execute_and_fetchall(sql,)
            pfam2data = {}
            for one_pfam_entry in pfam_data_raw:
                pfam2data[one_pfam_entry[0]] = one_pfam_entry[1:]




            taxon_dicos = []
            for taxon in taxon_list:

                sql = 'select distinct signature_accession,count(*) as n ' \
                      ' from interpro where analysis="Pfam" and taxon_id=%s ' \
                      ' group by signature_accession;' % (taxon)

                accession2count_taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


                for accession in pfam2data:
                    if accession not in accession2count_taxon:
                        accession2count_taxon[accession] = 0
                taxon_dicos.append(accession2count_taxon)

            envoi_comp = True

    else:  
        form = comp_metabo_form()

    return render(request, 'chlamdb/pfam_comp.html', my_locals(locals()))


def orthogroup_comparison(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    comp_metabo_form = make_metabo_from(db)
    if request.method != 'POST': 
        form = comp_metabo_form(request.POST)
        return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))

    try:
        all_targets = form.get_choices()
    except:
        # TODO: add error message
        return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))

    genomes = db.get_genomes_description().description.to_dict()
    og_count = db.get_og_count(all_targets)
    og_count.columns = [genomes[int(col)] for col in og_count.columns]
    annotations = db.get_genes_from_og(orthogroups=og_count.index.tolist(),
            taxon_ids=all_targets, terms=["product"])

    products = annotations.groupby("orthogroup")["product"].apply(list)
    n_orthogroups = len(og_count.index)

    og_data = []
    for og, items in og_count.iterrows():
        piece = [format_orthogroup(og)]
        piece.append(items.tolist())
        if og in products.index:
            piece.append(format_lst_to_html(products.loc[og]))
        else:
            piece.append("-")
        og_data.append(piece)

    genomes_list = og_count.columns.tolist()
    envoi_comp = True
    return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))


def ko_comparison(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    comp_metabo_form = make_metabo_from(db)

    if request.method != "POST":
        form = comp_metabo_form()
        return render(request, 'chlamdb/ko_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/ko_comp.html', my_locals(locals()))

    include = form.get_choices()
    mat_include = db.get_ko_count(include).unstack(level=0, fill_value=0)
    mat_include.columns = [col for col in mat_include["count"].columns.values]

    ko2annot = db.get_ko_desc(mat_include.index.tolist())
    df_ttl   = db.get_ko_count(mat_include.index.tolist(), search_on="ko_id")
    ko2total_count = df_ttl.groupby("KO").sum()["count"].to_dict()
    ko2counts = mat_include.to_dict()
    ko2counts = {}
    ko2_print = {}
    for key, values in mat_include.iterrows():
        ko2counts[key] = values.values.tolist()
        ko2_print[key] = format_ko(key)

    hsh_gen_desc = db.get_genomes_description().description.to_dict()
    taxon_list = [hsh_gen_desc[int(col)] for col in mat_include.columns.values]
    n_ko = len(mat_include.index.tolist())
    envoi_comp = True
    return render(request, 'chlamdb/ko_comp.html', my_locals(locals()))


def faq(request):
    a =2
    return render(request, 'chlamdb/FAQ.html', my_locals(locals()))

def phylogeny_intro(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100*x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x/pow(10,6), 2))

    data_table_header = ["Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot", "n_contigs", "length", "coding_density"]].values.tolist()

    # plot phylo only of not already in assets
    tree = db.get_reference_phylogeny()
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    if not R is None:
        t1.set_outgroup(R)
    t1.ladderize()

    e_tree = EteTree(t1)
    header_params = {"rotation": -30}
    stacked_face_params = {"margin_right": 8, "margin_left": 5}

    tree_params = [
            # serie_name, header, color and is_relative
            ["length", "Size (Mbp)", "#91bfdb", True],
            ["gc", "GC %", "#fc8d59", False],
            ["coding_density", "Coding density %", "#99d594", False],
            ["completeness", "Completeness", "#d7191c", False],
            ["contamination", "Contamination", "black", False]]
        
    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours = [col, "white"],
                relative=is_relative, header=header,
                header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/phylogeny_intro.html', my_locals(locals()))

def genomes_intro(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100*x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x/pow(10,6), 2))

    data_table_header = ["Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot", "n_contigs", "length", "coding_density"]].values.tolist()

    # plot phylo only of not already in assets
    tree = db.get_reference_phylogeny()
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    if not R is None:
        t1.set_outgroup(R)
    t1.ladderize()

    e_tree = EteTree(t1)
    header_params = {"rotation": -30}
    stacked_face_params = {"margin_right": 8, "margin_left": 5}

    tree_params = [
            # serie_name, header, color and is_relative
            ["length", "Size (Mbp)", "#91bfdb", True],
            ["gc", "GC %", "#fc8d59", False],
            ["coding_density", "Coding density %", "#99d594", False],
            ["completeness", "Completeness", "#d7191c", False],
            ["contamination", "Contamination", "black", False]]
        
    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours = [col, "white"],
                relative=is_relative, header=header,
                header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/genomes_intro.html', my_locals(locals()))

    