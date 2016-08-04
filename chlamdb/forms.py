#-*- coding: utf-8 -*-
from django import forms
#import models
from manipulate_biosqldb import load_db, get_biodatabase_list
#from blah import *
#from models import GenDB
#from models import Genome
from django.forms import ModelForm

server = load_db()

sql ="select accession from bioentry where biodatabase_id = 22"
result = server.adaptor.execute_and_fetchall(sql, )
accession_list = [i[0] for i in result]
accession_choices = []
for i in accession_list:
    accession_choices.append((i,i))


def get_accessions(database_name, all=False, plasmid=False):

    import manipulate_biosqldb
    server = manipulate_biosqldb.load_db()
    if not plasmid:
        print "no plasmid"
        sql ='SELECT bioentry.taxon_id, bioentry.description FROM bioentry ' \
             'inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id ' \
             'where biodatabase.name ="%s"and bioentry.description not like "%%%%plasmid%%%%" and bioentry.description not like "%%%%phage%%%%"' \
             'order by bioentry.description' % database_name #
    else:
        print 'plasmid'
        sql ='SELECT bioentry.accession, bioentry.description FROM bioentry ' \
             'inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id ' \
             'where biodatabase.name ="%s"' \
             'order by bioentry.description' % database_name
    result = server.adaptor.execute_and_fetchall(sql, )
    accession_list = [i for i in result]
    print "acc", accession_list
    accession_choices = []



    for accession in accession_list:
        accession_choices.append((accession[0], accession[1]))


    import re
    accessions = {}
    for i, accession in enumerate(accession_choices):
        print i, accession
        description = accession[1]
        description = re.sub(", complete genome\.", "", description)
        description = re.sub(", complete genome", "", description)
        description = re.sub(", complete sequence\.", "", description)
        description = re.sub("strain ", "", description)
        description = re.sub("str\. ", "", description)
        description = re.sub(" complete genome sequence\.", "", description)
        description = re.sub(" complete genome\.", "", description)
        description = re.sub(" chromosome", "", description)
        description = re.sub(" DNA", "S.", description)
        description = re.sub("Merged record from ", "", description)
        description = re.sub(", wgs", "", description)
        description = re.sub("Candidatus ", "", description)
        description = re.sub(".contig.0_1, whole genome shotgun sequence.", "", description)
        accession_choices[i] = (accession[0], description)
        accessions[description] = accession[0]

    accession_choices = []
    for description in sorted(accessions.keys()):
        accession_choices.append([accessions[description], description])

    if all:
        accession_choices = [["all", "all"]] + accession_choices

    return accession_choices

biodatabases = [ i for i in get_biodatabase_list(server)]

choices = []
for i in biodatabases:
    choices.append((i,i))

choices = tuple(choices)


def make_contact_form(server, database_name):

    #accession_choices = get_accessions(database_name)


    class ContactForm(forms.Form):

        accession = forms.CharField(max_length=100)
        #biodatabase = forms.ChoiceField(choices=choices)
        #plot_region = forms.NullBooleanField(widget=forms.CheckboxInput())
        #region_size = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
        #genomes = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "30" }), required = False)
    return ContactForm

def make_plot_form(database_name):

    accession_choices = get_accessions(database_name)

    class PlotForm(forms.Form):
        accession = forms.CharField(max_length=100)
        #location_start = forms.CharField(max_length=9, label="sart (bp)", required=False)
        #location_stop = forms.CharField(max_length=9, label="end (bp)", required=False)
        region_size = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
        genomes = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "15" }), required = False)
    return PlotForm



def make_interpro_from(database_name):

    accession_choices = get_accessions(database_name)

    '''
    database_choices =  [["Coils", "Coils"],
                        ["Gene3D", "Gene3D"],
                        ["Hamap", "Hamap"],
                        ["Pfam", "Pfam"],
                        ["Phobius", "Phobius"],
                        ["PIRSF", "PIRSF"],
                        ["PRINTS", "PRINTS"],
                        ["ProDom", "ProDom"],
                        ["ProSitePatterns", "ProSitePatterns"],
                        ["ProSiteProfiles", "ProSiteProfiles"],
                        ["SignalP_EUK", "SignalP_EUK"],
                        ["SignalP_GRAM_NEGATIVE", "SignalP_GRAM_NEGATIVE"],
                        ["SignalP_GRAM_POSITIVE", "SignalP_GRAM_POSITIVE"],
                        ["SMART", "SMART"],
                        ["SUPERFAMILY", "SUPERFAMILY"],
                        ["TIGRFAM", "TIGRFAM"]]
    '''
    class InterproForm(forms.Form):
        SEARCH_CHOICES = (('description', 'Description'), ('GO','GO Number'), ('EC','EC Number'), ('interpro_accession','Interpro Accession'))
        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % (len(accession_choices)/2) }), required = False)
        search_type = forms.ChoiceField(choices=SEARCH_CHOICES)
        search_term = forms.CharField(max_length=100)

    return InterproForm


def make_metabo_from(database_name):

    accession_choices = get_accessions(database_name)


    class MetaboForm(forms.Form):
        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % (len(accession_choices)/4) }), required = False)

    return MetaboForm


def make_venn_from(database_name):

    accession_choices = get_accessions(database_name)

    class VennForm(forms.Form):
        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'20' }), required = False)

        def clean_venn(self):
            value = self.cleaned_data['targets']
            if len(value) > 6:
                raise forms.ValidationError("You can't select more than 6 items.")
            return value
    return VennForm







class BiodatabaseForm(forms.Form):
    biodatabase = forms.ChoiceField(choices=choices, required = False)

    def save(self):
        self.biodatabase = self.cleaned_data["biodatabase"]




def make_extract_region_form(database_name):

    accession_choices = get_accessions(database_name, plasmid=True)

    extraction_choices = [['annotation', ''],['sequence', ''],['sequence_trans', '']]

    class ExtractRegionForm(forms.Form):



        genome = forms.ChoiceField(choices=accession_choices)
        region = forms.CharField(max_length=100, label="Region start, stop", initial = "1, 8000", required = False)
        extract = forms.ChoiceField(choices=extraction_choices, widget=forms.RadioSelect, label='')
        #get_annotation = forms.NullBooleanField(widget=forms.CheckboxInput())
        #get_sequence = forms.NullBooleanField(widget=forms.CheckboxInput())

    return ExtractRegionForm


def make_priam_form(database_name):

    accession_choices = get_accessions(database_name, plasmid=True)

    class PriamForm(forms.Form):
        genome = forms.ChoiceField(choices=accession_choices)

    return PriamForm


def make_circos_form(database_name):

    accession_choices = get_accessions(database_name)
    
    class CircosForm(forms.Form):
        circos_reference = forms.ChoiceField(choices=accession_choices)
        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'20'}), required = False)
        #get_region = forms.NullBooleanField(widget=forms.CheckboxInput())
        #region = forms.CharField(max_length=100, label="Region start, stop", initial = "1, 8000", required = False)
        
        def save(self):
            self.reference = self.cleaned_data["reference"]
            self.get_region = self.cleaned_data["get_region"]
            self.region = self.cleaned_data["region"]
    return CircosForm

def make_extract_form(database_name, plasmid=False):

    if not plasmid:
        accession_choices = get_accessions(database_name)
    else:
        accession_choices = get_accessions(database_name, plasmid=True)

    class ExtractForm(forms.Form):
        FREQ_CHOICES = ((0, 0),(1, 1), (2,2), (3,3), (4,4), (5,5), (6,6), (7,7), (8,8), (9,9), (10,10))


        orthologs_in = forms.MultipleChoiceField(label='', choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "17" }), required = False)
        no_orthologs_in = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "17" }), required = False, label="")

        new_choices = [['None', 'None']] + accession_choices
        frequency = forms.ChoiceField(choices=FREQ_CHOICES, label='')
        reference = forms.ChoiceField(choices=new_choices,label="")

    return ExtractForm

def make_module_overview_form(database_name):
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql = 'select distinct module_sub_cat from enzyme.kegg_module;'
    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [(i[0],i[0]) for i in categories]

    class ModuleCatChoice(forms.Form):
        category = forms.ChoiceField(choices=CHOICES)
    return ModuleCatChoice


class SearchForm(forms.Form):

    SEARCH_CHOICES = (('gene', 'gene'), ('product','product'), ('locus_tag','locus_tag'))

    search_type = forms.ChoiceField(choices=SEARCH_CHOICES)
    search_term = forms.CharField(max_length=100)
    #biodatabase = forms.ChoiceField(choices=choices)

def make_circos_orthology_form(biodb):

    accession_choices =  get_accessions(biodb, plasmid=True, all=True)

    class Circos_orthology(forms.Form):
        accession = forms.CharField(max_length=100)
        targets = forms.ChoiceField(choices=accession_choices)
    return Circos_orthology

def make_blast_form(biodb):

    accession_choices =  get_accessions(biodb, plasmid=True, all=True)

    class BlastForm(forms.Form):
        blast = forms.ChoiceField(choices=[("blastn_ffn", "blastn_ffn"),
                                           ("blastn_fna", "blastn_fna"),
                                           ("blastp", "blastp"),
                                           ("blastx", "blastx"),
                                           ("tblastn", "tblastn")])
        target = forms.ChoiceField(choices=accession_choices)
        blast_input = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 20}))
        #biodatabase = forms.ChoiceField(choices=choices)

    return BlastForm


class AnnotForm(forms.Form):
    orthogroups = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 10}))
    #biodatabase = forms.ChoiceField(choices=choices)



def make_motif_form(database_name):

    accession_choices = get_accessions(database_name, all=True)

    class MotifForm(forms.Form):
        '''

         choix fuzznuc / fuzzpro
        '''
        DATA_CHOICES = (('nucleotide', 'nucleotide'), ('protein','protein'))
        #search = forms.ChoiceField(choices=DATA_CHOICES)
        n_missmatch = forms.CharField(max_length=3, label="N. mismatches", required=False, initial=0)
        #input_file = forms.FileField(
        #    label='Select a file',
        #    help_text='max. 42 megabytes', required=False
        #)
        motif_input = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 20}),
                                      help_text="ex: Y-G-G-[LIV]-T-{I}-{N}-x(2)-N (PROSITE style patterns)",
                                      required=False)

        search_in = forms.ChoiceField(choices=accession_choices)
    return MotifForm


class PCRForm(forms.Form):
    '''
     primersearch
    '''

    input_file = forms.FileField(
        label='Select a file',
        help_text='max. 42 megabytes'
    )
    primersearch_input = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 20}))

def make_circos2genomes_form(database_name):

    accession_choices = get_accessions(database_name)

    class Circos2genomesForm(forms.Form):
        
        locus_list = forms.CharField(max_length=1000, required = False)
        reference_genome = forms.ChoiceField(choices=accession_choices)
        query_genome = forms.ChoiceField(choices=accession_choices)

        def save(self):
            self.reference_genome = self.cleaned_data["reference_genome"]
            self.query_genome = self.cleaned_data["query_genome"]

    return Circos2genomesForm


def make_mummer_form(database_name):

    accession_choices = get_accessions(database_name)

    class Circos2genomesForm(forms.Form):

        reference_genome = forms.ChoiceField(choices=accession_choices)
        query_genome = forms.ChoiceField(choices=accession_choices)

        def save(self):
            self.reference_genome = self.cleaned_data["reference_genome"]
            self.query_genome = self.cleaned_data["query_genome"]

    return Circos2genomesForm



def make_crossplot_form(database_name):


    accession_choices = get_accessions(database_name)

    class CrossplotForm(forms.Form):

        accession = forms.CharField(max_length=100)
        region_size = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
        reference_genome = forms.ChoiceField(choices=accession_choices)
        query_genome = forms.ChoiceField(choices=accession_choices)

        def save(self):
            self.reference_genome = self.cleaned_data["reference_genome"]
            self.query_genome = self.cleaned_data["query_genome"]

    return CrossplotForm


class ConnexionForm(forms.Form):
    username = forms.CharField(label="User Name", max_length=30)
    password = forms.CharField(label="Password", widget=forms.PasswordInput)
    biodatabase = forms.ChoiceField(choices=choices, required = False)


class DBForm(ModelForm):

    accession = forms.CharField(max_length=100)
    region_size = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
    #target_region = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
    '''
    class Meta:
        model=GenDB
        fields = ["ref_genome", "query_genome"]
    '''

