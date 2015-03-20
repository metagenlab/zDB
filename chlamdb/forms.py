#-*- coding: utf-8 -*-
from django import forms
#import models
from manipulate_biosqldb import load_db, get_biodatabase_list
from blah import *
from models import GenDB
from models import Genome
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
        sql ='SELECT bioentry.taxon_id, bioentry.description FROM bioentry ' \
             'inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id ' \
             'where biodatabase.name ="%s"and bioentry.description not like "%%%%plasmid%%%%"' % database_name #
    else:
        sql ='SELECT bioentry.taxon_id, bioentry.description FROM bioentry ' \
             'inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id ' \
             'where biodatabase.name ="%s"' % database_name
    result = server.adaptor.execute_and_fetchall(sql, )
    accession_list = [i for i in result]
    print "acc", accession_list
    accession_choices = []
    if all:
        accession_choices.append(("all", "all"))

    for accession in accession_list:
        accession_choices.append((accession[0], accession[1]))


    import re
    for i, accession in enumerate(accession_choices):
        print i, accession
        description = accession[1]
        description = re.sub(", complete genome\.", "", description)
        description = re.sub(", complete sequence\.", "", description)
        description = re.sub("strain ", "", description)
        description = re.sub("str\. ", "", description)
        description = re.sub(" complete genome sequence\.", "", description)
        description = re.sub(" complete genome\.", "", description)
        description = re.sub(" chromosome", "", description)
        description = re.sub(" DNA", "S.", description)
        accession_choices[i] = (accession[0], description)

    return accession_choices

biodatabases = [ i for i in get_biodatabase_list(server)]

choices = []
for i in biodatabases:
    choices.append((i,i))

choices = tuple(choices)


def make_contact_form(server, database_name):

    accession_choices = get_accessions(database_name)


    class ContactForm(forms.Form):

        accession = forms.CharField(max_length=100)
        #biodatabase = forms.ChoiceField(choices=choices)
        plot_region = forms.NullBooleanField(widget=forms.CheckboxInput())
        region_size = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
        genomes = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "30" }), required = False)
    return ContactForm

def make_plot_form(database_name):


    accession_choices = get_accessions(database_name)


    class PlotForm(forms.Form):

        location_start = forms.CharField(max_length=9, label="sart (bp)", required=False)
        location_stop = forms.CharField(max_length=9, label="end (bp)", required=False)
        genome = forms.ChoiceField(choices=accession_choices, required=True)
    return PlotForm


class BiodatabaseForm(forms.Form):
    biodatabase = forms.ChoiceField(choices=choices, required = False)

    def save(self):
        self.biodatabase = self.cleaned_data["biodatabase"]



def make_circos_form(database_name):

    accession_choices = get_accessions(database_name)
    
    class CircosForm(forms.Form):
        reference = forms.ChoiceField(choices=accession_choices)
        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % len(accession_choices) }), required = False)
        get_region = forms.NullBooleanField(widget=forms.CheckboxInput())
        region = forms.CharField(max_length=100, label="Region start, stop", initial = "1, 8000", required = False)
        
        def save(self):
            self.reference = self.cleaned_data["reference"]
            self.get_region = self.cleaned_data["get_region"]
            self.region = self.cleaned_data["region"]
            

            
    return CircosForm

def make_extract_form(database_name):

    accession_choices = get_accessions(database_name)

    class ExtractForm(forms.Form):
        orthologs_in = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "20" }), required = False, label="")
        no_orthologs_in = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "20" }), required = False, label="")



    return ExtractForm


class SearchForm(forms.Form):

    SEARCH_CHOICES = (('gene', 'gene'), ('product','product'))

    search_type = forms.ChoiceField(choices=SEARCH_CHOICES)
    search_term = forms.CharField(max_length=100)
    #biodatabase = forms.ChoiceField(choices=choices)





class BlastForm(forms.Form):

    blast_input = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 20}))
    #biodatabase = forms.ChoiceField(choices=choices)



def make_motif_form(database_name):

    accession_choices = get_accessions(database_name, all=True)

    class MotifForm(forms.Form):
        '''

         choix fuzznuc / fuzzpro
        '''
        DATA_CHOICES = (('nucleotide', 'nucleotide'), ('protein','protein'))
        search = forms.ChoiceField(choices=DATA_CHOICES)
        n_missmatch = forms.CharField(max_length=3, label="N. mismatches", required=False, initial=0)
        input_file = forms.FileField(
            label='Select a file',
            help_text='max. 42 megabytes', required=False
        )
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
    class Meta:
        model=GenDB
        fields = ["ref_genome", "query_genome"]


