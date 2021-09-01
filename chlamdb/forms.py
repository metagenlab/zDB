#-*- coding: utf-8 -*-
from django import forms
#import models
from chlamdb.biosqldb.manipulate_biosqldb import load_db, get_biodatabase_list
#from blah import *
#from models import GenDB
#from models import Genome
from django.forms import ModelForm
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, Row, Column, MultiField, Div, Fieldset
from crispy_forms.bootstrap import AppendedText
from django.conf import settings

choices = []

class GenerateRandomUserForm(forms.Form):
    total_user = forms.IntegerField(
        label='Number of users',
        required=True,
    )


def get_accessions(db, all=False, plasmid=False):
    result            = db.get_genomes_description(lst_plasmids=True)
    accession_choices = []
    index             = 0
    reverse_index     = []

    # cannot use taxid because plasmids have the same
    # taxids as the chromosome
    for taxid, data in result.iterrows():
        accession_choices.append((index, data.description))
        index += 1
        reverse_index.append((taxid, False))
        if plasmid and data.has_plasmid==1:
            accession_choices.append((index, data.description + " plasmid"))
            reverse_index.append((taxid, True))
            index += 1
            is_plasmid = True

    if all:
        accession_choices = [["all", "all"]] + accession_choices
    return accession_choices , reverse_index


def get_accessions_BLAST(db, all=False, plasmid=False):
    result            = db.get_genomes_description(lst_plasmids=True)
    accession_choices = []
    index          = 1
    reverse_index     = []

    # cannot use taxid because plasmids have the same
    # taxids as the chromosome
    for taxid, data in result.iterrows():
        accession_choices.append((index, data.description))
        reverse_index.append((taxid, False))
        try:
            if plasmid and data.has_plasmid==1:
                accession_choices.append((str(index) + " plasmid", data.description + " plasmid"))
                reverse_index.append((taxid, True))
                is_plasmid = True
                index += 1
        except:
            index += 1

    if all:
        accession_choices = [["all", "all"]] + accession_choices
    return accession_choices #, reverse_index  #understand why there is this 'reverse_index' in the original def (error in def blast)


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
        choices = (("yes","yes"),("no", "best hit only"))
        accession = forms.CharField(max_length=100, label="Protein accession (e.g. CT_015)")
        #location_start = forms.CharField(max_length=9, label="sart (bp)", required=False)
        #location_stop = forms.CharField(max_length=9, label="end (bp)", required=False)
        region_size = forms.CharField(max_length=5, label="Region size (bp)", initial = 8000, required = False)
        genomes = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'1', "class":"selectpicker", "data-live-search":"true", "multiple data-max-options":"8", "multiple data-actions-box":"true"}), required = False)
        all_homologs = forms.ChoiceField(choices=choices, initial="no")

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            #self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            #self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset("Compare genomes",
                                                Column(
                                                Row(Column("accession", css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                                                    Column("region_size", css_class='form-group col-lg-6 col-md-6 col-sm-12')),
                                                    Row('genomes', style="padding-left: 15px"),
                                                Column(Row('all_homologs'),
                                                Submit('submit', 'Compare'), css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                                                css_class="col-lg-8 col-md-8 col-sm-12")
                                                )
                                        )

            super(PlotForm, self).__init__(*args, **kwargs)

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


        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % (20) }), required = False)
        search_type = forms.ChoiceField(choices=SEARCH_CHOICES)
        search_term = forms.CharField(max_length=100)


    return InterproForm


def make_metabo_from(db, add_box=False):

    accession_choices, rev_index = get_accessions(db)

    class MetaboForm(forms.Form):
        targets = forms.MultipleChoiceField(choices=accession_choices,
                widget=forms.SelectMultiple(attrs={'size':'%s' % (20), "class":"selectpicker", "data-live-search":"true", "multiple data-actions-box":"true"}),
                required = False)

        if add_box:
            input_box = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 10}))

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-1 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-4 col-md-6 col-sm-6'
            if not add_box:
                self.helper.layout = Layout(
                                            Fieldset("Compare genomes",
                                                     Column(
                                                           Row('targets'),
                                                           Submit('submit', 'Submit'),
                                                           css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                                                    )
                                            )

            else:
                self.helper.layout = Layout(
                                            Fieldset("Compare genomes",
                                                     Column(
                                                           Row('targets'),
                                                           Row('input_box'),
                                                           Submit('submit', 'Submit'),
                                                           css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                                                    )
                                            )
            super(MetaboForm, self).__init__(*args, **kwargs)

        def get_choices(self):
            targets = self.cleaned_data["targets"]
            taxids = []
            for index in targets:
                taxid, _ = rev_index[int(index)]
                taxids.append(taxid)
            return taxids

    return MetaboForm


def make_venn_from(db, plasmid=False, label="Orthologs", limit=None):

    accession_choices, rev_index = get_accessions(db, plasmid=plasmid)

    class VennForm(forms.Form):
        if limit is None:
            targets = forms.MultipleChoiceField(choices=accession_choices,
                    widget=forms.SelectMultiple(attrs={'size':'1', "class":"selectpicker", "data-live-search":"true", "multiple data-actions-box":"true"}), required = True)
        else:
            targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'1', "class":"selectpicker", "data-live-search":"true", "multiple data-max-options":"%s" % limit ,"multiple data-actions-box":"true"}), required = True)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-1 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-4 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset("Compare genomes",
                                                 Column(
                                                       Row('targets'),
                                                       Submit('submit', 'Compare %s' % label,  style="margin-top:15px" ),
                                                       css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                                                )
                                        )

            super(VennForm, self).__init__(*args, **kwargs)

        def get_taxids(self):
            indices = self.cleaned_data["targets"]
            taxids  = []
            for index in indices:
                taxid, _ = rev_index[int(index)]
                taxids.append(taxid)
            return taxids

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

def make_blast_form(biodb):

    accession_choices =  get_accessions_BLAST(biodb, plasmid=True, all=True)

    print(accession_choices)
    class BlastForm(forms.Form):
        blast = forms.ChoiceField(choices=[("blastn_ffn", "blastn_ffn"),
                                           ("blastn_fna", "blastn_fna"),
                                           ("blastp", "blastp"),
                                           ("blastx", "blastx"),
                                           ("tblastn", "tblastn")])

        max_number_of_hits = forms.ChoiceField(choices=[("10", "10"),
                                           ("5", "5"),
                                           ("20", "20"),
                                           ("30", "30"),
                                           ("all", "all")])
        evalue= forms.CharField(widget=forms.TextInput({'placeholder': '10'}))

        target = forms.ChoiceField(choices=accession_choices, widget=forms.Select(attrs={"class":"selectpicker", "data-live-search":"true", }))
        blast_input = forms.CharField(widget=forms.Textarea(attrs={'cols': 50, 'rows': 5}))



def make_extract_region_form(database_name):

    accession_choices = get_accessions(database_name, plasmid=True)

    extraction_choices = [['annotation', 'annotation'],['sequence', 'sequence'],['sequence_trans', 'translation']]

    class ExtractRegionForm(forms.Form):



        genome = forms.ChoiceField(choices=accession_choices, required = True)
        region = forms.CharField(max_length=100, label="Region start, stop", initial = "1, 8000", required = True)
        extract = forms.ChoiceField(choices=extraction_choices, required = True)
        #get_annotation = forms.NullBooleanField(widget=forms.CheckboxInput())
        #get_sequence = forms.NullBooleanField(widget=forms.CheckboxInput())

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            #self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            #self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset(
                                                Row('Extract region'),
                                                Row('genome'),
                                                Row('region'),
                                                Row('extract'),
                                                Submit('submit', 'Submit'),
                                                css_class="col-lg-5 col-md-6 col-sm-6")
                                        )

            super(ExtractRegionForm, self).__init__(*args, **kwargs)

    return ExtractRegionForm


def make_priam_form(database_name):

    accession_choices = get_accessions(database_name, plasmid=True)

    class PriamForm(forms.Form):
        genome = forms.ChoiceField(choices=accession_choices)

    return PriamForm


def make_species_curation_form(database_name, species_id):
    
    server, db = load_db(database_name)
    
    sql = 'select phylum, `order`, family, genus, species from species_curated_taxonomy where species_id=%s;' % (species_id)
    print(sql)
    data = server.adaptor.execute_and_fetchall(sql,)[0]
    print(data)
    class TaxonomyCurationForm(forms.Form):
        phylum = forms.CharField(max_length=600, required = True, initial=data[0])
        order = forms.CharField(max_length=600, required = True, initial=data[1])
        family = forms.CharField(max_length=600, required = True, initial=data[2])
        genus = forms.CharField(max_length=600, required = True, initial=data[3])
        species = forms.CharField(max_length=600, required = True, initial=data[4])

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            #self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            #self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset("Taxonomy curation",
                                                Column(
                                                    Row('phylum', style="padding-left:15px"),
                                                    Row('order', style="padding-left:15px"),
                                                    Row('family', style="padding-left:15px"),
                                                    Row('genus', style="padding-left:15px"),
                                                    Row('species', style="padding-left:15px"),
                                                Submit('submit', 'Save'), css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                                                css_class="col-lg-8 col-md-8 col-sm-12")
                                        )

            super(TaxonomyCurationForm, self).__init__(*args, **kwargs)

    return TaxonomyCurationForm



def make_circos_form(database_name):

    accession_choices = get_accessions(database_name)

    class CircosForm(forms.Form):
        circos_reference = forms.ChoiceField(choices=accession_choices)
        targets = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'1', "class":"selectpicker", "data-live-search":"true", "multiple data-max-options":"8",}), required=False)
        #get_region = forms.NullBooleanField(widget=forms.CheckboxInput())
        #region = forms.CharField(max_length=100, label="Region start, stop", initial = "1, 8000", required = False)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            #self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            #self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset(
                                                Row("Circos"),
                                                Row('circos_reference'),
                                                Row('targets'),
                                                Submit('submit_circos', 'Submit',  style="padding-left:15px"),
                                                css_class="col-lg-5 col-md-6 col-sm-6")
                                        )

            super(CircosForm, self).__init__(*args, **kwargs)

        def save(self):
            self.reference = self.cleaned_data["reference"]
            self.get_region = self.cleaned_data["get_region"]
            self.region = self.cleaned_data["region"]


    return CircosForm


def make_genome_selection_form(database_name):

    accession_choices = get_accessions(database_name)

    class GenomeForm(forms.Form):
        genome = forms.ChoiceField(choices=accession_choices)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            #self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            #self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset(
                                                Row("Genome"),
                                                Row('genome'),
                                                Submit('submit', 'Submit'),
                                                css_class="col-lg-5 col-md-6 col-sm-6")
                                        )

            super(GenomeForm, self).__init__(*args, **kwargs)


    return GenomeForm


def make_kegg_form(database_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql_pathways = 'select description,description from enzyme_locus2ko t1 ' \
                   ' inner join  enzyme_pathway2ko t2 ' \
                   ' on t1.ko_id = t2.ko_id ' \
                   ' inner join  enzyme_kegg_pathway t3 ' \
                   ' on t3.pathway_id=t2.pathway_id group by description;'

    pathway_choices = server.adaptor.execute_and_fetchall(sql_pathways,)

    sql_modules = 'select description,description from  enzyme_locus2ko t1 ' \
                  ' inner join  enzyme_module2ko t2 on t1.ko_id = t2.ko_id ' \
                  ' inner join  enzyme_kegg_module t3 on t3.module_id=t2.module_id group by description;'

    module_choices = server.adaptor.execute_and_fetchall(sql_modules,)

    class KeggForm(forms.Form):

        pathway_choice = forms.MultipleChoiceField(label='', choices=pathway_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "17" }), required = False)
        module_choice = forms.MultipleChoiceField(choices=module_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "17" }), required = False, label="")

    return KeggForm


def make_extract_form(db, plasmid=False, label="Orthologs"):
    accession_choices, rev_index = get_accessions(db, plasmid=plasmid)

    class ExtractForm(forms.Form):
        checkbox_accessions = forms.BooleanField(required = False,
                label="Distinguish plasmids from chromosomes")
        checkbox_single_copy = forms.BooleanField(required = False,
                label="Only consider single copy %s" % label)

        orthologs_in = forms.MultipleChoiceField(label=f"{label} conserved in",
                choices=accession_choices,
                widget=forms.SelectMultiple(attrs={'size':'%s' % "17", "class":"selectpicker", "data-live-search":"true"}), required = True)
        no_orthologs_in = forms.MultipleChoiceField(label="%s absent from (optional)" % label,
                choices=accession_choices,
                widget=forms.SelectMultiple(attrs={'size':'%s' % "17", "class":"selectpicker remove-example", "data-live-search":"true"}), required = False)

        new_choices = [['None', 'None']] + accession_choices

        frequency_choices = ((i, i) for i in range(len(accession_choices)))
        frequency = forms.ChoiceField(choices=frequency_choices,
                label='Missing data (optional)', required = False)

        def extract_choices(self, indices):
            keep_plasmids = self.cleaned_data["checkbox_accessions"]
            taxids = []
            plasmids = None
            if keep_plasmids:
                plasmids = []

            for index in indices:
                taxid, is_plasmid = rev_index[index]
                if keep_plasmids and is_plasmid:
                    plasmids.append(taxid)
                elif not keep_plasmids or not is_plasmid:
                    taxids.append(taxid)
            return taxids, plasmids


        def get_n_missing(self):
            return int(self.cleaned_data["frequency"])

        def get_include_choices(self):
            return self.extract_choices((int(i) for i in self.cleaned_data["orthologs_in"]))

        def get_exclude_choices(self):
            return self.extract_choices((int(i) for i in self.cleaned_data["no_orthologs_in"]))

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            #self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            #self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                    Fieldset("Compare genomes",
                            Column(
                                Row('checkbox_accessions', style="padding-left:15px"),
                                Row('checkbox_single_copy', style="padding-left:15px"),
                            Row(Column("orthologs_in", css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                                Column("no_orthologs_in", css_class='form-group col-lg-6 col-md-6 col-sm-12')),
                            Column(Row('frequency'),
                            Submit('submit', 'Compare %s' % label,   style="margin-top:15px"), css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                            css_class="col-lg-8 col-md-8 col-sm-12")
                            )
                    )

            super(ExtractForm, self).__init__(*args, **kwargs)
    return ExtractForm


def locus_int_form(database_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql = 'select distinct category from custom_tables_annot_table;'
    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [(i[0],i[0]) for i in categories]
    CHOICES.append(("all","all"))

    show_identity_CHOICES = (
        ('ientity', 'Show identity percentages'),
        ('color', 'Show coloured heatmap'),
    )


    class IntCatChoice(forms.Form):
        category = forms.ChoiceField(choices=CHOICES)
        #Show_identity_values = forms.ChoiceField(choices=show_identity_CHOICES,
        #                                        widget=forms.RadioSelect)
    return IntCatChoice


def hmm_sets_form_circos(database_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql = 'select distinct name from hmm.hmm_sets;' #% database_name
    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [(i[0],i[0]) for i in categories]
    CHOICES.append(("all","all"))

    accession_choices = get_accessions(database_name)

    class HmmSetChoice(forms.Form):
        hmm_set = forms.ChoiceField(choices=CHOICES)
        genome = forms.ChoiceField(choices=accession_choices)
        score_cutoff = forms.CharField(max_length=3, label="Bitscoire cutoff", initial = 10, required = False)
        query_coverage_cutoff = forms.CharField(max_length=3, label="Query coverage cutoff", initial = 0.5, required = False)

    return HmmSetChoice



def hmm_sets_form(database_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql = 'select distinct name from hmm.hmm_sets;' #% database_name
    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [(i[0],i[0]) for i in categories]
    CHOICES.append(("all","all"))

    class HmmSetChoice(forms.Form):
        hmm_set = forms.ChoiceField(choices=CHOICES)

        score_cutoff = forms.CharField(max_length=3, label="Bitscoire cutoff", initial = 10, required = False)
        query_coverage_cutoff = forms.CharField(max_length=3, label="Query coverage cutoff", initial = 0.5, required = False)

    return HmmSetChoice


def transporters_superfam_form(database_name, show_taxon=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql = 'select distinct description from transporters.transporter_table t1 ' \
          'inner join transporters.tc_table t2 on t1.superfamily=t2.tc_id;' #% database_name

    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [("all","all")]
    CHOICES += [(i[0],i[0]) for i in categories]
    if show_taxon:
        accession_choices = get_accessions(database_name)
    class TransporterSuperfamilyChoice(forms.Form):

        if show_taxon:
            genome = forms.ChoiceField(choices=accession_choices)

        transporter_superfamily = forms.MultipleChoiceField(choices=CHOICES, widget=forms.SelectMultiple(attrs={'size':'20' }), required = False)

        score_cutoff = forms.CharField(max_length=3, label="Bitscore cutoff", initial = 10, required = False)
        evalue_cutoff = forms.CharField(max_length=10, label="Evalue cutoff", initial = 0.005, required = False)
        query_coverage_cutoff = forms.CharField(max_length=3, label="Query coverage cutoff", initial = 0.5, required = False)
        hit_coverage_cutoff = forms.CharField(max_length=3, label="Query coverage cutoff", initial = 0.5, required = False)

    return TransporterSuperfamilyChoice


def blast_sets_form(database_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)

    sql = 'select distinct name from blast.blast_sets;' #% database_name
    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [(i[0],i[0]) for i in categories]
    CHOICES.append(("all","all"))

    class BlastSetChoice(forms.Form):
        blast_set = forms.MultipleChoiceField(choices=CHOICES, widget=forms.SelectMultiple(attrs={'size':'20' }), required = False)

        score_cutoff = forms.CharField(max_length=3, label="Bitscore cutoff", initial = 10, required = False)
        query_coverage_cutoff = forms.CharField(max_length=3, label="Query coverage cutoff", initial = 0.5, required = False)
        hit_coverage_cutoff = forms.CharField(max_length=3, label="Query coverage cutoff", initial = 0.5, required = False)

    return BlastSetChoice

def make_locus2network_form(database_name):
    accessions = get_accessions(database_name)

    class Network(forms.Form):
        genome = forms.ChoiceField(choices=accessions)
        locus_list = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 10}), required = False)

    return Network

def make_pairwiseid_form(database_name):

    plot_CHOICES = (('blast_identity', 'blast_identity'), ('msa_identity','msa_identity'), ('score','score'))

    accessions = get_accessions(database_name)
    accessions2 = [['None', 'None']] + accessions

    class PairwiseID(forms.Form):
        plot = forms.ChoiceField(choices=plot_CHOICES)
        genome_1 = forms.ChoiceField(choices=accessions)
        genome_2 = forms.ChoiceField(choices=accessions)
        genome_3 = forms.ChoiceField(choices=accessions2)
        genome_4 = forms.ChoiceField(choices=accessions2)
        genome_5 = forms.ChoiceField(choices=accessions2)
        genome_6 = forms.ChoiceField(choices=accessions2)

    return PairwiseID


def make_pairwiseCDS_length_form(database_name):

    accessions = get_accessions(database_name)
    accessions2 = [['None', 'None']] + accessions

    class PairwiseID(forms.Form):
        genome_1 = forms.ChoiceField(choices=accessions)
        genome_2 = forms.ChoiceField(choices=accessions)
        genome_3 = forms.ChoiceField(choices=accessions2)
        genome_4 = forms.ChoiceField(choices=accessions2)

    return PairwiseID

def heatmap_form(database_name):

    plot_CHOICES = (('blast_identity', 'blast_identity'), ('core_msa_identity','core_msa_identity'),
                    ('n_RBBH','n_RBBH'), ('n_shared_orthogroups','n_shared_orthogroups'))

    accessions = get_accessions(database_name)

    class Heatmap(forms.Form):
        plot = forms.ChoiceField(choices=plot_CHOICES)
        targets = forms.MultipleChoiceField(choices=accessions, widget=forms.SelectMultiple(attrs={'size':'20' }), required = False)

    return Heatmap

def make_module_overview_form(db, sub_sub_cat=False):

    if sub_sub_cat:
        categories = db.get_module_sub_categories()
    else:
        categories = db.get_module_categories()

    CHOICES = [(cat_id, cat) for cat_id, cat in categories]
    class ModuleCatChoice(forms.Form):
        category = forms.ChoiceField(choices=CHOICES)

    return ModuleCatChoice


def make_pathway_overview_form(database_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(database_name)
    sql = 'select distinct pathway_category from enzyme_kegg_pathway where pathway_category_short not in ("drug","disease","organismal") ' \
          ' order by pathway_category;'

    categories = server.adaptor.execute_and_fetchall(sql,)
    CHOICES = [(i[0],i[0]) for i in categories]


    class PathwayCatChoice(forms.Form):
        category = forms.ChoiceField(choices=CHOICES)
    return PathwayCatChoice


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


def make_blastnr_best_non_top_phylum_form(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    accession_choices = get_accessions(biodb, all=False)

    class Blastnr_top(forms.Form):
        accessions = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "15" }), required = False)

        CHOICES=[('all','all'),
         ('specific','species specific')]
        selection = forms.ChoiceField(choices=CHOICES)
    return Blastnr_top


def make_blastnr_form(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)



    sql ='SELECT bioentry.bioentry_id, bioentry.description FROM bioentry ' \
             'inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id ' \
             'where biodatabase.name ="%s"' \
             'order by bioentry.description' % biodb
    result = server.adaptor.execute_and_fetchall(sql, )
    accession_list = [i for i in result]
    accession_choices = []

    for accession in accession_list:
        accession_choices.append((accession[0], accession[1]))

    sql = 'show columns from blastnr_blastnr_taxonomy;'
    ranks = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    rank_choices = []

    sql2 = 'select max(hit_number) from blastnr_blastnr ;' % biodb

    max_nr_hits = server.adaptor.execute_and_fetchall(sql2,)[0][0]

    for rank in ranks:
        rank_choices.append((rank, rank))

    class Blastnr_top(forms.Form):
        accession = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "15" , 'class':"selectpicker", "data-width":"200px",}), required = False)
        rank = forms.ChoiceField(choices=rank_choices)
        CHOICES=[('BBH','BBH'),
         ('Majority','Majority Rule')]
        type = forms.ChoiceField(choices=CHOICES)
        top_number = forms.CharField(max_length=3, label="Majority over top n hits", initial = max_nr_hits, required = False)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset(
                                                Row("BLAST taxonomy"),
                                                Row('accession'),
                                                Row('rank'),
                                                Row('type'),
                                                Row('top_number'),
                                                Submit('submit', 'Submit'),
                                                css_class="col-lg-5 col-md-6 col-sm-6")
                                        )

            super(Blastnr_top, self).__init__(*args, **kwargs)
    return Blastnr_top

def make_interpro_taxonomy(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    accession_choices = get_accessions(biodb,all=True)

    # p_eukaryote | p_archae | p_virus
    TAXO_CHOICES = (('p_eukaryote', 'eukaryote'), ('p_archae','archae'), ('p_virus','virus'))

    class Interpro_taxonomy(forms.Form):

        target_taxons = forms.MultipleChoiceField(choices=accession_choices, widget=forms.SelectMultiple(attrs={'size':'%s' % "15" }), required = False)
        percentage_cutoff = forms.CharField(max_length=3, label="Percentage Cutoff", initial=90, required=True)
        kingdom = forms.ChoiceField(choices=TAXO_CHOICES)

    return Interpro_taxonomy


class BlastProfileForm(forms.Form):

    fasta_file = forms.FileField(label='Select a file')

    blast = forms.ChoiceField(choices=[("blastn_fna", "blastn_fna"),
                                       ("blastp", "blastp"),
                                       ("tblastn", "tblastn")], required = False)




def make_blast_form(biodb):

    accession_choices =  get_accessions_BLAST(biodb, plasmid=False, all=True)

    print(accession_choices)
    class BlastForm(forms.Form):
        blast = forms.ChoiceField(choices=[("blastn_ffn", "blastn_ffn"),
                                           ("blastn_fna", "blastn_fna"),
                                           ("blastp", "blastp"),
                                           ("blastx", "blastx"),
                                           ("tblastn", "tblastn")])

        max_number_of_hits = forms.ChoiceField(choices=[("10", "10"),
                                           ("5", "5"),
                                           ("20", "20"),
                                           ("30", "30"),
                                           ("all", "all")])
        evalue= forms.CharField(widget=forms.TextInput({'placeholder': '10'}))

        target = forms.ChoiceField(choices=accession_choices, widget=forms.Select(attrs={"class":"selectpicker", "data-live-search":"true"}))
        blast_input = forms.CharField(widget=forms.Textarea(attrs={'cols': 50, 'rows': 5}))


        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-4 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-6 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                                        Fieldset(
                                                Row("BLAST"),
                                                Row('target'),
                                                Row('blast_input'),
                                                css_class="col-lg-5 col-md-6 col-sm-6")
                                        )

            super(BlastForm, self).__init__(*args, **kwargs)

    return BlastForm

def make_comment_from(biodb, locus_tag):
    import chlamdb.manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select * from manual_annotation where locus_tag="%s"' % locus_tag
    data = server.adaptor.execute_and_fetchall(sql,)
    if len(data) == 0:
        description = '-'
    else:
        description = data[0][1]

    class CommentForm(forms.Form):

        locus = forms.CharField(max_length=200, initial="%s" % locus_tag, required = True)
        comment = forms.CharField(widget=forms.Textarea(attrs={'cols': 80, 'rows': 3}), initial="%s" % description, required = True)


    return CommentForm


class LocusInt(forms.Form):
    category = forms.CharField(max_length=600, required = True)
    gene = forms.CharField(max_length=200, required = True)
    locus_tag = forms.CharField(max_length=400, required = True)
    description = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 2,'style': 'float:left;position:relative;left:25px;font-size: 10px;width:400px;'}) , required = True)
    reference = forms.CharField(widget=forms.Textarea(attrs={'cols': 60, 'rows': 2, 'style': 'float:left;position:relative;left:-5px;font-size: 10px;width:400px;'}), required = True)



class AnnotForm(forms.Form):
    orthogroups = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 10}))

def get_LocusAnnotForm(database_name):

    accession_choices = get_accessions(database_name)

    class LocusAnnotForm(forms.Form):
        locus_list = forms.CharField(widget=forms.Textarea(attrs={'cols': 10, 'rows': 10}))
        circos_target = forms.ChoiceField(choices=accession_choices, label='Reference genome')
    return LocusAnnotForm

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
