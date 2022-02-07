#-*- coding: utf-8 -*-

from django import forms

from crispy_forms.helper import FormHelper
from crispy_forms.layout import Layout, Submit, Row, Column, Fieldset

choices = []


def get_accessions(db, all=False, plasmid=False):
    result            = db.get_genomes_description(lst_plasmids=plasmid)
    accession_choices = []
    index             = 0
    reverse_index     = []

    # cannot use taxid because plasmids have the same
    # taxids as the chromosome
    for taxid, data in result.iterrows():
        accession_choices.append((index, data.description))
        index += 1

        if plasmid:
            reverse_index.append((taxid, False))
        else:
            reverse_index.append(taxid)
        if plasmid and data.has_plasmid==1:
            accession_choices.append((index, data.description + " plasmid"))
            reverse_index.append((taxid, True))
            index += 1
            is_plasmid = True

    if all:
        accession_choices = [["all", "all"]] + accession_choices
    return accession_choices, reverse_index


def make_plot_form(db):
    accession_choices, reverse_index = get_accessions(db)

    class PlotForm(forms.Form):
        choices = (("yes", "all homologs"),("no", "best hits only"))
        accession = forms.CharField(max_length=100,
                label="Protein accession (e.g. CT_015)", required=True)
        region_size = forms.CharField(max_length=5,
                label="Region size (bp)", initial = 8000, required = False)
        genomes = forms.MultipleChoiceField(choices=accession_choices,
                widget=forms.SelectMultiple(attrs={'size':'1', "class":"selectpicker",
                    "data-live-search":"true", "multiple data-max-options":"8",
                    "multiple data-actions-box":"true"}),
                required = False)
        all_homologs = forms.ChoiceField(choices=choices, initial="no", label="Homologs")

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.layout = Layout(
                Fieldset("Compare genomes",
                    Column(
                    Row(Column("accession", css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                        Column("region_size", css_class='form-group col-lg-6 col-md-6 col-sm-12')),
                        Row('genomes', style="padding-left: 15px"),
                    Column(Row('all_homologs', css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                    Submit('submit', 'Compare'), css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                    css_class="col-lg-8 col-md-8 col-sm-12"))
            )

        def get_accession(self):
            return self.cleaned_data["accession"]

        def get_all_homologs(self):
            return self.cleaned_data.get("all_homologs", "") == "yes"

        def get_region_size(self):
            return int(self.cleaned_data.get("region_size", 8000))

        def get_genomes(self):
            indices = self.cleaned_data["genomes"]
            taxids  = []
            for index in indices:
                taxid = reverse_index[int(index)]
                taxids.append(taxid)
            return taxids

    return PlotForm


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
                taxid = rev_index[int(index)]
                taxids.append(taxid)
            return taxids

    return MetaboForm


def make_venn_from(db, plasmid=False, label="Orthologs", limit=None,
        action="venn_orthogroup"):

    accession_choices, rev_index = get_accessions(db, plasmid=plasmid)

    class VennForm(forms.Form):
        attrs = {'size':'1', "class":"selectpicker",
                "data-live-search":"true",
                "data-actions-box":"true"}
        if not limit is None:
            attrs["data-max-options"] = f"{limit}"
        targets = forms.MultipleChoiceField(choices=accession_choices,
                    widget=forms.SelectMultiple(attrs=attrs), required = True)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-1 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-4 col-md-6 col-sm-6'
            self.helper.form_action = action
            self.helper.layout = Layout(
                Fieldset("Compare genomes",
                     Column(
                           Row('targets'),
                           Submit('submit', 'Compare %s' % label,  style="margin-top:15px" ),
                           css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                    )
                )

        def get_taxids(self):
            indices = self.cleaned_data["targets"]
            taxids  = []
            for index in indices:
                taxid = rev_index[int(index)]
                taxids.append(taxid)
            return taxids

        def clean_venn(self):
            value = self.cleaned_data['targets']
            if len(value) > 6:
                raise forms.ValidationError("You can't select more than 6 items.")
            return value

    return VennForm


def make_circos_form(database_name):

    accession_choices, reverse_index = get_accessions(database_name)

    class CircosForm(forms.Form):
        circos_reference = forms.ChoiceField(choices=accession_choices,
                widget=forms.Select(attrs={"class":"selectpicker", "data-live-search":"true"}))
        targets = forms.MultipleChoiceField(choices=accession_choices,
                widget=forms.SelectMultiple(attrs={'size':'1', "class":"selectpicker", "data-live-search":"true", "multiple data-max-options":"10"}), required=False)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            
            self.helper.layout = Layout(
                Fieldset(
                        Row("Circos"),
                        Row('circos_reference'),
                        Row('targets', style="margin-top:1em"),
                        Submit('submit_circos', 'Submit',
                            style="padding-left:15px; margin-top:15px; margin-bottom:15px "),
                        css_class="col-lg-5 col-md-6 col-sm-6")
                )
            super(CircosForm, self).__init__(*args, **kwargs)

        def save(self):
            self.reference = self.cleaned_data["reference"]
            self.get_region = self.cleaned_data["get_region"]
            self.region = self.cleaned_data["region"]

        def get_target_taxids(self):
            indices = self.cleaned_data["targets"]
            taxids  = []
            for index in indices:
                taxid = reverse_index[int(index)]
                taxids.append(taxid)
            return taxids

        def get_ref_taxid(self):
            indice = self.cleaned_data["circos_reference"]
            taxid = reverse_index[int(indice)]
            return taxid

    return CircosForm


def make_extract_form(db, action, plasmid=False, label="Orthologs"):
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
            self.helper.form_action = action
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


def make_module_overview_form(db, sub_sub_cat=False):

    if sub_sub_cat:
        categories = db.get_module_sub_categories()
    else:
        categories = db.get_module_categories()

    CHOICES = [(cat_id, cat) for cat_id, cat in categories]
    class ModuleCatChoice(forms.Form):
        category = forms.ChoiceField(choices=CHOICES)

    return ModuleCatChoice


def make_single_genome_form(db):
    accession_choices, rev_index = get_accessions(db)

    class SingleGenomeForm(forms.Form):
        genome = forms.ChoiceField(choices=accession_choices)
        def get_genome(self):
            target = self.cleaned_data["genome"]
            return rev_index[int(target)]

    return SingleGenomeForm


def make_blast_form(biodb):
    accession_choices, rev_index =  get_accessions(biodb, all=True)

    class BlastForm(forms.Form):
        DEFAULT_E_VALUE = 10
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
        target = forms.ChoiceField(choices=accession_choices,
                widget=forms.Select(attrs={"class":"selectpicker", "data-live-search":"true"}))
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

        def get_target(self):
            target = self.cleaned_data["target"]
            if target=="all":
                return target
            return rev_index[int(target)]

    return BlastForm
