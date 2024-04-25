# -*- coding: utf-8 -*-
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from crispy_forms.helper import FormHelper
from crispy_forms.layout import Column, Fieldset, Layout, Row, Submit
from django import forms
from django.core.exceptions import ValidationError
from django.core.validators import MaxLengthValidator, MinLengthValidator


def get_accessions(db, all=False, plasmid=False):
    result = db.get_genomes_description(lst_plasmids=plasmid)
    accession_choices = []
    index = 0
    reverse_index = []

    # cannot use taxid because plasmids have the same
    # taxids as the chromosome
    for taxid, data in result.iterrows():
        accession_choices.append((index, data.description))
        index += 1

        if plasmid:
            reverse_index.append((taxid, False))
        else:
            reverse_index.append(taxid)
        if plasmid and data.has_plasmid == 1:
            accession_choices.append((index, data.description + " plasmid"))
            reverse_index.append((taxid, True))
            index += 1

    if all:
        accession_choices = [["all", "all"]] + accession_choices
    return accession_choices, reverse_index


def make_plot_form(db):
    import pandas as pd
    accession_choices, reverse_index = get_accessions(db)

    # selct random locus present in at least 50% of genomes
    og_df = pd.DataFrame(db.get_all_orthogroups())
    og_df.columns = ["og", "size"]
    try:
        random_group = og_df.query(
            f"size >= {len(accession_choices) - 1}").iloc[-1]
        locus_list = db.get_genes_from_og(
            [str(random_group.og)], taxon_ids=None, terms=["locus_tag"])
        locus = locus_list.locus_tag.to_list()[0]
    except Exception:
        locus = 'n/a'

    class PlotForm(forms.Form):
        choices = (("yes", "all homologs"), ("no", "best hits only"))
        accession = forms.CharField(max_length=100,
                                    required=True,
                                    label=f"locus_tag (e.g. {locus})")
        region_size = forms.IntegerField(min_value=5000,
                                         max_value=100000,
                                         label="Region size (bp)",
                                         initial=8000,
                                         required=True)
        genomes = forms.MultipleChoiceField(choices=accession_choices,
                                            widget=forms.SelectMultiple(attrs={
                                                'size': '1',
                                                "class": "selectpicker",
                                                "data-live-search": "true",
                                                "multiple data-max-options": "8",
                                                "multiple data-actions-box": "true"}),
                                            required=False)
        all_homologs = forms.ChoiceField(
            choices=choices, initial="no", label="Homologs")

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()

            self.helper.form_method = 'post'
            self.helper.layout = Layout(Fieldset(
                "",
                Column(
                    Row(Column(
                            "accession",
                            css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                        Column(
                            "region_size",
                            css_class='form-group col-lg-6 col-md-6 col-sm-12')
                        ),
                    Row('genomes', style="padding: 15px"),
                    Column(
                        Row('all_homologs',
                            css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                        Submit('submit', 'Compare plot regions'),
                        css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                    css_class="col-lg-8 col-md-8 col-sm-12")
                )
            )
            self.db = db

        def clean_accession(self):
            accession = self.cleaned_data["accession"]
            prot_info = db.get_proteins_info(
                ids=[accession], search_on="locus_tag", as_df=True)
            if prot_info.empty:
                raise ValidationError("Accession not found", code="invalid")
            # Accession is now the sequence id
            return int(prot_info.index[0])

        def get_seqid(self):
            return self.cleaned_data["accession"]

        def get_all_homologs(self):
            return self.cleaned_data.get("all_homologs", "") == "yes"

        def get_region_size(self):
            return int(self.cleaned_data.get("region_size", 8000))

        def get_genomes(self):
            indices = self.cleaned_data["genomes"]
            taxids = []
            for index in indices:
                taxid = reverse_index[int(index)]
                taxids.append(taxid)
            return taxids

    return PlotForm


def make_metabo_from(db, type_choices=None):

    accession_choices, rev_index = get_accessions(db)

    class MetaboForm(forms.Form):
        targets = forms.MultipleChoiceField(
            choices=accession_choices,
            widget=forms.SelectMultiple(
                attrs={'size': '%s' % (20),
                       "class": "selectpicker",
                       "data-live-search": "true",
                       "multiple data-actions-box": "true"}
                ),
            required=True
        )

        if type_choices:
            comp_type = forms.ChoiceField(
                choices=type_choices,
                required=True
            )

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-1 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-4 col-md-6 col-sm-6'
            rows = [Row('targets')]
            if type_choices:
                rows.append(Row('comp_type'))

            self.helper.layout = Layout(

                Fieldset("",
                         Column(
                             *rows,
                             Submit('submit', 'Submit', style="margin-top:15px"),
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
                   limit_type="upper", action=""):

    accession_choices, rev_index = get_accessions(db, plasmid=plasmid)

    class VennForm(forms.Form):
        attrs = {'size': '1', "class": "selectpicker",
                 "data-live-search": "true",
                 "data-actions-box": "true"}
        help_text = ""
        targets_validators = []
        if limit is not None and limit_type == "upper":
            attrs["data-max-options"] = f"{limit}"
            help_text = f"Select a maximum of {limit} genomes"
            targets_validators.append(
                MaxLengthValidator(
                    limit,
                    message=f"Please select at most {limit} genomes")
            )
        elif limit is not None and limit_type == "lower":
            help_text = f"Select a minimum of {limit} genomes"
            targets_validators.append(
                MinLengthValidator(
                    limit,
                    message=f"Please select at least {limit} genomes")
            )

        targets = forms.MultipleChoiceField(
            choices=accession_choices,
            widget=forms.SelectMultiple(attrs=attrs),
            help_text=help_text,
            required=True,
            validators=targets_validators)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-1 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-4 col-md-6 col-sm-6'
            self.helper.form_action = action
            self.helper.layout = Layout(

                Fieldset("",
                         Column(
                             Row('targets'),
                             Submit('submit', 'Compare %s' %
                                    label, style="margin-top:15px"),
                             css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                         )
            )

        def get_taxids(self):
            indices = self.cleaned_data["targets"]
            taxids = []
            for index in indices:
                taxid = rev_index[int(index)]
                taxids.append(taxid)
            return taxids

    return VennForm


def make_circos_form(database_name):

    accession_choices, reverse_index = get_accessions(database_name)

    class CircosForm(forms.Form):
        circos_reference = forms.ChoiceField(
            choices=accession_choices,
            widget=forms.Select(attrs={
                "class": "selectpicker",
                "data-live-search": "true"
                }
            ))
        targets = forms.MultipleChoiceField(
            choices=accession_choices,
            widget=forms.SelectMultiple(attrs={'size': '1',
                                               "class": "selectpicker",
                                               "data-live-search": "true",
                                               "data-max-options": "10"}),
            required=False)

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = "post"
            self.helper.form_action = "circos"

            self.helper.layout = Layout(
                Fieldset(
                    " ",
                    Row('circos_reference'),
                    Row('targets', style="margin-top:1em"),
                    Submit('submit_circos', 'Submit',
                           style="padding-left:15px; margin-top:15px; margin-bottom:15px "),
                    css_class="col-lg-5 col-md-6 col-sm-6")
            )

        def save(self):
            self.reference = self.cleaned_data["reference"]
            self.get_region = self.cleaned_data["get_region"]
            self.region = self.cleaned_data["region"]

        def get_target_taxids(self):
            indices = self.cleaned_data["targets"]
            taxids = []
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
        checkbox_accessions = forms.BooleanField(
            required=False,
            label="Distinguish plasmids from chromosomes")

        checkbox_single_copy = forms.BooleanField(
            required=False,
            label="Only consider single copy %s" % label)

        orthologs_in = forms.MultipleChoiceField(
            label=f"{label} conserved in",
            choices=accession_choices,
            widget=forms.SelectMultiple(attrs={'size': '%s' % "17",
                                               "class": "selectpicker",
                                               "data-live-search": "true"}),
            required=True)

        no_orthologs_in = forms.MultipleChoiceField(
            label="%s absent from (optional)" % label,
            choices=accession_choices,
            widget=forms.SelectMultiple(attrs={'size': '%s' % "17",
                                               "class": "selectpicker remove-example",
                                               "data-live-search": "true"}),
            required=False)

        new_choices = [['None', 'None']] + accession_choices

        frequency_choices = ((i, i) for i in range(len(accession_choices)))
        frequency = forms.ChoiceField(
            choices=frequency_choices,
            label='Missing data (optional)',
            required=False)

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
                elif is_plasmid:
                    continue
                else:
                    taxids.append(taxid)
            return taxids, plasmids

        def get_n_missing(self):
            return int(self.cleaned_data["frequency"])

        def get_include_choices(self):
            return self.extract_choices((int(i) for i in self.cleaned_data["orthologs_in"]))

        def get_exclude_choices(self):
            return self.extract_choices((int(i) for i in self.cleaned_data["no_orthologs_in"]))

        def __init__(self, *args, **kwargs):
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.form_action = action
            self.helper.layout = Layout(Fieldset(
                " ",
                Column(
                    Row('checkbox_accessions', style="padding-left:15px"),
                    Row('checkbox_single_copy', style="padding-left:15px"),
                    Row(
                        Column(
                            "orthologs_in",
                            css_class='form-group col-lg-6 col-md-6 col-sm-12'),
                        Column(
                            "no_orthologs_in",
                            css_class='form-group col-lg-6 col-md-6 col-sm-12')
                        ),
                    Column(
                        Row('frequency'),
                        Submit('submit',
                               'Compare %s' % label,
                               style="margin-top:15px"),
                        css_class='form-group col-lg-12 col-md-12 col-sm-12'
                        ),
                    css_class="col-lg-8 col-md-8 col-sm-12"
                    )
                )
            )

            super(ExtractForm, self).__init__(*args, **kwargs)

        def clean(self):
            cleaned_data = super(ExtractForm, self).clean()
            self.included_taxids, self.included_plasmids = self.get_include_choices()
            self.excluded_taxids, self.excluded_plasmids = self.get_exclude_choices()
            self.n_missing = self.get_n_missing()
            self.n_included = len(self.included_taxids)
            if self.included_plasmids is not None:
                self.n_included += len(self.included_plasmids)
            if self.n_missing >= self.n_included:
                err = ValidationError(
                    "This must be smaller than the number of included genomes.",
                    code="invalid")
                self.add_error("frequency", err)
            return cleaned_data

    return ExtractForm


def make_module_overview_form(db, sub_sub_cat=False):
    attrs = {'size': '1', "class": "selectpicker",
             "data-live-search": "true",
             "data-actions-box": "true"}
    if sub_sub_cat:
        categories = db.get_module_sub_categories()

        CHOICES = [(cat_id, cat) for cat_id, cat in categories]

        class ModuleCatChoice(forms.Form):
            subcategory = forms.ChoiceField(
                choices=CHOICES, widget=forms.Select(attrs=attrs))

    else:
        categories = db.get_module_categories()

        CHOICES = [(cat_id, cat) for cat_id, cat in categories]

        class ModuleCatChoice(forms.Form):
            category = forms.ChoiceField(
                choices=CHOICES, widget=forms.Select(attrs=attrs))

    return ModuleCatChoice


def make_pathway_overview_form(db):

    pathway_id = db.get_pathways()
    choices = [(path_id, path_desc) for path_id, path_desc in pathway_id]

    class ModuleCatChoice(forms.Form):
        pathway = forms.ChoiceField(choices=choices, widget=forms.Select(
            attrs={"class": "selectpicker", "data-live-search": "true"}))

    return ModuleCatChoice


def make_single_genome_form(db):
    accession_choices, rev_index = get_accessions(db)

    class SingleGenomeForm(forms.Form):
        genome = forms.ChoiceField(choices=accession_choices, widget=forms.Select(
            attrs={"class": "selectpicker", "data-live-search": "true"}))

        def get_genome(self):
            target = self.cleaned_data["genome"]
            return rev_index[int(target)]

    return SingleGenomeForm


def make_blast_form(biodb):
    accession_choices, rev_index = get_accessions(biodb, all=True)

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
        target = forms.ChoiceField(
            choices=accession_choices,
            widget=forms.Select(attrs={"class": "selectpicker",
                                       "data-live-search": "true"})
            )

        input_help = "This can be either an amino-acid or a nucleotide "\
                     "sequence, or a set (one or more) of fasta sequences."
        blast_input = forms.CharField(
            widget=forms.Textarea(attrs={"placeholder": input_help, "rows": 10}))

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.layout = Layout(
                Fieldset(
                    "",
                    Row(
                        Column("blast", css_class='col-lg-4 col-md-4 col-sm-12'),
                        Column("max_number_of_hits", css_class='col-lg-4 col-md-4 col-sm-12'),
                        Column('target', css_class='col-lg-4 col-md-4 col-sm-12'),
                    ),
                    Row(Column('blast_input', css_class='col-lg-12 col-md-12 col-sm-12')),
                    Submit('submit', 'Submit',
                           style="padding-left:15px; margin-top:1em; margin-bottom:15px "),
                    css_class="col-lg-10 col-md-10 col-sm-12")
            )
            super(BlastForm, self).__init__(*args, **kwargs)

        def _get_records(self):
            input_sequence = self.cleaned_data['blast_input']

            if '>' in input_sequence:
                self.no_query_name = False
                try:
                    records = [i for i in SeqIO.parse(
                        StringIO(input_sequence), 'fasta')]
                    for record in records:
                        if len(record.seq) == 0:
                            raise ValidationError(
                                "Empty sequence in input", code="invalid")

                except Exception:
                    raise ValidationError(
                        "Error while parsing the fasta query", code="invalid")
            else:
                self.no_query_name = True
                input_sequence = "".join(input_sequence.split()).upper()
                records = [SeqRecord(Seq(input_sequence))]
            return records

        def _check_sequence_contents(self, records):
            dna = set("ATGCNRYKMSWBDHV")
            prot = set('ACDEFGHIKLMNPQRSTVWYXZJOU')
            sequence_set = set()
            for rec in records:
                sequence_set = sequence_set.union(set(rec.seq.upper()))
            check_seq_DNA = sequence_set - dna
            check_seq_prot = sequence_set - prot

            blast_type = self.cleaned_data["blast"]
            if check_seq_prot and blast_type in ["blastp", "tblastn"]:
                plural = len(check_seq_prot) > 1
                wrong_chars = ", ".join(check_seq_prot)
                errmsg = (f"Unexpected character{'s' if plural else ''}"
                          f" in amino-acid query: {wrong_chars}")
                raise ValidationError(errmsg, code="invalid")

            elif check_seq_DNA and blast_type in ["blastn", "blastn_ffn",
                                                  "blast_fna", "blastx"]:
                wrong_chars = ", ".join(check_seq_DNA)
                plural = len(check_seq_DNA) > 1
                errmsg = (f"Unexpected character{'s' if plural else ''}"
                          f" in nucleotide query: {wrong_chars}")
                raise ValidationError(errmsg, code="invalid")

        def clean_blast_input(self):
            self.records = self._get_records()
            self._check_sequence_contents(self.records)
            return self.cleaned_data['blast_input']

        def get_target(self):
            target = self.cleaned_data["target"]
            if target == "all":
                return target
            return rev_index[int(target)]

    return BlastForm


def get_groups(db):
    return [(el[0], el[0]) for el in db.get_groups()]


def make_gwas_form(biodb):
    import pandas as pd

    group_choices = get_groups(biodb)

    class GwasForm(forms.Form):

        file_help = "CSV file containing 2 columns: taxon IDs or taxon names "\
                    "in the first column and presence (1) or absence (0) "\
                    "of the trait in the second."
        phenotype_file = forms.FileField(help_text=file_help, required=False)

        groups_help = "Select groups having a given phenotype. Genomes not in"\
                      " any of the selected groups will be considered as not"\
                      " having the phenotype."
        groups = forms.MultipleChoiceField(
            choices=group_choices,
            widget=forms.SelectMultiple(attrs={
                'size': '1',
                "class": "selectpicker",
                "data-live-search": "true",
                "multiple data-actions-box": "true"}),
            required=False,
            help_text=groups_help)

        max_number_of_hits = forms.TypedChoiceField(
            choices=[("all", "all"),
                     ("10", "10"),
                     ("20", "20"),
                     ("50", "50"),
                     ("100", "100")],
            coerce=lambda x: int(x) if x != "all" else None)

        bonferroni_cutoff = forms.FloatField(
            initial=0.1,
            widget=forms.NumberInput(
                attrs={"max": 1, "min": 0, "step": 0.01, "default": 0.1})
            )

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.helper = FormHelper()
            self.helper.form_method = 'post'
            self.helper.label_class = 'col-lg-1 col-md-6 col-sm-6'
            self.helper.field_class = 'col-lg-4 col-md-6 col-sm-6'
            self.helper.layout = Layout(
                Fieldset(
                    "",
                    Row("phenotype_file"),
                    Row("groups", style="margin-top:1em"),
                    Row('bonferroni_cutoff', style="margin-top:1em"),
                    Row('max_number_of_hits', style="margin-top:1em"),
                    Submit('submit', 'Submit',
                           style="padding-left:15px; margin-top:1em; "
                                 "margin-bottom:15px "),
                    css_class="col-lg-5 col-md-6 col-sm-6")
            )
            self.db = biodb
            super(GwasForm, self).__init__(*args, **kwargs)

        def clean(self):
            cleaned_data = super(GwasForm, self).clean()
            self.phenotype = self.get_phenotype()
            return cleaned_data

        def get_phenotype(self):
            has_groups = bool(self.cleaned_data["groups"])
            has_file = bool(self.cleaned_data["phenotype_file"])
            if (has_groups and has_file) or not (has_groups or has_file):
                msg = 'You have to provide either "Groups" or '\
                      '"Phenotype file" but not both.'
                raise ValidationError(msg, code="invalid")

            genomes = self.db.get_genomes_description()
            if self.cleaned_data["phenotype_file"]:
                phenotype = pd.read_csv(self.cleaned_data["phenotype_file"],
                                        header=None, names=["taxids", "trait"])
                phenotype["trait"] = phenotype["trait"].astype(bool)
                if all(phenotype.taxids.isin(genomes.index)):
                    return phenotype
                elif all(phenotype.taxids.isin(genomes.description)):
                    mapping = {genome.description: taxid
                               for taxid, genome in genomes.iterrows()}
                    phenotype.taxids = phenotype.taxids.apply(lambda x: mapping[x])
                else:
                    err = ValidationError(
                        "File could not be parsed and matched to genomes.",
                        code="invalid")
                    self.add_error("phenotype_file", err)
                return phenotype
            else:
                taxids_with_phenotype = set(self.db.get_taxids_for_groups(
                    self.cleaned_data["groups"]))
                if not set(genomes.index).difference(taxids_with_phenotype):
                    err = ValidationError(
                        "Your selection is invalid as it contains all genomes.",
                        code="invalid")
                    self.add_error("groups", err)
                phenotype = [[taxid, 1 if taxid in taxids_with_phenotype else 0]
                             for taxid in genomes.index]
                return pd.DataFrame(phenotype, columns=["taxids", "trait"])

    return GwasForm


class CustomPlotsForm(forms.Form):

    help_text = """Entry IDs can be COG, KO, Pfam, VF, AMR or orthogroups.
    IDs should be coma separated.
    You can add custom labels by specifying them together with the entry ID
    separated by a colon (i.e. entryID:label)."""

    entries = forms.CharField(
        widget=forms.Textarea(attrs={'cols': 50, 'rows': 5}),
        required=True, label="Entry IDs", help_text=help_text)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.helper = FormHelper()

        self.helper.form_method = 'post'
        self.helper.layout = Layout(Fieldset(
            "",
            Column(
                Row(Column(
                        "entries",
                        css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                    ),
                Row(Submit('submit', 'Make plot'),
                    css_class='form-group col-lg-12 col-md-12 col-sm-12'),
                css_class="col-lg-8 col-md-8 col-sm-12")
            )
        )

    def get_entries(self):
        raw_entries = self.cleaned_data["entries"].split(",")
        entries = []
        entry2label = {}
        for entry in raw_entries:
            entry = entry.strip()
            if ":" in entry:
                entry, label = entry.split(":", 1)
                entry2label[entry] = label
            entries.append(entry)
        return entries, entry2label
