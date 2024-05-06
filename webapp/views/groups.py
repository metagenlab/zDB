from collections import namedtuple

from chlamdb.forms import make_group_add_form
from django.shortcuts import redirect, render
from django.urls import reverse
from django.views import View

from views.mixins import BaseViewMixin, GenomesTableMixin
from views.object_type_metadata import GroupMetadata


class GroupsOverview(BaseViewMixin, View, GroupMetadata):

    template = 'chlamdb/groups_overview.html'
    view_name = "groups"

    @property
    def metadata(self):
        return namedtuple("metadata", "description")(self.overview_description)

    def get(self, request, *args, **kwargs):
        groups = [self.format_entry(group[0], to_url=True)
                  for group in self.db.get_groups()]
        return render(request, self.template, self.get_context(groups=groups))


class GroupDetails(BaseViewMixin, View, GroupMetadata, GenomesTableMixin):

    template = 'chlamdb/group_details.html'
    view_name = "groups"
    genome_source_object = "group"

    @property
    def metadata(self):
        return namedtuple("metadata", "description")(self.group)

    def get(self, request, group, *args, **kwargs):
        self.group = group
        taxids = self.db.get_taxids_for_groups([group])
        genome_table = self.get_genomes_table(tuple(taxids))
        return render(request, self.template, self.get_context(group=self.group, genome_table=genome_table))


class GroupDelete(BaseViewMixin, View, GroupMetadata, GenomesTableMixin):

    view_name = "groups_delete"
    genome_source_object = "group"

    def post(self, request, group, *args, **kwargs):
        self.group = group
        self.db.delete_group(group)
        return redirect(reverse("groups"))


class GroupAdd(BaseViewMixin, View, GroupMetadata):

    template = 'chlamdb/group_add.html'
    view_name = "groups"
    genome_source_object = "group"
    description = "Add form"

    @property
    def metadata(self):
        return namedtuple("metadata", "description")(self.description)

    def dispatch(self, request, *args, **kwargs):
        self.form_class = make_group_add_form(self.db)
        return super(GroupAdd, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, self.template, self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            return render(request, self.template, self.get_context())

        group_name = self.form.cleaned_data["group_name"]
        genomes = self.form.get_taxids()

        self.db.load_data_into_table("groups", [(group_name,)])
        self.db.load_data_into_table(
            "taxon_in_group", [(group_name, taxid) for taxid in genomes])
        self.db.commit()
        return redirect(reverse("groups", args=[group_name]))
