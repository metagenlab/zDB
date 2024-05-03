from collections import namedtuple

from django.shortcuts import render
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
    view_name = "group"
    genome_source_object = "group"

    @property
    def metadata(self):
        return namedtuple("metadata", "description")(self.group)

    def get(self, request, group, *args, **kwargs):
        self.group = group
        taxids = self.db.get_taxids_for_groups([group])
        genome_table = self.get_genomes_table(tuple(taxids))
        return render(request, self.template, self.get_context(genome_table=genome_table))
