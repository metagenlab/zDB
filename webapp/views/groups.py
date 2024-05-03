from collections import namedtuple

from django.shortcuts import render
from django.views import View

from views.mixins import BaseViewMixin
from views.object_type_metadata import GroupMetadata


class GroupsOverview(BaseViewMixin, View, GroupMetadata):

    template = 'chlamdb/groups_overview.html'
    view_name = "groups"

    @property
    def metadata(self):
        return namedtuple("metadata", "description")(self.overview_description)

    def get(self, request, *args, **kwargs):
        groups = [group[0] for group in self.db.get_groups()]
        return render(request, self.template, self.get_context(groups=groups))
