from django.shortcuts import render
from django.views import View
from views.mixins import BaseViewMixin
from views.utils import format_genome


class GenomicIsland(BaseViewMixin, View):
    template = "chlamdb/genomic_island.html"
    view_name = "genomic_island"

    def get(self, request, entry_id, *args, **kwargs):
        self.gis_id, self.bioentry_id, self.start_pos, self.end_pos = (
            self.db.get_genomic_island(entry_id)
        )
        context = self.get_context(
            entry_id=self.gis_id,
            start_pos=self.start_pos,
            end_pos=self.end_pos,
            description=self.description,
        )
        return render(request, self.template, context)

    @property
    def description(self):
        taxon_id, accession, description = self.db.server.adaptor.execute_and_fetchall(
            "SELECT taxon_id, accession, description from bioentry WHERE bioentry_id=?",
            [self.bioentry_id],
        )[0]
        return (
            f"Genomic island: "
            f"{format_genome((taxon_id, description))} "
            f"({accession}: {self.start_pos} - {self.end_pos})"
        )
