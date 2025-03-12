from django.shortcuts import render
from django.views import View
from views.mixins import BaseViewMixin
from views.utils import format_genome
from views.utils import format_genomic_island
from views.utils import genomic_region_df_to_js
from views.utils import locusx_genomic_region


class GenomicIsland(BaseViewMixin, View):
    template = "chlamdb/genomic_island.html"
    view_name = "genomic_island"

    def get(self, request, entry_id, *args, **kwargs):
        self.gis_id, self.bioentry_id, self.start_pos, self.end_pos = (
            self.db.get_genomic_island(entry_id)
        )
        all_infos, wd_start, wd_end, contig_size, contig_topology = (
            locusx_genomic_region(
                self.db,
                bioentry=self.bioentry_id,
                window_start=self.start_pos,
                window_stop=self.end_pos,
            )
        )
        genomic_region = genomic_region_df_to_js(
            all_infos, wd_start, wd_end, contig_size, contig_topology
        )
        window_size = wd_end - wd_start
        context = self.get_context(
            entry_id=self.gis_id,
            start_pos=self.start_pos,
            end_pos=self.end_pos,
            description=self.description,
            genomic_region=genomic_region,
            window_size=window_size,
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
            f"({format_genomic_island(self.gis_id, accession, self.start_pos, self.end_pos)})"
        )
