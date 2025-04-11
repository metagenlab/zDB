from django.shortcuts import render
from django.views import View
from views.mixins import GiViewMixin
from views.utils import format_genome
from views.utils import format_genomic_island
from views.utils import genomic_region_df_to_js
from views.utils import locusx_genomic_region


class GenomicIsland(GiViewMixin, View):
    template = "chlamdb/genomic_island.html"
    view_name = "genomic_island"

    def get(self, request, entry_id, *args, **kwargs):
        self.data = self.get_gi_descriptions([entry_id], transformed=False)
        bioentry = int(self.data.iloc[0]["bioentry.bioentry_id"])
        cluster_id = int(self.data.iloc[0]["cluster_id"])

        cluster_descr = self.get_hit_descriptions([cluster_id]).iloc[0]

        self.data = self.transform_data(self.data).iloc[0]
        all_infos, wd_start, wd_end, contig_size, contig_topology = (
            locusx_genomic_region(
                self.db,
                bioentry=bioentry,
                window_start=self.data.start_pos,
                window_stop=self.data.end_pos,
            )
        )
        genomic_region = genomic_region_df_to_js(
            all_infos, wd_start, wd_end, contig_size, contig_topology
        )

        window_size = wd_end - wd_start
        context = self.get_context(
            organism=self.data.taxon_id,
            gis_id=self.data.gis_id,
            cluster_id=self.data.cluster_id,
            cluster_average_size=cluster_descr.length,
            bioentry=bioentry,
            start_pos=self.data.start_pos,
            end_pos=self.data.end_pos,
            island_size=self.data.end_pos - self.data.start_pos,
            description=self.data.gis_id,
            genomic_region=genomic_region,
            window_size=window_size,
        )
        return render(request, self.template, context)
