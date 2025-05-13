from django.shortcuts import render
from django.views import View
from views.mixins import GiViewMixin
from views.utils import genomic_region_df_to_js
from views.utils import locusx_genomic_region
from views.utils import optional2status


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
        seqids = all_infos.index.unique().tolist()
        to_highlight = {}
        if optional2status.get("vf", False):
            vfs = self.db.vf.get_hits_from_seqids(seqids, columns=("seqid",))
            to_highlight.update(
                {
                    el: "purple"
                    for el in self.db.get_proteins_info(
                        vfs.seqid.to_list(), to_return=["locus_tag"], as_df=True
                    ).get("locus_tag", [])
                }
            )
        if optional2status.get("amr", False):
            amrs = self.db.get_amr_hits_from_seqids(seqids, columns=("seqid",))
            to_highlight.update(
                {
                    el: "magenta"
                    for el in self.db.get_proteins_info(
                        amrs.seqid.to_list(), to_return=["locus_tag"], as_df=True
                    ).get("locus_tag", [])
                }
            )
        window_size = wd_end - wd_start
        context = self.get_context(
            organism=self.data.taxon_id,
            gis_id=self.data.gis_id,
            cluster_id=self.data.cluster_id,
            cluster_average_size=cluster_descr.length,
            bioentry=self.data.bioentry,
            start_pos=self.data.start_pos,
            end_pos=self.data.end_pos,
            island_size=self.data.end_pos - self.data.start_pos,
            description=self.data.gis_id,
            genomic_region=genomic_region,
            to_highlight=to_highlight,
            window_size=window_size,
        )
        return render(request, self.template, context)
