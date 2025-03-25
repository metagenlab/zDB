from django.shortcuts import render
from django.views import View
from views.analysis_view_metadata import EntryListMetadata
from views.mixins import AmrViewMixin
from views.mixins import BaseViewMixin
from views.mixins import CogViewMixin
from views.mixins import GiViewMixin
from views.mixins import KoViewMixin
from views.mixins import PfamViewMixin
from views.mixins import VfViewMixin


class EntryListViewBase(View, BaseViewMixin):
    _specific_colname_to_header_mapping = {"freq": "Frequency (n genomes)"}
    _metadata_cls = EntryListMetadata

    def get(self, request):
        # retrieve taxid list
        genomes_data = self.db.get_genomes_infos()
        self.taxids = [str(i) for i in genomes_data.index.to_list()]

        table_data = self.get_table_data()

        context = self.get_context(
            table_headers=self.table_headers,
            table_data_accessors=self.table_data_accessors,
            table_data=table_data,
        )
        return render(request, "chlamdb/entry_list.html", context)

    def get_table_data(self):
        # retrieve hits
        all_hits = self.get_hit_counts(self.taxids, search_on="taxid", indexing="taxid")
        # retrieve descriptions
        descriptions = self.get_hit_descriptions(all_hits.index.to_list())

        # count frequency and n genomes
        pfam_count = all_hits.sum(axis=1)
        pfam_freq = all_hits[all_hits > 0].count(axis=1)

        # combine into df
        combined_df = (
            descriptions.merge(
                pfam_count.rename("count"), left_index=True, right_index=True
            )
            .merge(pfam_freq.rename("freq"), left_index=True, right_index=True)
            .sort_values(["count", "freq"], ascending=False)
        )

        combined_df = combined_df.where(combined_df.notna(), "-")
        return combined_df

    @property
    def table_data_accessors(self):
        return super(EntryListViewBase, self).table_data_accessors + ["count", "freq"]


class PfamEntryListView(EntryListViewBase, PfamViewMixin):
    pass


class KoEntryListView(EntryListViewBase, KoViewMixin):
    pass


class GiEntryListView(EntryListViewBase, GiViewMixin):
    def get_table_data(self):
        return self.get_hit_descriptions(None)

    @property
    def table_data_accessors(self):
        return super(EntryListViewBase, self).table_data_accessors


class CogEntryListView(EntryListViewBase, CogViewMixin):
    pass


class AmrEntryListView(EntryListViewBase, AmrViewMixin):
    pass


class VfEntryListView(EntryListViewBase, VfViewMixin):
    pass
