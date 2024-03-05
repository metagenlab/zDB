

import pandas as pd
from django.shortcuts import render
from django.views import View

from views.mixins import (AmrViewMixin, BaseViewMixin, CogViewMixin,
                          KoViewMixin, PfamViewMixin, VfViewMixin)
from views.utils import format_ko, format_ko_modules, format_ko_path


class EntryListViewBase(View, BaseViewMixin):

    _specific_colname_to_header_mapping = {"freq": "Frequency (n genomes)"}

    def get(self, request):
        # retrieve taxid list
        genomes_data = self.db.get_genomes_infos()
        self.taxids = [str(i) for i in genomes_data.index.to_list()]

        table_data = self.get_table_data()

        context = self.get_context(
            table_headers=self.table_headers,
            table_data_accessors=self.table_data_accessors,
            table_data=table_data)
        return render(request, 'chlamdb/entry_list.html', context)

    def get_table_data(self):
        # retrieve hits
        all_hits = self.get_hit_counts(self.taxids,
                                       search_on="taxid",
                                       indexing="taxid")
        # retrieve descriptions
        descriptions = self.get_hit_descriptions(all_hits.index.to_list())

        # count frequency and n genomes
        pfam_count = all_hits.sum(axis=1)
        pfam_freq = all_hits[all_hits > 0].count(axis=1)

        # combine into df
        combined_df = descriptions.merge(pfam_count.rename('count'),
                                         left_index=True,
                                         right_index=True)\
                                  .merge(pfam_freq.rename('freq'),
                                         left_index=True,
                                         right_index=True)\
                                  .sort_values(["count", "freq"],
                                               ascending=False)

        combined_df = combined_df.where(combined_df.notna(), "-")
        return combined_df

    @property
    def table_data_accessors(self):
        return super(EntryListViewBase, self).table_data_accessors + ["count", "freq"]


class PfamEntryListView(EntryListViewBase, PfamViewMixin):

    pass


class KoEntryListView(EntryListViewBase, KoViewMixin):

    def get_table_data(self):
        # retrieve entry list
        ko_all = self.db.get_ko_hits(self.taxids,
                                     search_on="taxid",
                                     indexing="taxid")
        # retrieve annotations
        ko_desc = self.db.get_ko_desc(ko_all.index.to_list())
        ko_mod = self.db.get_ko_modules(ko_all.index.to_list())
        ko_path = self.db.get_ko_pathways(ko_all.index.to_list())

        # count frequency and n genomes
        combined_df = pd.DataFrame(ko_all.sum(axis=1).rename('count'))
        ko_freq = ko_all[ko_all > 0].count(axis=1).to_dict()

        combined_df["ko"] = [format_ko(ko, as_url=True) for ko in combined_df.index]
        combined_df["modules"] = [format_ko_modules(ko_mod, ko) if ko in ko_mod else '-' for ko in combined_df.index]
        combined_df["description"] = [ko_desc[ko] for ko in combined_df.index]
        combined_df["pathways"] = [format_ko_path(ko_path, ko) if ko in ko_path else '-' for ko in combined_df.index]
        combined_df["freq"] = [ko_freq[ko] for ko in combined_df.index]

        return combined_df.sort_values(["count", "freq"], ascending=False)


class CogEntryListView(EntryListViewBase, CogViewMixin):

    pass


class AmrEntryListView(EntryListViewBase, AmrViewMixin):

    pass


class VfEntryListView(EntryListViewBase, VfViewMixin):

    pass
