

import pandas as pd
from django.conf import settings
from django.shortcuts import render
from django.views import View
from lib.db_utils import DB

from views.mixins import AmrAnnotationsMixin
from views.utils import (format_amr, format_cog, format_hmm_url, format_ko,
                         format_ko_modules, format_ko_path, format_pfam,
                         my_locals, page2title)


class EntryListViewBase(View):

    entry_type = None
    table_headers = None
    table_data_accessors = None

    def get(self, request):
        page_title = page2title[f"entry_list_{self.entry_type}"]
        self.db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
        # retrieve taxid list
        genomes_data = self.db.get_genomes_infos()
        self.taxids = [str(i) for i in genomes_data.index.to_list()]

        table_data = self.get_table_data()

        context = my_locals({
            "page_title": page_title,
            "entry_type": self.entry_type,
            "table_headers": self.table_headers,
            "table_data_accessors": self.table_data_accessors,
            "table_data": table_data,
            })
        return render(request, 'chlamdb/entry_list.html', context)

    def get_table_data(self):
        raise NotImplementedError()


class PfamEntryListView(EntryListViewBase):

    entry_type = "pfam"
    table_headers = ["Accession", "Description", "Count", "Frequency (n genomes)"]
    table_data_accessors = ["accession", "def", "count", "freq"]

    def get_table_data(self):
        # retrieve entry list
        pfam_all = self.db.get_pfam_hits(self.taxids,
                                         search_on="taxid",
                                         indexing="taxid")
        # retrieve annotations
        pfam_annot = self.db.get_pfam_def(pfam_all.index.to_list())

        # count frequency and n genomes
        pfam_count = pfam_all.sum(axis=1)
        pfam_freq = pfam_all[pfam_all > 0].count(axis=1)
        pfam_annot["accession"] = [format_pfam(pfam, to_url=True)
                                   for pfam in pfam_annot.index]

        # combine into df
        combined_df = pfam_annot.merge(pfam_count.rename('count'),
                                       left_index=True,
                                       right_index=True)\
                                .merge(pfam_freq.rename('freq'),
                                       left_index=True,
                                       right_index=True)\
                                .sort_values(["count"],
                                             ascending=False)

        return combined_df


class KoEntryListView(EntryListViewBase):

    entry_type = "ko"
    table_headers = ["Accession", "Description", "Modules", "Pathways",
                     "Count", "Frequency (n genomes)"]
    table_data_accessors = ["accession", "description", "modules", "pathways",
                            "count", "freq"]

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

        combined_df["accession"] = [format_ko(ko, as_url=True) for ko in combined_df.index]
        combined_df["modules"] = [format_ko_modules(ko_mod, ko) if ko in ko_mod else '-' for ko in combined_df.index]
        combined_df["description"] = [ko_desc[ko] for ko in combined_df.index]
        combined_df["pathways"] = [format_ko_path(ko_path, ko) if ko in ko_path else '-' for ko in combined_df.index]
        combined_df["freq"] = [ko_freq[ko] for ko in combined_df.index]

        return combined_df.sort_values(["count", "freq"], ascending=False)


class CogEntryListView(EntryListViewBase):

    entry_type = "cog"
    table_headers = ["Accession", "Function", "Description", "Count",
                     "Frequency (n genomes)"]
    table_data_accessors = ["accession", "function", "description", "count",
                            "freq"]

    def get_table_data(self):
        # retrieve entry list
        cog_all = self.db.get_cog_hits(self.taxids,
                                       search_on="taxid")
        # retrieve annotations
        cogs_summaries = self.db.get_cog_summaries(
            cog_all.index.tolist(), as_df=True, only_cog_desc=True)

        # count frequency and n genomes
        cog_count = cog_all.sum(axis=1)
        cog_freq = cog_all[cog_all > 0].count(axis=1)
        cogs_summaries["accession"] = [format_cog(cog, as_url=True)
                                       for cog in cogs_summaries.index]

        # combine into df
        combined_df = cogs_summaries.merge(
            cog_count.rename('count'),
            left_index=True,
            right_index=True).merge(
                cog_freq.rename('freq'),
                left_index=True,
                right_index=True).sort_values(["count"], ascending=False)
        return combined_df


class AmrEntryListView(EntryListViewBase, AmrAnnotationsMixin):

    entry_type = "amr"
    table_headers = ["Gene", "Description", "Scope", "Type", "Class",
                     "Subclass", "HMM", "Count", "Frequency (n genomes)"]
    table_data_accessors = ["accession", "seq_name", "scope", "type", "class",
                            "subclass", "hmm", "count", "freq"]

    def get_table_data(self):
        # retrieve entry list
        amr_all = self.db.get_amr_hit_counts(self.taxids,
                                             search_on="taxid",
                                             indexing="taxid")
        # retrieve annotations
        amr_annotations = self.db.get_amr_descriptions(amr_all.index.tolist())

        self.aggregate_amr_annotations(amr_annotations)

        # count frequency and n genomes
        combined_df = pd.DataFrame(amr_all.sum(axis=1).rename('count'))
        freq = amr_all[amr_all > 0].count(axis=1).to_dict()

        # prepare accession
        amr_annotations["accession"] = amr_annotations["gene"].apply(
            format_amr, to_url=True)

        # link hmms
        amr_annotations["hmm"] = amr_annotations["hmm_id"].apply(
            format_hmm_url)

        combined_df["freq"] = [freq[gene] for gene in combined_df.index]
        for col in self.table_data_accessors[:-2]:
            combined_df[col] = [
                amr_annotations[amr_annotations.gene == gene].iloc[0][col]
                for gene in combined_df.index]

        combined_df = combined_df.where(combined_df.notna(), "-")
        return combined_df.sort_values(["count", "freq"], ascending=False)
