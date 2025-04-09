class BaseQueries:
    where_clause_mapping = {}

    def __init__(self, db):
        self.db = db

    def __getattr__(self, key):
        return getattr(self.db, key)

    def get_number_of_entries(self):
        query = f"SELECT COUNT(*) FROM {self.description_table}"
        return self.server.adaptor.execute_and_fetchall(query)[0][0]

    def get_hit_descriptions(self, hit_ids, columns=None, as_dataframe=True):
        if columns is None:
            columns = self.default_columns
        if hit_ids is None:
            where = ""
        else:
            entries = self.gen_placeholder_string(hit_ids)
            where = f"WHERE descr.{self.id_col} IN ({entries})"

        col_selection = ", ".join([f"descr.{col}" for col in columns])
        query = (
            f"SELECT {col_selection} FROM {self.description_table} as descr {where};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, hit_ids)
        if as_dataframe:
            return self.to_pandas_frame(results, columns)
        return results

    def gen_where_clause(self, search_on, entries):
        entries = self.gen_placeholder_string(entries)

        if search_on == "bioentry":
            where_clause = f" bioentry.bioentry_id IN ({entries}) "
        elif search_on == "seqid":
            where_clause = f" hsh.seqid IN ({entries}) "
        elif search_on == self.id_col:
            where_clause = f" {self.hit_table}.{self.id_col} IN ({entries}) "
        elif search_on in self.where_clause_mapping:
            where_clause = f" {self.where_clause_mapping[search_on]} IN ({entries}) "
        elif search_on == "taxid":
            where_clause = f" bioentry.taxon_id IN ({entries}) "
        else:
            raise RuntimeError(f"Searching on {search_on} is not supported")
        return where_clause

    def get_hits(self, taxids):
        where_clause = self.gen_where_clause("taxid", taxids)
        columns = [
            "bioentry.taxon_id",
            f"{self.hit_table}.{self.id_col}",
            "evalue",
            "prot_name",
            "vfid",
            "category",
            "vf_category_id",
        ]
        query = (
            f"SELECT {', '.join(columns)} "
            f"FROM {self.hit_table} "
            f"INNER JOIN {self.description_table} ON {self.description_table}.{self.id_col} = {self.hit_table}.{self.id_col} "
            f"INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = {self.hit_table}.hsh "
            "INNER JOIN seqfeature ON hsh.seqid = seqfeature.seqfeature_id "
            "INNER JOIN bioentry ON seqfeature.bioentry_id = bioentry.bioentry_id "
            f"WHERE {where_clause} "
        )
        results = self.server.adaptor.execute_and_fetchall(query, taxids)
        df = self.to_pandas_frame(results, [col.split(".")[-1] for col in columns])
        return df

    def get_hit_counts(
        self, ids, indexing="taxid", search_on="taxid", keep_taxid=False, plasmids=None
    ):
        if indexing == "seqid":
            index = "seqfeature.seqfeature_id"
            if keep_taxid:
                index += ", bioentry.taxon_id "
        elif indexing == "taxid":
            index = "bioentry.taxon_id"
        else:
            raise RuntimeError(f"Indexing method not supported: {indexing}")

        where_clause = self.gen_where_clause(search_on, ids)

        plasmid_join = ""
        if plasmids is not None:
            index += ", CAST(is_plasmid.value AS int) "
            subclause = self.gen_where_clause(search_on, plasmids)
            where_clause = (
                f"({where_clause} AND is_plasmid.value=0) "
                f" OR ({subclause} AND is_plasmid.value=1)"
            )
            plasmid_join = (
                "INNER JOIN bioentry_qualifier_value AS is_plasmid ON "
                "  is_plasmid.bioentry_id=bioentry.bioentry_id "
                "INNER JOIN term AS plasmid_term ON plasmid_term.term_id=is_plasmid.term_id "
                '  AND plasmid_term.name="plasmid"'
            )

        query = (
            f"SELECT {index}, {self.id_col}, COUNT(*) "
            f"FROM {self.hit_table} "
            f"{self.join_bioentry}"
            f"{plasmid_join}"
            f"WHERE {where_clause} "
            f"GROUP BY {index}, {self.id_col};"
        )

        all_ids = ids
        if plasmids:
            all_ids += plasmids
        results = self.server.adaptor.execute_and_fetchall(query, all_ids)

        # ugly code, open for improvements
        if indexing == "taxid" or indexing == "bioentry":
            column_names = [indexing, self.id_col, "count"]
            index = [indexing, self.id_col]
            if plasmids is not None:
                column_names.insert(1, "plasmid")
                index.insert(1, "plasmid")
            df = self.to_pandas_frame(results, column_names)
            if df.empty:
                return df
            df = df.set_index(index).unstack(level=0, fill_value=0)

            if plasmids is not None:
                return df.unstack(level=0, fill_value=0)
            else:
                df.columns = [col for col in df["count"].columns.values]
                if search_on == "taxid":
                    for taxid in ids:
                        if taxid not in df.columns:
                            df[taxid] = 0

        elif indexing == "seqid":
            if plasmids is not None:
                raise RuntimeError("Not implemented for now")
            header = ["seqid", self.id_col]
            if keep_taxid:
                header.append("taxid")
                results = (
                    (seqid, obj_id, taxid) for seqid, taxid, obj_id, count in results
                )
            else:
                results = ((seqid, obj_id) for seqid, obj_id, count in results)

            df = self.to_pandas_frame(results, header)
            df = df.set_index(["seqid"])
        return df


class GIQueries(BaseQueries):
    description_table = "genomic_islands"
    id_col = "gis_id"
    default_columns = [id_col, "bioentry_id", "start_pos", "end_pos"]

    def get_containing_genomic_islands(self, bioentry_id, start, stop):
        sql = f"SELECT {self.id_col}, start_pos, end_pos FROM {self.description_table} WHERE bioentry_id=? AND (? BETWEEN start_pos AND end_pos OR ? BETWEEN start_pos AND end_pos)"
        return self.server.adaptor.execute_and_fetchall(sql, [bioentry_id, start, stop])


class VFQueries(BaseQueries):
    hit_table = "vf_hits"
    description_table = "vf_defs"
    id_col = "vf_gene_id"

    where_clause_mapping = {"vf": id_col}
    default_columns = [id_col, "prot_name", "vfid", "category", "vf_category_id"]

    join_bioentry = (
        f"INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = {hit_table}.hsh "
        "INNER JOIN seqfeature ON hsh.seqid = seqfeature.seqfeature_id "
        "INNER JOIN bioentry ON seqfeature.bioentry_id = bioentry.bioentry_id "
    )
