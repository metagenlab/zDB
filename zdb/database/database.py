from metagenlab_libs import db_utils


class DB(db_utils.DB):

    @classmethod
    def load_db(cls, db_file, params):
        db = super(DB, cls).load_db(db_file, params)
        return cls(db.server, db.db_name)

    def load_amr_hits(self, data):
        sql = (
            "CREATE TABLE amr_hits (hsh INTEGER, gene varchar(20), seq_name tinytext, "
            "scope char(4), type varchar(10), subtype varchar(10), class tinytext, "
            "subclass tinytext, coverage FLOAT, identity FLOAT, closest_seq tinytext, "
            "closest_seq_name tinytext, hmm_id varchar(20));"
        )
        self.server.adaptor.execute(sql,)

        sql = (
            "CREATE INDEX amr_hits_idx ON amr_hits(hsh);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("amr_hits", data)

    def get_amr_hits(self, ids):
        """
        For now we limit that search to AMR type
        """

        columns = ("gene", "scope", "type", "class", "subclass", "coverage",
                   "identity", "closest_seq", "hmm_id")

        query = (
            f"SELECT {', '.join(f'amr.{col}' for col in columns)} "
            "FROM amr_hits AS amr "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh = amr.hsh "
            f"WHERE hsh.seqid IN ({', '.join(str(el) for el in ids)});"
        )

        results = self.server.adaptor.execute_and_fetchall(query)

        df = DB.to_pandas_frame(results, columns)
        return df
