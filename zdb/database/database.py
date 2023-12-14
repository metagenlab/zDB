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
            "closest_seq_name tinytext);"
        )
        self.server.adaptor.execute(sql,)

        sql = (
            "CREATE INDEX amr_hits_idx ON amr_hits(hsh);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("amr_hits", data)
