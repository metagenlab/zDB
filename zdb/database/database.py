from metagenlab_libs import db_utils


class DB(db_utils.DB):

    @classmethod
    def load_db(cls, db_file, params):
        db = super(DB, cls).load_db(db_file, params)
        return cls(db.server, db.db_name)
