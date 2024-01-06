
import os

from collections import namedtuple

from whoosh.fields import SchemaClass, TEXT, KEYWORD, ID
from whoosh.qparser import MultifieldParser
from whoosh import index


class SearchBarSchema(SchemaClass):
    entry_type = ID(stored=True)
    name = KEYWORD(stored=True)
    description = TEXT(stored=True)
    organism = TEXT(stored=True)
    locus_tag = KEYWORD(stored=True)
    og = KEYWORD(stored=True)


field_list = ["entry_type", "name", "description", "organism", "locus_tag"]
SearchResult = namedtuple("SearchResult", field_list)


class EntryTypes:
    COG  = "C"
    Gene = "G"
    KO   = "K"
    PFAM = "D"
    Module = "M"
    Pathway = "P"


class ChlamdbIndex:

    def new_index(name):
        chlamdb_index = ChlamdbIndex()
        chlamdb_index.index = index.create_in(name, SearchBarSchema)
        chlamdb_index.writer = chlamdb_index.index.writer()
        return chlamdb_index

    def use_index(name):
        chlamdb_index = ChlamdbIndex()
        chlamdb_index.index = index.open_dir(name)
        return chlamdb_index

    def add(self, **kwargs):
        self.writer.add_document(**kwargs)

    def search(self, user_query, limit=10):
        parser = MultifieldParser(field_list, self.index.schema)
        query = parser.parse(user_query)

        for result in self.index.searcher().search(query, limit=limit):
            hsh_res = {}
            for item in field_list:
                hsh_res[item] = result.get(item, None)
            yield SearchResult(**hsh_res)

    def done_adding(self):
        self.writer.commit(optimize=True)
