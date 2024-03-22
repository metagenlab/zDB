
from collections import namedtuple

import pandas as pd
from whoosh import index
from whoosh.fields import ID, KEYWORD, TEXT, SchemaClass
from whoosh.qparser import MultifieldParser


class SearchBarSchema(SchemaClass):
    entry_type = ID(stored=True)
    name = KEYWORD(stored=True)
    description = TEXT(stored=True)
    organism = TEXT(stored=True)
    locus_tag = KEYWORD(stored=True)
    og = KEYWORD(stored=True)


field_list = ["entry_type", "name", "description", "organism", "locus_tag"]
SearchResult = namedtuple("SearchResult", field_list)
class EntryType():

    @property
    def _name_prefix(self):
        return self.entry_type

    @property
    def entry_type(self):
        return self.object_type[0].upper()

    def add_to_index(self, index, entry_id, descr):
        index.writer.add_document(entry_type=self.entry_type,
                                  name=self.get_name(entry_id),
                                  description=descr)


class KoEntry(EntryType):

    object_type = "ko"
    _name_format_spec = '05'

    def get_name(self, entry_id):
        if pd.isna(entry_id):
            return None
        return f"{self._name_prefix}{int(entry_id):{self._name_format_spec}}"


class CogEntry(KoEntry):

    object_type = "cog"
    _name_prefix = "COG"
    _name_format_spec = '04'


class ModuleEntry(EntryType):

    object_type = "module"
    _name_format_spec = '05d'

    def get_name(self, entry_id):
        return f"{self._name_prefix}{entry_id:{self._name_format_spec}}"


class PathwayEntry(ModuleEntry):

    object_type = "pathway"
    _name_prefix = "map"


class PfamEntry(KoEntry):

    entry_type = "D"
    object_type = "pfam"
    _name_prefix = "PF"


class AmrEntry(EntryType):

    entry_type = "R"
    object_type = "amr"

    def get_name(self, entry_id):
        return entry_id


class VfEntry(EntryType):

    object_type = "vf"

    def get_name(self, entry_id):
        return entry_id


class GeneEntry():

    entry_type = "G"
    object_type = "locus"

    def add_to_index(self, index, locus_tag, gene, product, organism, og):
        index.writer.add_document(entry_type=self.entry_type,
                                  locus_tag=locus_tag,
                                  name=gene,
                                  description=product,
                                  organism=organism,
                                  og=og)


entry_classes = [KoEntry, CogEntry, ModuleEntry, PathwayEntry, PfamEntry,
                 AmrEntry, VfEntry, GeneEntry]
entry_type_to_object_type = {cls.entry_type: cls.object_type
                             for cls in entry_classes}

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
