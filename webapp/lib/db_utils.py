import os
import sys

from BioSQL import BioSeqDatabase
from Bio.SeqRecord import SeqRecord #added

from Bio.Seq import Seq
from Bio.SeqUtils import GC

import sqlite3
import pandas as pd


# This file defines a class DB, that encapsulates all the SQL requests
# necessary to create the zDB database.
# In the future, the goal is to import all database queries needed by the
# zDB website as methods of this class.
#
# This improves code readability by removing SQL queries from the main python
# code and more importantly, it would allow to change of database without having
# to modify zDB's code.


# to litteral
# encases the string into quotes
def quote(v):
    return f"\"{v}\""


class NoPhylogenyException(Exception):
    """ Used when the user tries to get a phylogeny for an orthogroup
     too small to have one
    """
    pass


class DB:

    def __init__(self, server, db_name):
        self.server = server
        self.db_name = db_name
        self.conn_ncbi_taxonomy = None
        self.conn_refseq = None
        # this will need to be changed in case a MySQL database is used
        self.placeholder = "?"

    # the next two methods are necessary for DB objects to be used
    # in 'with' blocks.
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.server.close()

    def create_indices_on_cds(self):
        sql_index1 = 'create index ftgcga on feature_tables_genomes_cds(genome_accession)'    
        sql_index2 = 'create index ftgctx on feature_tables_genomes_cds(taxon_id)'
        sql_index3 = 'create index ftgrga on feature_tables_genomes_rrna(taxon_id)'
        sql_index4 = 'create index ftgrtx on feature_tables_genomes_rrna(genome_accession)'    
        sql_index5 = 'create index ftcain on feature_tables_cds_accessions(id_name)'
        sql_index6 = 'create index ftcait on feature_tables_cds_accessions(id_type)'
        
        self.server.adaptor.execute(sql_index1,)
        self.server.adaptor.execute(sql_index2,)
        self.server.adaptor.execute(sql_index3,)
        self.server.adaptor.execute(sql_index4,)
        self.server.adaptor.execute(sql_index5,)
        self.server.adaptor.execute(sql_index6,)


    def get_taxid_from_accession(self, accessions, look_on="locus_tag"):
        plc = self.gen_placeholder_string(accessions)
        if look_on=="locus_tag":
            select = " acc.value , t1.taxon_id "
            query = (
                " INNER JOIN seqfeature AS fet ON fet.bioentry_id = t1.bioentry_id "
                " INNER JOIN seqfeature_qualifier_value AS acc "
                "     ON acc.seqfeature_id = fet.seqfeature_id "
                " INNER JOIN term AS gene_term ON gene_term.term_id = fet.type_term_id "
                "     AND gene_term.name = \"CDS\" "
                " INNER JOIN term AS lc ON lc.term_id = acc.term_id AND "
                "     lc.name = \"locus_tag\" "
                f" WHERE acc.value IN ({plc})"
            )
        elif look_on=="contig":
            select = "t1.name, t1.taxon_id "
            query = f" WHERE t1.name IN ({plc})"
        else:
            raise RuntimeError("Unknown option " + look_on + " expect either locus_tag or contig")

        sql = (
            f"SELECT {select} "
            "FROM bioentry t1 "
            f"{query};"
        )
        results = self.server.adaptor.execute_and_fetchall(sql, accessions)
        header = ["accession", "taxid"]
        return DB.to_pandas_frame(results, header).set_index("accession")


    def get_taxid_from_seqid(self, seqids):
        query = ",".join("?" for _ in seqids)
        sql = (
            "SELECT fet.seqfeature_id, entry.taxon_id "
            "FROM seqfeature AS fet "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=fet.bioentry_id "
            f"WHERE fet.seqfeature_id IN ({query});"
        )
        results = self.server.adaptor.execute_and_fetchall(sql, seqids)
        return {seqid: taxon_id for seqid, taxon_id in results}


    def create_seq_hash_to_seqid(self, to_load):
        sql = (
            f"CREATE TABLE sequence_hash_dictionnary (hsh INTEGER, seqid INTEGER," 
            " PRIMARY KEY(hsh, seqid), "
            " FOREIGN KEY(seqid) REFERENCES seqfeature(seqfeature_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("sequence_hash_dictionnary", to_load)

        sql = "CREATE INDEX shd_hsh on sequence_hash_dictionnary (hsh)"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX shd_seqid on sequence_hash_dictionnary (seqid)"
        self.server.adaptor.execute(sql)

    def get_seq_hash_to_seqid(self, to_load):
        query = (
            "SELECT * FROM sequence_hash_dictionnary;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    
    def get_refseq_taxonomy(self, taxids):
        placeholder = self.gen_placeholder_string(taxids)
        query = (
            "SELECT taxid, value "
            "FROM taxonomy_mapping "
            f"WHERE taxid IN ({placeholder});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, taxids)
        header = ["taxid", "taxonomic_name"]
        return DB.to_pandas_frame(results, header).set_index("taxid")


    def get_refseq_hits(self, seqids):
        placeholder = self.gen_placeholder_string(seqids)

        query = (
            "SELECT ref_homolog.sseqid, length, evalue, bitscore, gapopen, pident "
            "FROM sequence_hash_dictionnary AS hsh "
            "INNER JOIN diamond_refseq AS ref_homolog ON hsh.hsh=ref_homolog.seq_hash "
            f"WHERE seqid IN ({placeholder});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, seqids)
        return DB.to_pandas_frame(results, ["match_id", "length", "evalue", "bitscore", "gaps", "pident"])


    def get_refseq_matches_info(self, match_ids, search_on="match_id"):
        placeholder = self.gen_placeholder_string(match_ids)

        if search_on=="accession":
            column = "accession"
        elif search_on=="match_id":
            column = "match_id"
        else:
            raise RuntimeError("Unsupported search term")

        query = (
            "SELECT match_id, accession, organism, description "
            "FROM diamond_refseq_match_id "
            f"WHERE {column} IN ({placeholder});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, match_ids)
        df = DB.to_pandas_frame(results, ["match_id", "accession", "organism", "description"])
        return df.set_index("match_id")


    # Returns all refseq hits associated with a given orthogroup
    def get_diamond_match_for_og(self, og, sort_by_evalue=False):
        sort = ""
        if sort_by_evalue:
            sort = "ORDER BY eval ASC"
        query = (
            "SELECT match_id.accession, MIN(hit.evalue) AS eval "
            "FROM og_hits "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.seqid = og_hits.seqid "
            "INNER JOIN diamond_refseq AS hit ON hit.seq_hash = hsh.hsh "
            "INNER JOIN diamond_refseq_match_id AS match_id ON hit.sseqid = match_id.match_id "
            f"WHERE og_hits.orthogroup={og} "
            "GROUP BY match_id.accession "
            f"{sort};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return DB.to_pandas_frame(results, ["accession", "evalue"])


    def get_all_orthogroups(self):
        query = ("SELECT orthogroup_id, og_size "
            "FROM gene_phylogeny;"
        )
        return self.server.adaptor.execute_and_fetchall(query)


    def get_n_orthogroups(self, only_core=False):
        where = ""
        if only_core:
            where = "WHERE is_core=1"

        query = (
            "SELECT COUNT(*)"
            "FROM gene_phylogeny "
            f"{where};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return results[0][0]


    def get_all_sequences_for_orthogroup(self, orthogroup):
        from Bio import SeqRecord
        from Bio.Seq import Seq

        query = (
            "SELECT locus.value, seq.value "
            "FROM og_hits AS ortho "
            "INNER JOIN seqfeature_qualifier_value AS seq ON seq.seqfeature_id=ortho.seqid "
            " INNER JOIN term AS seq_term "
            "   ON seq_term.term_id=seq.term_id AND seq_term.name=\"translation\""
            "INNER JOIN seqfeature_qualifier_value locus ON locus.seqfeature_id=ortho.seqid "
            " INNER JOIN term as locus_term "
            "   ON locus_term.term_id=locus.term_id AND locus_term.name=\"locus_tag\" "
            f"WHERE ortho.orthogroup = {orthogroup};"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        tab = []
        for result in results:
            record = SeqRecord.SeqRecord(Seq(result[1]), id=result[0])
            tab.append(record)
        return tab


    def get_n_swissprot_homologs(self, seqid):
        query = (
            "SELECT COUNT(*) "
            "FROM sequence_hash_dictionnary AS hsh "
            "INNER JOIN swissprot_hits AS hits ON hits.hsh=hsh.hsh "
            f"WHERE hsh.seqid =?;"
        )
        return self.server.adaptor.execute_and_fetchall(query, [seqid])[0][0]


    def get_swissprot_homologs(self, seqids, indexing=None):
        """
        If indexing is accession: only returns the accession and definition of the
         swissprot hits.
        """

        plcder = self.gen_placeholder_string(seqids)
        sel = (
            "hsh.seqid, defs.swissprot_id, defs.definition, "
            " defs.taxid, defs.organism, defs.gene, hits.evalue, hits.score, "
            "    hits.perc_id, hits.leng, hits.gaps, defs.pe "
        )
        cols = ["seqid", "accession", "definition", "taxid", "organism", "gene",
                "evalue", "bitscore", "perc_id", "match_len", "gaps", "pe"]
        groupby = ""

        if indexing=="accession":
            sel = "defs.swissprot_id, defs.definition "
            groupby = "GROUP BY defs.swissprot_id "
            cols = ["accession", "definition"]

        query = (
            f"SELECT {sel}"
            "FROM sequence_hash_dictionnary AS hsh "
            "INNER JOIN swissprot_hits AS hits ON hits.hsh=hsh.hsh "
            "INNER JOIN swissprot_defs AS defs ON defs.prot_id=hits.prot_id "
            f"WHERE hsh.seqid IN ({plcder}) {groupby};"
        )
        args = seqids
        results = self.server.adaptor.execute_and_fetchall(query, args)
        return DB.to_pandas_frame(results, cols)


    def create_swissprot_tables(self):
        query = (
            "CREATE TABLE swissprot_hits ("
            " hsh INT, prot_id INT, evalue INT, score INT, perc_id INT, gaps INT, leng INT"
            ");"
        )
        self.server.adaptor.execute(query,)
        query = (
            "CREATE INDEX sphi ON swissprot_hits(hsh);"
        )
        self.server.adaptor.execute(query,)
        query = (
            "CREATE TABLE swissprot_defs ("
            " prot_id INT, swissprot_id TEXT, definition TEXT, taxid INT, organism TEXT, gene TEXT, pe INT"
            ");"
        )
        self.server.adaptor.execute(query,)
        query = (
            "CREATE INDEX spdi ON swissprot_defs(prot_id);"
        )
        self.server.adaptor.execute(query,)


    def load_swissprot_hits(self, data):
        self.load_data_into_table("swissprot_hits", data)


    def load_swissprot_defs(self, data):
        self.load_data_into_table("swissprot_defs", data)


    def create_diamond_refseq_match_id(self):
        query = (
            "CREATE TABLE diamond_refseq_match_id ( "
            "match_id INT, accession TEXT, organism TEXT, description TEXT, length INT, "
            "PRIMARY KEY(match_id));"
        )
        self.server.adaptor.execute(query,)

    def load_diamond_refseq_match_id(self, data):
        self.load_data_into_table("diamond_refseq_match_id", data)

    def create_diamond_refseq_match_id_indices(self):
        query = (
            "CREATE INDEX drmii ON diamond_refseq_match_id(match_id);"
        )
        self.server.adaptor.execute(query,)

    def create_refseq_hits_taxonomy(self):
        sql = (
            "CREATE TABLE IF NOT EXISTS refseq_hits_taxonomy(taxid INT, superkingdom INT, phylum INT, "
            "class INT, order_id INT, family INT, genus INT, specie INT, PRIMARY KEY(taxid));"
        )
        self.server.adaptor.execute(sql)

    def create_refseq_hits_taxonomy_indices(self):
        sql = (
            "CREATE INDEX rhti ON refseq_hits_taxonomy(taxid)"
        )
        self.server.adaptor.execute(sql)

    def create_taxonomy_mapping(self, hsh):
        sql = (
            "CREATE TABLE taxonomy_mapping(taxid INT, "
            "rank TEXT, value TEXT, PRIMARY KEY(taxid));"
        )
        self.server.adaptor.execute(sql)
        lst_values = []
        for key, (rank, value) in hsh.items():
            lst_values.append( (key, rank, value) )
        self.load_data_into_table("taxonomy_mapping", lst_values)

    def load_data_into_table(self, table, data):
        if len(data)==0:
            return

        fmt_string = ", ".join("?" for i in range(len(data[0])))
        sql_string = f"INSERT into {table} VALUES ({fmt_string});"
        self.server.adaptor.executemany(sql_string, data)

    def load_refseq_hits_taxonomy(self, data):
        self.load_data_into_table("refseq_hits_taxonomy", data)

    def get_all_taxids(self):
        query = (
            "SELECT DISTINCT taxid from diamond_refseq_match_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        taxids = []
        for line in results:
            taxids.append(int(line[0]))
        return taxids

    # Utility class to make user code more readable
    class Taxon:
        hsh_taxo_key = {
                "superkingdom" : (1, 2),
                "phylum" :       (3, 4),
                "class" :        (5, 6),
                "order" :        (7, 8),
                "family" :       (9, 10),
                "genus" :        (11, 12),
                "species" :      (13, 14)
        }

        def __init__(self, line):
            self.raw_result = line

        def phylum(self):
            idx = DB.Taxon.hsh_taxo_key["phylum"][0]
            return self.raw_result[idx]

        def taxid(self):
            return self.raw_result[0]

        def update_hash(self, curr_hash):
            for rank, (name_idx, taxid_idx) in DB.Taxon.hsh_taxo_key.items():
                if self.raw_result[taxid_idx] in curr_hash:
                    continue
                curr_hash[self.raw_result[taxid_idx]] = (rank, self.raw_result[name_idx])

        def get_all_taxids(self):
            all_taxids = [self.taxid()]
            for taxon_name, (idx_name, idx_taxid) in DB.Taxon.hsh_taxo_key.items():
                all_taxids.append(self.raw_result[idx_taxid])
            return all_taxids

    def get_linear_taxonomy(self, args, taxids):
        conn_refseq = sqlite3.connect(args["databases_dir"] + "/ncbi-taxonomy/linear_taxonomy.db")
        cursor = conn_refseq.cursor()

        query_string = ",".join("?" for i in taxids)
        query = (
            "SELECT tax_id, `superkingdom`, superkingdom_taxid, "
            " `phylum`, phylum_taxid, `class`, class_taxid, "
            " `order`, order_taxid, `family`, family_taxid, "
            " `genus`, genus_taxid, `species`, species_taxid "
            f"FROM ncbi_taxonomy WHERE tax_id IN ({query_string});"
        )
        results = cursor.execute(query, taxids).fetchall()
        return (DB.Taxon(line) for line in results)

    def get_accession_to_taxid(self, accession, params):
        # reuse the connection to the database if already open
        if self.conn_ncbi_taxonomy == None:
            dtb_path = params["databases_dir"] + "/ncbi-taxonomy/prot_accession2taxid.db"
            self.conn_ncbi_taxonomy = sqlite3.connect(dtb_path)

        cursor = self.conn_ncbi_taxonomy.cursor()
        query_string = ",".join(["\"" + a + "\"" for a in accession])
        query = (
            f"SELECT accession, taxid FROM accession2taxid "
            f"WHERE accession IN ({query_string});" 
        )
        results = cursor.execute(query,).fetchall()
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results

    # NOTE: need to check which indices are necessary to add to this table
    def load_cog_hits(self, data):
        sql = (
            "CREATE TABLE cog_hits (hsh INTEGER, cog_id INT, evalue FLOAT);"
            # " FOREIGN KEY(cog_id) REFERENCES cog_names(cog_id)); "
        )
        self.server.adaptor.execute(sql,)

        sql = "CREATE INDEX chi_hsh ON cog_hits(hsh);"
        self.load_data_into_table("cog_hits", data)

    def load_cog_fun_data(self, data):
        sql = (
            "CREATE TABLE cog_functions (function TEXT, description TEXT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("cog_functions", data)


    def load_cog_ref_data(self, data):
        sql = (
            "CREATE TABLE cog_names (cog_id INTEGER, function TEXT, description TEXT, "
            "PRIMARY KEY(cog_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("cog_names", data)


    def load_ko_pathway(self, data):
        sql = (
            "CREATE TABLE ko_pathway_def ("
            "pathway_id INTEGER, desc TEXT, "
            "PRIMARY KEY(pathway_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_pathway_def", data)
        sql = "CREATE INDEX kpd_i ON ko_pathway_def(pathway_id);"
        self.server.adaptor.execute(sql)


    def load_ko_module_classes(self, data):
        sql = (
            "CREATE TABLE ko_class( "
            "class_id INTEGER, descr TEXT, PRIMARY KEY(class_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_class", data)
        sql = "CREATE INDEX kc_i ON ko_class(class_id);"
        self.server.adaptor.execute(sql)


    def load_ko_module(self, data):
        sql = (
            "CREATE TABLE ko_module_def ("
            "module_id INTEGER, desc TEXT, definition TEXT, is_signature_module BOOL, class INT, subclass INT, "
            "PRIMARY KEY(module_id), FOREIGN KEY(class) REFERENCES ko_class(class_id), "
            "FOREIGN KEY(subclass) REFERENCES ko_class(class_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_module_def", data)
        sql = "CREATE INDEX kmd_i ON ko_module_def(module_id);"
        self.server.adaptor.execute(sql)

    def load_ko_to_pathway(self, data):
        sql = (
            "CREATE TABLE ko_to_pathway ("
            "ko_id INT, pathway_id INT, "
            "PRIMARY KEY(ko_id, pathway_id), "
            "FOREIGN KEY(ko_id) REFERENCES ko_def(ko_id), "
            "FOREIGN KEY(pathway_id) REFERENCES ko_pathway_def(pathway_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_to_pathway", data)
        sql = "CREATE INDEX ktpk_i ON ko_to_pathway(ko_id);"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX ktpp_i ON ko_to_pathway(pathway_id);"
        self.server.adaptor.execute(sql)

    def load_ko_to_module(self, data):
        sql = (
            "CREATE TABLE ko_to_module ("
            "ko_id INT, module_id INT, "
            "PRIMARY KEY(ko_id, module_id), "
            "FOREIGN KEY(ko_id) REFERENCES ko_def(ko_id), "
            "FOREIGN KEY(module_id) REFERENCES ko_module_def(module_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_to_module", data)
        sql = "CREATE INDEX ktmk_i ON ko_to_module(ko_id);"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX ktmm_i ON ko_to_module(module_id);"
        self.server.adaptor.execute(sql)


    # Note: EC to add separately?
    def load_ko_def(self, data):
        sql = (
            "CREATE TABLE ko_def ( "
            "ko_id INTEGER, descr TEXT, PRIMARY KEY(ko_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_def", data)
        sql = "CREATE INDEX kdko_i ON ko_def(ko_id);"
        self.server.adaptor.execute(sql)


    def load_ko_hits(self, data):
        sql = (
            "CREATE TABLE ko_hits ("
            "hsh INTEGER, ko_id INT, threshold FLOAT, score FLOAT, evalue FLOAT, "
            "PRIMARY KEY(hsh, ko_id),"
            "FOREIGN KEY(ko_id) REFERENCES ko_def(ko_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("ko_hits", data)
        sql = (
            "CREATE INDEX khi_i ON ko_hits(ko_id);"
        )
        self.server.adaptor.execute(sql)

    def create_pfam_def_table(self, entries):
        sql = (
            "CREATE TABLE pfam_table( "
            "pfam_id INTEGER, definition TEXT, "
            "PRIMARY KEY(pfam_id));"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("pfam_table", entries)
        sql = (
            "CREATE INDEX pfam_table_idx ON pfam_table(pfam_id);"
        )
        self.server.adaptor.execute(sql)

    def create_pfam_hits_table(self):
        sql = (
            "CREATE TABLE pfam_hits( "
            "hsh INTEGER, pfam_id INTEGER, start INTEGER, end INTEGER,"
            "PRIMARY KEY(hsh, pfam_id, start, end) "
            ");"
        )
        self.server.adaptor.execute(sql)
        sql = (
            "CREATE INDEX pfam_hits_idx ON pfam_hits(hsh);"
        )
        self.server.adaptor.execute(sql)
        sql = (
            "CREATE INDEX pfam_hits_pfam ON pfam_hits(pfam_id);"
        )
        self.server.adaptor.execute(sql)


    def get_all_modules_definition(self, allow_signature=False):
        where = ""
        if not allow_signature:
            where = "WHERE is_signature_module = 0"

        query = (
            "SELECT module_id, definition FROM ko_module_def "
            f"{where};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]


    def get_module_categories(self, module_ids=None):
        if module_ids!=None:
            selection_str = ",".join(str(mod_id) for mod_id in module_ids)
            selection = f"AND ko_module_def.module_id IN ({selection_str})"
        else:
            selection = ""

        query = (
            "SELECT class_id, descr "
            "FROM ko_module_def "
            "INNER JOIN ko_class ON class_id = class "
            f"WHERE is_signature_module = 0 {selection}"
            "GROUP BY (descr);"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]


    def get_module_sub_categories(self, module_ids=None):
        if module_ids!=None:
            selection_str = ",".join(str(mod_id) for mod_id in module_ids)
            selection = f"AND ko_module_def.module_id IN ({selection_str})"
        else:
            selection = ""

        query = (
            "SELECT class_id, descr "
            "FROM ko_module_def "
            "INNER JOIN ko_class ON class_id = subclass "
            f"WHERE is_signature_module = 0 {selection}"
            "GROUP BY (descr);"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]


    def get_ko_desc(self, ko_ids):
        if ko_ids is None:
            where = (
                "INNER JOIN ko_hits AS hit ON hit.ko_id=ko.ko_id "
                "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh=hit.hsh "
                "GROUP BY ko.ko_id"
            )
        else:
            entries = self.gen_placeholder_string(ko_ids)
            where = f"WHERE ko.ko_id IN ({entries})"

        query = (
            "SELECT ko.ko_id, ko.descr "
            "FROM ko_def as ko "
            f"{where};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    def get_pathways (self):

        query = (
            "SELECT ko_to_pathway.pathway_id, ko_pathway_def.desc "
            "FROM ko_hits "
            "JOIN ko_to_pathway ON ko_hits.ko_id = ko_to_pathway.ko_id "
            "JOIN ko_pathway_def ON ko_pathway_def.pathway_id = ko_to_pathway.pathway_id "
            "GROUP BY ko_pathway_def.desc;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        return [(line[0], line[1]) for line in results]



    def get_ko_pathways(self, ids, search_on="ko", as_df=False):
        if search_on!="ko" and search_on!="pathway":
            raise RuntimeError("Search term not supported: "+search_on)

        if not ids is None:
            entries = ",".join("?" for i in ids)

        if search_on=="ko":
            where = "ktp.ko_id"
        elif search_on=="pathway":
            where  = "ktp.pathway_id"

        query = (
            "SELECT ktp.ko_id, path.pathway_id, path.desc "
            "FROM ko_to_pathway AS ktp "
            "INNER JOIN ko_pathway_def AS path ON path.pathway_id = ktp.pathway_id "
            f"WHERE {where} IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ids)

        if as_df:
            return DB.to_pandas_frame(results, ["ko", "pathway", "description"])
        hsh_results = {}
        for line in results:
            ko_id = line[0]
            data = hsh_results.get(ko_id, [])
            data.append((line[1], line[2]))
            hsh_results[ko_id] = data
        return hsh_results


    def get_module_kos(self, module_id):
        query = (
            "SELECT ktm.ko_id "
            "FROM ko_to_module AS ktm "
            f"WHERE ktm.module_id = {module_id};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        lst_results = []
        for line in results:
            lst_results.append(line[0])
        return lst_results


    # compact: do not return the module description if true
    def get_ko_modules(self, ko_ids, as_pandas=False, compact=False):
        entries = ",".join("?" for i in ko_ids)
        if compact:
            supp_query = ""
            ids = ["ko_id", "module_id"]
        else:
            supp_query = ", mod.desc"
            ids = ["ko_id", "module_id", "desc"]

        query = (
            f"SELECT ktm.ko_id, mod.module_id {supp_query} "
            "FROM ko_to_module AS ktm "
            "INNER JOIN ko_module_def AS mod ON mod.module_id = ktm.module_id "
            f"WHERE ktm.ko_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)
        hsh_results = {}
        if as_pandas:
            return DB.to_pandas_frame(results, ids)
        for line in results:
            ko_id = line[0]
            data = hsh_results.get(ko_id, [])

            if compact:
                data.append(line[1])
            else:
                data.append((line[1], line[2]))
            hsh_results[ko_id] = data
        return hsh_results


    def get_module_completeness(self, taxids):
        plchd = self.gen_placeholder_string(ids)
        query = (
            "SELECT module_id, taxid "
            "FROM module_completeness "
            f"WHERE taxid IN ({plchd});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, taxids)
        return DB.to_pandas_frame(results, ["module_id", "taxid"])


    def create_module_completeness_table(self):
        query = (
            f"CREATE TABLE module_completeness ( "
            f"module_id INT, "
            f"taxid INT,"
            f"PRIMARY KEY(module_id, taxid));"
        )
        self.server.adaptor.execute(query)
        sql = "CREATE INDEX mdi_taxid ON module_completeness(taxid);"
        self.server.adaptor.execute(sql)


    def load_module_completeness(self, modules):
        if len(modules)==0:
            return
        self.load_data_into_table("module_completeness", modules)


    def get_ko_count_cat(self, subcategory=None, taxon_ids=None,
            subcategory_name=None, category=None, index=True):

        # TODO: will need to re-implement with a search_on syntax
        sel = " TRUE "
        join = ""

        args = []
        if taxon_ids != None:
            sel_str = ",".join("?" for _ in taxon_ids)
            sel = f" entry.taxon_id IN ({sel_str}) "
            args = taxon_ids

        if subcategory != None:
            where = f" module.subclass = ? AND is_signature_module = 0 AND {sel}"
            args = [subcategory] + args
        elif subcategory_name != None:
            join = "INNER JOIN ko_class AS class ON module.subclass = class.class_id "
            where = f" class.descr = ? AND is_signature_module = 0 AND {sel}"
            args = [subcategory_name] + args
        elif not category is None:
            join = "INNER JOIN ko_class AS class ON module.class = class.class_id "
            where = f" module.class = ? AND is_signature_module = 0 AND {sel}"
            args = [category]+args
        else:
            where = sel

        query = (
            "SELECT entry.taxon_id, module.module_id, ktm.ko_id, COUNT(*) "
            "FROM ko_module_def AS module "
            "INNER JOIN ko_to_module AS ktm ON module.module_id = ktm.module_id "
            "INNER JOIN ko_hits AS hit ON hit.ko_id = ktm.ko_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hit.hsh = hsh.hsh "
            "INNER JOIN seqfeature AS fet ON fet.seqfeature_id = hsh.seqid "
            "INNER JOIN bioentry AS entry ON fet.bioentry_id=entry.bioentry_id "
            f"{join}"
            f"WHERE {where}"
            "GROUP BY entry.taxon_id, module.module_id, ktm.ko_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, args)
        columns = ["taxon_id", "module_id", "KO", "count"]
        df = DB.to_pandas_frame(results, columns)
        if index == False:
            return df
        return df.set_index(["taxon_id", "module_id", "KO"])


    def get_modules_info(self, ids, search_on="module", as_pandas=False):
        if not search_on is None and search_on != "module" and \
                search_on != "category" and search_on != "subcategory":
            raise RuntimeError(f"Unsupported search term: {search_on}")

        if not ids is None and search_on is None:
            raise RuntimeError("Invalid combination of ids and search_on being None")

        if not ids is None:
            fmt = ",".join("?" for i in ids)

        where_term = ""
        if search_on=="module":
            where_term = f"module_id IN ({fmt})"
        elif search_on=="category":
            where_term = f"cat.class_id IN ({fmt})"
        elif search_on=="subcategory":
            where_term = f"subcat.class_id IN ({fmt})"
        elif search_on is None:
            where_term = "1"

        query = (
            "SELECT module_id, desc, definition, cat.descr, subcat.descr "
            "FROM ko_module_def AS def "
            "INNER JOIN ko_class AS subcat ON subcat.class_id = def.subclass "
            "INNER JOIN ko_class AS cat ON cat.class_id = def.class "
            f"WHERE is_signature_module = 0 AND {where_term};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ids)

        if as_pandas:
            return DB.to_pandas_frame(results, ["module_id", "descr", "definition", "cat", "subcat"])
        return [(line[0], line[1], line[2], line[3], line[4]) for line in results]


    def gen_ko_where_clause(self, search_on, entries):
        entries = self.gen_placeholder_string(entries)
        if search_on=="seqid":
            where_clause = f" hsh.seqid IN ({entries}) "
        elif search_on=="ko":
            where_clause = f" hit.ko_id IN ({entries}) "
        elif search_on=="taxid":
            where_clause = f" entry.taxon_id IN ({entries}) "
        else:
            raise RuntimeError(f"Searching on {search_on} is not supported")

        return where_clause

    def get_ko_hits(self, ids, search_on="taxid",
            keep_taxid=False, plasmids=None, indexing="taxid"):
        """
        Note: if search_on = ko, only the seqid are returned by default, if the keep_taxids
        is set, taxid are returned in an additional column.
        """

        where_clause = self.gen_ko_where_clause(search_on, ids)

        header = None
        select = "SELECT seqid.seqfeature_id, hit.ko_id "
        group_by = ""
        if indexing == "seqid":
            select =  "seqid.seqfeature_id, hit.ko_id "
            header = ["seqid", "ko"]
        elif indexing == "taxid":
            select = "entry.taxon_id, hit.ko_id, COUNT(*) "
            group_by = "GROUP BY entry.taxon_id, hit.ko_id "
            header = ["taxid", "ko", "count"]
        else:
            raise RuntimeError(f"Indexing by {search_on} not supported")

        if keep_taxid and not indexing=="taxid":
            select += ", entry.taxon_id "
            header.append("taxid")
        elif keep_taxid:
            raise RuntimeError(("Are you mocking me? Taxid is already returned! "
                "Please remove this pesky keep_taxid flag or use another search term!"))

        plasmid_join = ""
        if not plasmids is None:
            select += ", CAST(is_plasmid.value AS int) "
            subclause = self.gen_ko_where_clause(search_on, plasmids)
            where_clause = (
                    f"({where_clause} AND is_plasmid.value=0) "
                    f" OR ({subclause} AND is_plasmid.value=1)"
            )
            plasmid_join = (
                "INNER JOIN bioentry_qualifier_value AS is_plasmid ON "
                "  is_plasmid.bioentry_id=entry.bioentry_id "
                "INNER JOIN term AS plasmid_term ON plasmid_term.term_id=is_plasmid.term_id "
                "  AND plasmid_term.name=\"plasmid\""
            )
            if indexing!="seqid":
                group_by += ", CAST(is_plasmid.value AS int)"
            header.append("plasmid")

        query = (
            f"SELECT {select}"
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id=entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.seqid = seqid.seqfeature_id "
            "INNER JOIN ko_hits AS hit ON hit.hsh=hsh.hsh "
            f"{plasmid_join}"
            f"WHERE {where_clause} "
            f"{group_by};"
        )
        all_ids = ids
        if not plasmids is None:
            all_ids += plasmids
        results = self.server.adaptor.execute_and_fetchall(query, all_ids)

        df = DB.to_pandas_frame(results, header)
        if df.empty:
            return df

        if indexing=="taxid":
            index = ["taxid", "ko"]
            if not plasmids is None:
                index += ["plasmid"]
            df = df.set_index(index).unstack(level=0, fill_value=0)

            if not plasmids is None:
                return df.unstack(level=1, fill_value=0)
            else:
                df.columns = [col for col in df["count"].columns.values]
        elif indexing=="seqid":
            df = df.set_index(["seqid"])
        return df


    def get_ko_count(self, search_entries, keep_seqids=False, search_on="taxid", as_multi=True):
        if search_on=="ko_id":
            where_clause = "hit.ko_id"
        elif search_on=="taxid":
            where_clause = "entry.taxon_id"
        else:
            raise RuntimeError(f"Searching on {search_on} not supported, must be taxid or ko_id")

        if keep_seqids:
            keep_sel = " hsh.seqid, "
            keep_grp = " , hsh.seqid"
            ids = ["taxid", "KO", "seqid", "count"]
        else:
            keep_sel = ""
            keep_grp = ""
            ids = ["taxid", "KO", "count"]

        entries = self.gen_placeholder_string(search_entries)
        query = (
            f"SELECT entry.taxon_id, hit.ko_id, {keep_sel} count(*) "
            "FROM seqfeature AS feature "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON feature.seqfeature_id = hsh.seqid "
            "INNER JOIN ko_hits AS hit ON hit.hsh = hsh.hsh "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=feature.bioentry_id "
            f"WHERE {where_clause} IN ({entries})"
            f"GROUP BY entry.taxon_id, hit.ko_id {keep_grp};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, search_entries)
        if not as_multi:
            return DB.to_pandas_frame(results, ids)
        return DB.to_pandas_frame(results, ids).set_index(["taxid", "KO"])


    def get_seqids_for_ko(self, ko_ids, only_seqids=False):
        if only_seqids:
            selection = ""
        else:
            selection = ", hits.ko_id "

        entries = ",".join("?" for i in ko_ids)
        query = (
            f"SELECT hsh.seqid {selection}"
            "FROM ko_hits AS hits "
            "INNER JOIN sequence_hash_dictionnary as hsh ON hsh.hsh = hits.hsh "
            f"WHERE ko_id IN ({entries}) GROUP BY hsh.seqid;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ko_ids)

        if only_seqids:
            return [line[0] for line in results]
        else:
            hsh_results = {}
            for line in results:
                assert line[0] not in hsh_results
                hsh_results[line[0]] = line[1]
            return hsh_results


    def get_hsh_locus_to_seqfeature_id(self, only_CDS=False):
        filtering = ""
        if only_CDS:
            filtering = "INNER JOIN term as t4 ON t4.term_id=t1.type_term_id AND t4.name = \"CDS\" "

        query = (
            "SELECT t2.value, t1.seqfeature_id "
            "FROM seqfeature as t1 "
            f"{filtering}"
            "INNER JOIN seqfeature_qualifier_value AS t2 ON t1.seqfeature_id=t2.seqfeature_id "
            "INNER JOIN term as t3 on t3.term_id=t2.term_id AND t3.name=\"locus_tag\" "
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results


    def get_og_phylogeny(self, og):
        query = (
            "SELECT tree FROM gene_phylogeny WHERE orthogroup_id=?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, [og])
        if len(results)!=1:
            raise NoPhylogenyException
        if len(results[0][0]) == 0:
            raise NoPhylogenyException
        return results[0][0]


    def get_locus_to_og(self):
        query = (
            f"SELECT locus.value, og.value"
            f"FROM seqfeature_qualifier_value as og "
            f"INNER JOIN term AS term_og ON term_og.term_id=og.term_id AND term_og.name=\"orthogroup\" "
            f"INNER JOIN seqfeature_qualifier_value as locus ON locus.seqfeature_id=og.seqfeature_id "
            f"INNER JOIN term AS term_locus ON term_locus.term_id=locus.term_id "
            f"  AND term_locus.name=\"locus_tag\";"
        )
        query_results = self.server.adaptor.execute_and_fetchall(query,)
        locus_to_og = {}
        for locus, og in query_results:
            locus_to_og[locus] = int(og)
        return locus_to_og

    def create_new_og_matrix(self):
        query = (
            f"CREATE TABLE orthology_identity ( "
            f"orthogroup INT, "
            f"id_1 INT, id_2 INT, identity FLOAT(5), "
            f"PRIMARY KEY(orthogroup, id_1, id_2), "
            f"FOREIGN KEY(id_1) REFERENCES seqfeature(seqfeature_id)"
            f"FOREIGN KEY(id_2) REFERENCES seqfeature(seqfeature_id))"
        )
        self.server.adaptor.execute(query, )

    # To be removed
    # For each taxon, get the number of protein present in each of the existing
    # ortholog group
    # Returns a table of maps, organized as follows:
    #  group 0: hsh taxid -> number of proteins of taxid in group 0
    #  group 1: hsh taxid -> number of proteins of taxid in group 1
    #  ... 
    #  group N: hsh taxid -> number of proteins of taxid in group N
    def get_orthogroup_count_table(self):
        query = (
            "SELECT value, taxon_id, count(*) "
            "FROM seqfeature_qualifier_value as ortho_table "
            "INNER JOIN term ON ortho_table.term_id = term.term_id AND term.name=\"orthogroup\" "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=ortho_table.seqfeature_id "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=feature.bioentry_id "
            "GROUP BY value, taxon_id; "
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        size = max([int(group) for group, taxid, cnt in results])
        arr_group = [None]*(size+1)
        for line in results:
            group, taxid, cnt = int(line[0]), line[1], int(line[2])
            hsh_taxid_to_count = arr_group[group]
            if hsh_taxid_to_count==None:
                hsh_taxid_to_count = {}
                arr_group[group] = hsh_taxid_to_count
            hsh_taxid_to_count[taxid] = cnt
        return arr_group

    def load_og_matrix(self, matrix):
        for orthogroup_id, id_1, id_2, identity in matrix:
            query = (
                f"INSERT INTO orthology_identity "
                f"VALUES ({orthogroup_id}, {id_1}, {id_2}, {identity});"
            )
            self.server.adaptor.execute(query,)


    def create_BBH_phylogeny_table(self, data):
        sql = (
            "CREATE TABLE BBH_phylogeny (orthogroup_id INTEGER, tree TEXT, "
            " PRIMARY KEY(orthogroup_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("BBH_phylogeny", data)

    
    def get_refseq_phylogeny(self, og_id):
        query = (
            "SELECT tree FROM BBH_phylogeny WHERE orthogroup_id = ?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, [og_id])
        if len(results) != 1:
            raise RuntimeError("No refseq best hit tree")
        return results[0][0]


    def create_gene_phylogeny_table(self, data):
        sql = (
            "CREATE TABLE gene_phylogeny (orthogroup_id INTEGER, tree TEXT, "
            " is_core INTEGER, og_size INTEGER, num_genomes INTEGER, "
            "PRIMARY KEY(orthogroup_id)); "
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("gene_phylogeny", data)


    def create_og_matrix_indices(self):
        sql_1 = "CREATE INDEX oio ON orthology_identity(orthogroup);"
        self.server.adaptor.execute(sql_1,)


    def load_og_hits(self, lst_to_load):
        query = (
            "CREATE TABLE og_hits ( "
            "seqid INTEGER, orthogroup INTEGER, "
            "PRIMARY KEY(seqid), "
            "FOREIGN KEY(seqid) REFERENCES seqfeature(seqfeature_id));"
        )
        self.server.adaptor.execute(query,)
        query = "INSERT INTO og_hits VALUES(?, ?);"
        self.server.adaptor.executemany(query, lst_to_load)
        query = "CREATE INDEX og_og_idx ON og_hits(seqid);"
        self.server.adaptor.execute(query)
        query = "CREATE INDEX og_seqid_idx ON og_hits(orthogroup);"
        self.server.adaptor.execute(query)


    def taxon_ids(self):
        query = "SELECT taxon_id FROM taxon";
        result = self.server.adaptor.execute_and_fetchall(query,)
        arr_taxon_ids = []
        for line in result:
            arr_taxon_ids.append(int(line[0]))
        return arr_taxon_ids


    def set_status_in_config_table(self, status_name, status_val):
        sql = f"update biodb_config set status={status_val} where name={quote(status_name)};"
        self.server.adaptor.execute(sql,)


    def get_config_table(self, ret_mandatory=False):
        sql = "SELECT * from biodb_config;"
        values = self.server.adaptor.execute_and_fetchall(sql)
        if ret_mandatory:
            return {val[0]: (val[1]=="mandatory", val[2]) for val in values}
        else:
            return {val[0]: val[2] for val in values}


    def create_biosql_database(self, args):
        self.server.new_database(self.db_name)

    # NOTE:
    # As the primary key in sequence_hash_dictionnary is a tuple of hsh and seqid, we would 
    # need to have an entry in diamond_refseq for every possible combination of (hsh, seqid) to 
    # be able to use FOREIGN KEY. To avoid that and spare some disk space, foreign keys were not
    # used..
    def create_refseq_hits_table(self):
        sql = (
            f"CREATE TABLE IF NOT EXISTS diamond_refseq(hit_count INT, seq_hash INTEGER, "
            f"sseqid INT, pident FLOAT, length INT, mismatch INT, "
            f"gapopen INT, qstart INT, qend INT, sstart INT, "
            f"send INT, evalue FLOAT, bitscore FLOAT, "
            "FOREIGN KEY(sseqid) REFERENCES diamond_refseq_match_id(match_id), "
            f"PRIMARY KEY (hit_count, seq_hash));"
        )
        self.server.adaptor.execute(sql,)

    def load_refseq_hits(self, data):
        self.load_data_into_table("diamond_refseq", data)

    def create_refseq_hits_indices(self):
        sql = (
            "CREATE INDEX dri ON diamond_refseq(seq_hash);"
        )
        self.server.adaptor.execute(sql,)

    def n_accession_to_n_tRNA(self):
        query = (
            "SELECT accession, count(*) FROM seqfeature AS seq "
            " INNER JOIN term AS t ON t.term_id = seq.type_term_id "
            " INNER JOIN bioentry AS entry ON seq.bioentry_id=entry.bioentry_id "
            " AND t.name=\"tRNA\" GROUP BY accession; "
        )
        n_tRNA_results = self.server.adaptor.execute_and_fetchall(query)
        results = {}
        for line in n_tRNA_results:
            results[line[0]] = int(line[1])
        return results

    
    def get_locus_to_genomes(self, locus_lst):
        # quick and dirty: may need to regroup it with another function
        values = ",".join("?" for _ in locus_lst)
        query = (
            "SELECT locus_tag.value, name.name "
            "FROM seqfeature_qualifier_value AS locus_tag "
            "INNER JOIN term AS locus_term ON locus_term.term_id=locus_tag.term_id "
            " AND locus_term.name=\"locus_tag\" "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=locus_tag.seqfeature_id "
            "INNER JOIN bioentry AS entry ON feature.bioentry_id=entry.bioentry_id "
            "INNER JOIN taxon_name AS name ON entry.taxon_id = name.taxon_id "
            f"WHERE locus_tag.value IN ({values});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, locus_lst)
        return {locus: description for locus, description in results}

    
    def get_term_id(self, term, create_if_absent=False, ontology="Annotation Tags"):
        query = f"SELECT term_id FROM term WHERE name=\"{term}\";"
        result = self.server.adaptor.execute_and_fetchall(query)
        gc_term_id = None
        if len(result) == 0 and create_if_absent:
            sql1 = (
                "SELECT ontology_id FROM ontology "
                f"WHERE name=\"{ontology}\";"
            )
            ontology_id = self.server.adaptor.execute_and_fetchall(sql1)[0][0]
            sql = (
                f"INSERT INTO term (name, ontology_id) VALUES (\"{term}\", {ontology_id});"
            )
            self.server.adaptor.execute(sql)
            gc_term_id = self.server.adaptor.cursor.lastrowid
        else:
            gc_term_id = result[0][0]
        return gc_term_id


    def get_genomes_description(self, lst_plasmids=True):
        """
        Returns the description of the genome as it has been read from the genbank
        files, indexed by taxon_id. The output also contains a flag has_plasmid
        indicating whether the genome contains a plasmid or not, if the lst_plasmid flag
        has been set.
        """

        has_plasmid_query = (
            "SELECT * "
            "FROM bioentry_qualifier_value AS has_plasmid "
            "INNER JOIN term AS pls_term ON pls_term.term_id=has_plasmid.term_id "
            " AND pls_term.name=\"plasmid\" AND has_plasmid.value=1 " 
            "INNER JOIN bioentry AS plasmid ON has_plasmid.bioentry_id=plasmid.bioentry_id "
            "WHERE plasmid.taxon_id=entry.taxon_id"
        )
        query = (
            "SELECT entry.taxon_id, txn_name.name,  "
            f" CASE WHEN EXISTS ({has_plasmid_query}) THEN 1 ELSE 0 END "
            "FROM bioentry AS entry "
            "INNER JOIN bioentry_qualifier_value AS orga ON entry.bioentry_id=orga.bioentry_id " 
            "INNER JOIN taxon_name as txn_name ON entry.taxon_id=txn_name.taxon_id "
            "INNER JOIN term AS orga_term ON orga.term_id=orga_term.term_id "
            " AND orga_term.name=\"organism\" "
            "GROUP BY entry.taxon_id;" 
        )
        descr = self.server.adaptor.execute_and_fetchall(query)
        columns = ["taxon_id", "description", "has_plasmid"]
        return DB.to_pandas_frame(descr, columns).set_index(["taxon_id"])


    def get_genomes_infos(self):
        """
        Note: for efficiency sake, it would be possible to order the different
         tables by taxon_id, walk through the table using zip in an iteator
         and pass the iterator to pandas.

         However, given the small size of the dataset, it's probably not worth
         the implementation effort.
        """

        query = (
            "SELECT entry.taxon_id, entry.taxon_id, COUNT(*) "
            " FROM seqfeature AS seq "
            " INNER JOIN term AS cds ON cds.term_id = seq.type_term_id AND cds.name=\"CDS\" "
            " INNER JOIN bioentry AS entry ON entry.bioentry_id = seq.bioentry_id " 
            " GROUP BY entry.taxon_id;"
        )
        n_prot_results = self.server.adaptor.execute_and_fetchall(query)

        query = (
            "SELECT entry.taxon_id, COUNT(*) "
            "FROM bioentry AS entry "
            "GROUP BY entry.taxon_id;"
        )
        n_contigs = self.server.adaptor.execute_and_fetchall(query)

        cols = ["taxon_id", "completeness", "contamination", "gc", "length", "coding_density"]
        query = (
            "SELECT taxon_id, completeness, contamination, gc, length, coding_density "
            " from genome_summary;"
        )
        all_other_results = self.server.adaptor.execute_and_fetchall(query)

        df_n_prot = DB.to_pandas_frame(n_prot_results, ["taxon_id", "id", "n_prot"]).set_index("taxon_id")
        df_n_contigs = DB.to_pandas_frame(n_contigs, ["taxon_id", "n_contigs"]).set_index("taxon_id")
        df_stats = DB.to_pandas_frame(all_other_results, cols).set_index("taxon_id")
        return df_n_prot.join(df_n_contigs).join(df_stats)


    def get_contigs_to_seqid (self, genome):
        query = (
            "SELECT bioentry.name, ids.seqfeature_id "
            "FROM seqfeature as ids "
            "JOIN bioentry ON bioentry.bioentry_id=ids.bioentry_id "
            "WHERE bioentry.taxon_id=?;"
        )
        contigs_seqids = self.server.adaptor.execute_and_fetchall(query, [genome])
        df_contigs_seqid = DB.to_pandas_frame(contigs_seqids, ["contig", "seqid"])
        return df_contigs_seqid.set_index("seqid")

    
    def get_accession_to_entry(self):
        query = (
            "SELECT accession, bioentry_id FROM bioentry;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = line[1]
        return hsh_results

    # Return the number of proteins, indexed either by
    # - accession
    # - bioentry_id
    def n_CDS(self, indexing="bioentry", taxons=None, bioentries=None):

        if taxons!=None and bioentries!=None:
            raise RuntimeError("Specifying both taxons and bioentries not supported")

        # This needs some testing. Also needs query sanitizing
        filtering = ""
        if taxons != None:
            j = ",".join([str(taxon) for taxon in taxons])
            filtering = (
                " INNER JOIN taxon AS t ON t.taxon_id=entry.taxon_id "
                f" AND ncbi_taxon_id IN ({j});"
            )
        elif bioentries != None:
            j = ",".join([str(entry) for entry in bioentries])
            filtering = f"WHERE entry.bioentry_id IN {j}"

        index = None
        if indexing == "bioentry":
            index = "entry.bioentry_id"
        elif indexing == "accession":
            index = "accession"
        elif indexing == "ncbi":
            index = "ncbi_taxon_id"
        else:
            raise RuntimeError(f"Unsupported indexing method: {indexing}")

        query = (
            f"SELECT {index}, COUNT(*) FROM seqfeature AS prot "
            " INNER JOIN term AS t ON t.term_id = prot.type_term_id AND t.name=\"CDS\" "
            " INNER JOIN bioentry AS entry ON entry.bioentry_id = prot.bioentry_id "
            f"{filtering}"
            f"GROUP BY {index};"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for index, str_count in results:
            hsh_results[index] = int(str_count)
        return hsh_results

    def get_cog_code_description(self):
        query = (
            "SELECT function, description FROM cog_functions;"
        )
        results = self.server.adaptor.execute_and_fetchall(query)
        hsh_results = {}
        for line in results:
            hsh_results[line[0].strip()] = line[1].strip()
        return hsh_results


    # Note: need to check whether this is a more correct way to proceed
    # as my Rhabdo genomes currently all have the same taxon id, the following
    # code wouldn't work
    # sql_n_genomes = 'select count(*) from (select distinct taxon_id from bioentry t1 ' \
    # ' inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id where t2.name="%s") A;' % biodb
    def get_n_genomes(self):
        query = "SELECT COUNT(*) FROM bioentry GROUP BY taxon_id;"
        result = self.server.adaptor.execute_and_fetchall(query)
        return result[0][0]


    def load_reference_phylogeny(self, tree):
        sql = "CREATE TABLE IF NOT EXISTS reference_phylogeny (tree TEXT);"
        self.server.adaptor.execute(sql,)

        sql = f"INSERT INTO reference_phylogeny VALUES (\"{tree}\");"
        self.server.adaptor.execute(sql,)

    def get_reference_phylogeny(self):
        sql = "SELECT tree FROM reference_phylogeny;"
        values = self.server.adaptor.execute_and_fetchall(sql)
        return values[0][0]


    def get_bioentry_id_for_record(self, record):
        locus_tag = record["gene"]["locus_tag"]
        sql = (
            "SELECT bioentry_id FROM seqfeature_qualifier_value AS tag"
            " INNER JOIN seqfeature AS seq ON seq.seqfeature_id = tag.seqfeature_id "
            f" where value = \"{locus_tag}\";" 
        )
        results = self.server.adaptor.execute_and_fetchall(sql)
        return results[0][0]



    def load_filenames(self, data):
        sql = (
            "CREATE TABLE filenames (taxon_id INTEGER, filename TEXT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("filenames", data)


    def load_genomes_info(self, data):
        sql = (
            "CREATE TABLE genome_summary (taxon_id INTEGER, completeness FLOAT, "
            " contamination FLOAT, gc INTEGER, length INTEGER, coding_density FLOAT);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("genome_summary", data)


    def get_cog_summaries(self, cog_ids, only_cog_desc=False, as_df=False):
        if cog_ids is None:
            where = (
                "INNER JOIN cog_hits AS hits ON hits.cog_id=names.cog_id "
                "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh=hits.hsh "
                "GROUP BY names.cog_id"
            )
        else:
            ids = ",".join(["?"] * len(cog_ids))
            where = f"WHERE names.cog_id IN ({ids})"

        query = (
            "SELECT names.cog_id, function, description "
            "FROM cog_names AS names "
            f"{where};"
        )

        results = self.server.adaptor.execute_and_fetchall(query, cog_ids)
        if only_cog_desc and not as_df:
            hsh_results = {}
            for line in results:
                hsh_results[line[0]] = (line[1], line[2])
            return hsh_results
        elif only_cog_desc:
            return DB.to_pandas_frame(results, ["cog", "function", "description"]).set_index(["cog"])

        funcs = "SELECT function, description FROM cog_functions;"
        functions = self.server.adaptor.execute_and_fetchall(funcs)
        hsh_func_to_description = {}
        for function, description in functions:
            hsh_func_to_description[function] = description

        hsh_results = { }
        for cog_id, function, cog_description in results:
            hsh_results[cog_id] = []
            for i in range(0, len(function)):
                func = function[i]
                func_descr = hsh_func_to_description[func]
                hsh_results[cog_id].append((func, func_descr, cog_description))
        return hsh_results


    def get_CDS_from_locus_tag(self, locus_tag):
        # NOTE: I did not add a join to filter on locus_tag,
        # it may be worth it performance-wise to pre-filter the 
        # database entries if this query were to become an issue.

        query = (
            "SELECT seqfeature_id "
            "FROM seqfeature_qualifier_values AS locus_tag "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=locus_tag.seqfeature_id "
            "INNER JOIN term AS t ON feature.type_term_id=t.term_id AND t.name = \"CDS\" "
            "WHERE locus_tag.value=?;"
        )
        ret = self.server.execute_and_fetchall(query, locus_tag)
        if ret==None or len(ret)==0:
            return None
        return ret[0][0]


    def get_seqid(self, locus_tag, feature_type=False):
        add_type = ""
        sel = "AND cds_term.name=\"CDS\""
        if feature_type:
            add_type = ", cds_term.name"
            sel = (
                "AND (cds_term.name=\"CDS\" OR cds_term.name=\"tRNA\" "
                " OR cds_term.name=\"tmRNA\" OR cds_term.name=\"rRNA\") "
            )
        is_pseudo = (
            ", CASE WHEN EXISTS ("
            "SELECT NULL FROM seqfeature_qualifier_value AS pseudo "
            "INNER JOIN term AS pseudo_term ON pseudo.term_id=pseudo_term.term_id "
            "  AND (pseudo_term.name = \"pseudo\" OR pseudo_term.name = \"pseudogene\") "
            "WHERE pseudo.seqfeature_id = locus_tag.seqfeature_id"
            ") THEN 1 ELSE 0 END "
        )

        query = (
            f"SELECT locus_tag.seqfeature_id {add_type} {is_pseudo}"
            "FROM seqfeature_qualifier_value AS locus_tag "
            "INNER JOIN seqfeature AS cds "
            " ON cds.seqfeature_id=locus_tag.seqfeature_id "
            "INNER JOIN term AS cds_term ON cds_term.term_id=cds.type_term_id "
            f" {sel}"
            "WHERE locus_tag.value = ?;"
        )

        values =  self.server.adaptor.execute_and_fetchall(query, [locus_tag,])
        if len(values)==0:
            raise RuntimeError("No such entry")
        return values[0]


    def get_DNA_sequence(self, bioentry_id, alphabet="dna"):
        query = (
            "SELECT seq FROM biosequence WHERE bioentry_id=? AND alphabet = ?;"
        )

        results = self.server.adaptor.execute_and_fetchall(query, [bioentry_id, alphabet])
        return Seq(results[0][0])


    # Note: ordering by seqid makes it faster to assemble informations
    # from several queries if the index is the same.
    def get_gene_loc(self, seqids, search_on="seqid", as_hash=True):
        seqids_query = ",".join(["?"] * len(seqids))

        query = (
            "SELECT seqfeature_id, strand, start_pos, end_pos "
            f"FROM location WHERE seqfeature_id IN ({seqids_query});"
        )

        results = self.server.adaptor.execute_and_fetchall(query, seqids)
        if not as_hash:
            return DB.to_pandas_frame(results, ["seqid", "strand", "start", "end"])

        hsh_results = {}
        for line in results:
            seqid = line[0]
            strand = line[1]
            start = line[2]
            end = line[3]
            hsh_results[seqid] = [strand, start, end]
        return hsh_results


    def get_CDS_type(self, ids):
        plchd = self.gen_placeholder_string(ids)

        is_pseudo_query = (
            "SELECT NULL FROM seqfeature_qualifier_value AS pseudo "
            "INNER JOIN term AS pseudo_term ON pseudo.term_id=pseudo_term.term_id "
            "  AND (pseudo_term.name = \"pseudo\" OR pseudo_term.name = \"pseudogene\") "
            "WHERE pseudo.seqfeature_id = fet.seqfeature_id"
        )

        query = (
            "SELECT fet.seqfeature_id, type.name,  "
            f" CASE WHEN EXISTS({is_pseudo_query}) THEN 1 ELSE 0 END "
            "FROM seqfeature AS fet "
            "INNER JOIN term AS type ON type.term_id=fet.type_term_id "
            f"WHERE fet.seqfeature_id IN ({plchd});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ids)
        return DB.to_pandas_frame(results, ["seqid", "type", "is_pseudo"]).set_index(["seqid"])


    # Returns the seqid, locus tag, protein id, product and gene for a given
    # list of seqids
    def get_proteins_info(self, ids, search_on="seqid", as_df=False,
            to_return=None, inc_non_CDS=False, inc_pseudo=False):
        """
        Args:
         - to_return: accepts a list of entries to return for each seqid.
                Allowed: locus_tag, protein_id, gene and product
        """
        entries = self.gen_placeholder_string(ids)

        if to_return is None:
            term_names = ["locus_tag", "protein_id", "gene", "product"]
        else:
            term_names = to_return
            for term_name in to_return:
                if term_name not in ["locus_tag", "protein_id", "gene", "product"]:
                    raise RuntimeError(f"Value not support {term_name}")

        term_names_query = ",".join(f"\"{name}\"" for name in term_names)

        add_cond = ""
        if inc_non_CDS:
            add_cond = (
                "OR cds_term.name=\"tRNA\" OR cds_term.name=\"tmRNA\" "
                "OR cds_term.name=\"rRNA\""
            )

        no_pseudo = ""
        if not inc_pseudo:
            no_pseudo = (
                "AND NOT EXISTS ("
                " SELECT NULL FROM seqfeature_qualifier_value AS cds"
                " INNER JOIN term AS pseudo_term ON cds.term_id=pseudo_term.term_id "
                "  AND (pseudo_term.name=\"pseudogene\" OR pseudo_term.name=\"pseudo\") "
                " WHERE cds.seqfeature_id=seq.seqfeature_id "
                ")"
            )

        where = ""
        if search_on=="taxid":
            sel = (
                "INNER JOIN bioentry AS entry ON seq.bioentry_id=entry.bioentry_id "
            )
            where = f" entry.taxon_id IN ({entries})"
        elif search_on=="seqid":
            sel = ""
            where = f"v.seqfeature_id IN ({entries})"
        elif search_on=="locus_tag":
            sel = (
                "INNER JOIN seqfeature_qualifier_value AS locus_tag ON locus_tag.seqfeature_id=v.seqfeature_id "
                "INNER JOIN term as locus_tag_term ON locus_tag.term_id=locus_tag_term.term_id "
                " AND locus_tag_term.name=\"locus_tag\" "
            )
            where = f" locus_tag.value IN ({entries}) "
        else:
            raise RuntimeError("Searching is only possible on seqid and taxid")

        query = (
            "SELECT v.seqfeature_id, t.name, v.value "
            "FROM seqfeature_qualifier_value AS v "
            "INNER JOIN term AS t ON t.term_id = v.term_id "
            "INNER JOIN seqfeature AS seq ON seq.seqfeature_id = v.seqfeature_id "
            "INNER JOIN term AS cds_term ON seq.type_term_id=cds_term.term_id "
            f" AND (cds_term.name=\"CDS\" {add_cond})"
            f"{sel}"
            f"WHERE {where} AND t.name IN ({term_names_query}) {no_pseudo};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, ids)

        if as_df:
            df = DB.to_pandas_frame(results, ["seqid", "name", "value"])
            if df.empty:
                return df
            df = df.set_index(["seqid", "name"]).unstack(level="name")
            df.columns = [col for col in df.value.columns.values]
            return df

        hsh_results = {}
        for line in results:
            seqid = line[0]
            term_name = line[1]
            value = line[2]
            if seqid not in hsh_results:
                hsh_results[seqid] = [None]*len(term_names)
            hsh_results[seqid][term_names.index(term_name)] = value
        return hsh_results


    def get_organism(self, ids, as_df=False, id_type="seqid", as_taxid=False):
        seqids_query = ",".join(["?"] * len(ids))

        val = ""
        if as_taxid:
            val = "entry.taxon_id"
        else:
            val = "organism.value"

        if id_type == "seqid":
            query = (
                f"SELECT feature.seqfeature_id, {val} "
                "FROM seqfeature AS feature "
                "INNER JOIN bioentry AS entry ON feature.bioentry_id = entry.bioentry_id "
                "INNER JOIN bioentry_qualifier_value AS organism "
                "    ON organism.bioentry_id = entry.bioentry_id "
                "INNER JOIN term AS organism_term ON organism.term_id = organism_term.term_id "
                " AND organism_term.name = \"organism\" "
                f"WHERE feature.seqfeature_id IN ({seqids_query});"
            )
        elif id_type == "bioentry":
            query = (
                f"SELECT organism.bioentry_id, {val} "
                "FROM bioentry_qualifier_value AS organism "
                "INNER JOIN term AS organism_term ON organism.term_id = organism_term.term_id "
                " AND organism_term.name = \"organism\" "
                f"WHERE organism.bioentry_id IN ({seqids_query});"
            )
        else:
            raise RuntimeError("Id_type should be either seqid or bioentry")

        results = self.server.adaptor.execute_and_fetchall(query, ids)
        if as_df:
            col_name = "taxid" if as_taxid else "organism"
            df = DB.to_pandas_frame(results, columns=["seqid", "taxid"])
            return df.set_index("seqid")

        hsh_results = {}
        for line in results:
            seqid = line[0]
            organism = line[1]
            hsh_results[seqid] = organism
        return hsh_results


    # Note:
    # First extracting the data in memory and creating the dataframe
    # from it is much faster than iterating over results and adding
    # elements separately.
    #
    # NOTE: may be interesting to use int8/16 whenever possible 
    # to spare memory.
    def to_pandas_frame(db_results, columns, types=None):
        return pd.DataFrame(db_results, columns=columns)


    def get_bioentries_in_taxon(self, bioentries=None):
        # NOTE: need to write the code for the base where bioentries is None
        # -> returns all the entries
        query_str = ",".join("?" for entry in bioentries)
        query = (
            "SELECT two.bioentry_id, one.taxon_id, one.bioentry_id "
            "FROM bioentry AS one "
            "INNER JOIN bioentry AS two ON one.taxon_id = two.taxon_id "
            f"WHERE one.bioentry_id IN ({query_str});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, bioentries)
        return DB.to_pandas_frame(results, ["bioentry", "taxon", "ref_genome_bioentry"])


    def get_og_identity(self, og=None, ref_seqid=None, pairs=None):
        """
         Hackish... should be refactored at some point
        """
        if og is None and ref_seqid is None and pairs is None:
            raise RuntimeError("Need at least og or ref_seqid specified")

        values = []
        conjonctions = []
        og_query = ""
        conj_op = " AND "
        if not og is None and pairs is None:
            values.append(og)
            conjonctions.append("orthogroup = ?")

        if not ref_seqid is None:
            conjonctions.append("(id_1 = ? OR id_2 = ?)")
            values.append(ref_seqid)
            values.append(ref_seqid)
        
        if not pairs is None:
            if og is None:
                raise RuntimeError("")
            elif not ref_seqid is None:
                raise RuntimeError("")
            for curr_og, (id1, id2) in zip(og, pairs):
                conj = f"(orthogroup=? AND ((id_1=? AND id_2=?) OR (id_1=? AND id_2=?)))"
                conjonctions.append(conj)
                values.extend([curr_og, id1, id2, id2, id1])
            conj_op = " OR "

        where = f"{conj_op}".join(conjonctions)
        query = (
            "SELECT id_1, id_2, identity "
            "FROM orthology_identity "
            f"WHERE {where};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, values)

        if not pairs is None:
            return DB.to_pandas_frame(results, ["seqid_x", "seqid_y", "identity"])
        filtered_values = []

        for id_1, id_2, identity in results:
            if id_1==ref_seqid:
                filtered_values.append((id_2, identity))
            else:
                filtered_values.append((id_1, identity))
        df = DB.to_pandas_frame(filtered_values, ["seqid", "identity"])
        return df.set_index(["seqid"])


    def get_og_count(self, lookup_terms, search_on="taxid", plasmids=None, keep_taxid=False):
        """
        This function returns a pandas dataframe containing the orthogroup
        count for a set of taxon_ids. The user can differentiate between the chromosome
        and the plasmid.

        plasmids: if set to None, the function will not differentiate between
            plasmids and chromosomes. Otherwise, should contain a list of taxids whose plasmids
            will be taken into account in the search.
        lookup_term: the terms that will be search. Either a list of taxids, of orthogroup or of seqids.
        search_on: can be either orthogroup, seqid or taxid. Specifies what lookup_term is.
        keep_taxid: for queries using the seqid lookup term, also return the taxid in the results

        The index of the table depends on whether the diff_plasmid flag was set. 
        If the flag is set, MultiIndex(taxid, is_plasmid). If it is not set Index(taxid)
        """

        if not plasmids is None and search_on!="taxid":
            raise RuntimeError("Plasmid search only supported for taxid")

        entries = self.gen_placeholder_string(lookup_terms)

        add_plasmid = ""
        if not plasmids is None:
            add_plasmid = "CAST(is_plasmid.value AS int), "

        select   = f"SELECT entry.taxon_id, {add_plasmid} orthogroup, COUNT(*) "
        group_by = f"GROUP BY entry.taxon_id, {add_plasmid} orthogroup"
        where_clause = None
        add_plasmid_join = ""
        if not plasmids is None:
            plasmids_placeholder = self.gen_placeholder_string(plasmids)
            add_plasmid = "is_plasmid.value"
            add_plasmid_join = (
                "INNER JOIN bioentry_qualifier_value AS is_plasmid ON "
                "  is_plasmid.bioentry_id=entry.bioentry_id "
                "INNER JOIN term AS plasmid_term ON plasmid_term.term_id=is_plasmid.term_id "
                "  AND plasmid_term.name=\"plasmid\""
            )
            where_clause = (
                f" (entry.taxon_id IN ({entries}) AND is_plasmid.value=0) "
                f"OR (entry.taxon_id IN ({plasmids_placeholder}) AND is_plasmid.value=1)"
            )
        elif search_on=="taxid":
            where_clause = f"entry.taxon_id IN ({entries})"
        elif search_on=="orthogroup":
            where_clause = f"og.orthogroup IN ({entries})"
        elif search_on=="seqid":
            where_clause = f"feature.seqfeature_id IN ({entries})"
            select       = "SELECT feature.seqfeature_id, og.orthogroup "
            group_by     = ""
            if keep_taxid:
                select += ", entry.taxon_id "
        else:
            raise RuntimeError(f"Unsupported search {search_on}, must use orthogroup, bioentry or seqid")

        query = (
            f"{select} "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS feature ON entry.bioentry_id = feature.bioentry_id "
            "INNER JOIN og_hits AS og ON og.seqid = feature.seqfeature_id "
            f"{add_plasmid_join} "
            f"WHERE {where_clause} "
            f"{group_by};"
        )
        all_terms = lookup_terms
        if not plasmids is None:
            all_terms = lookup_terms + plasmids
        results = self.server.adaptor.execute_and_fetchall(query, all_terms)

        header = None
        if not plasmids is None:
            header = ["taxid", "plasmid", "orthogroup", "count"]
        elif search_on=="taxid" or search_on=="orthogroup":
            header = ["taxid", "orthogroup", "count"]
        elif search_on=="seqid":
            header = ["seqid", "orthogroup"]
            if keep_taxid:
                header.append("taxid")

        df = DB.to_pandas_frame(results, header)
        if len(df.index) == 0:
            return df
        if not plasmids is None:
            df = df.set_index(["taxid", "plasmid", "orthogroup"]).unstack(level=0, fill_value=0)
            return df.unstack(level=0, fill_value=0)
        elif search_on=="taxid" or search_on=="orthogroup":
            df = df.set_index(["taxid", "orthogroup"]).unstack(level=0, fill_value=0)
            df.columns = [col for col in df["count"].columns.values]
        elif search_on=="seqid":
            df = df.set_index(["seqid"])
        return df

    def get_translation(self, seqid):
        query = (
            "SELECT translation.value  "
            "FROM seqfeature_qualifier_value AS translation "
            "INNER JOIN term AS transl_term ON transl_term.term_id=translation.term_id "
            " AND transl_term.name=\"translation\" "
            f"WHERE translation.seqfeature_id = ?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, [seqid])
        if results==None or len(results)==0:
            raise RuntimeError("No translation")
        return results[0][0]

    def get_genes_from_og(self, orthogroups, taxon_ids=None, terms=["gene", "product"]):
        for term in terms:
            if term not in ["length", "gene", "product", "locus_tag"]:
                raise RuntimeError(f"Term not supported: {term}")

        og_entries = self.gen_placeholder_string(orthogroups)
        query_args = orthogroups

        db_terms = [t for t in terms if t!="length"]
        if "length" in terms:
            db_terms.append("translation")
        sel_terms = ",".join("\"" + f"{i}" + "\"" for i in db_terms)
        join, sel = "", ""
        if taxon_ids != None:
            taxon_id_query = self.gen_placeholder_string(taxon_ids)
            join = (
                "INNER JOIN seqfeature AS seq ON og.seqid=seq.seqfeature_id "
                "INNER JOIN bioentry AS entry ON seq.bioentry_id=entry.bioentry_id "
            )
            sel = f"AND entry.taxon_id IN ({taxon_id_query})"
            query_args = orthogroups+taxon_ids
        
        query = (
            f"SELECT feature.seqfeature_id, og.orthogroup, t.name, "
            "   CASE WHEN t.name==\"translation\" THEN "
            "   LENGTH(value) ELSE value END "
            "FROM og_hits AS og "
            "INNER JOIN seqfeature_qualifier_value AS feature ON og.seqid = feature.seqfeature_id "
            f"INNER JOIN term AS t ON t.term_id=feature.term_id AND t.name IN ({sel_terms})"
            f"{join}"
            f"WHERE og.orthogroup IN ({og_entries}) {sel};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, query_args)
        df = DB.to_pandas_frame(results, columns=["seqid", "orthogroup", "term", "value"])
        df = df.set_index(["orthogroup", "seqid", "term"]).unstack("term", fill_value=None)
        if "value" in df.columns:
            df.columns = ["length" if col=="translation" else col for col in df["value"].columns.values]

        df = df.reset_index(level=0)
        for t in terms:
            if t not in df.columns:
                df[t] = None
        return df


    def get_cog_counts_per_category(self, taxon_ids):
        entries = self.gen_placeholder_string(taxon_ids)
        query = (
            "SELECT entry.taxon_id, cog.function "
            "FROM seqfeature AS feature "
            "INNER JOIN bioentry AS entry ON feature.bioentry_id=entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh on hsh.seqid = feature.seqfeature_id "
            "INNER JOIN cog_hits AS hit ON hit.hsh = hsh.hsh "
            "INNER JOIN cog_names AS cog ON cog.cog_id = hit.cog_id "
            f"WHERE entry.taxon_id IN ({entries});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, taxon_ids)
        hsh_results = {}
        for line in results:
            entry_id, func = line
            if entry_id not in hsh_results:
                hsh_results[entry_id] = {func : 1}
            else:
                cnt = hsh_results[entry_id].get(func, 0)
                hsh_results[entry_id][func] = cnt+1
        return hsh_results
    
    def get_bioentry_qualifiers(self, bioentry_id):
        query = (
            "SELECT t.name, value.value "
            "FROM bioentry_qualifier_value AS value "
            "INNER JOIN term AS t ON t.term_id=value.term_id "
            "WHERE value.bioentry_id = ?;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, [bioentry_id])
        header = ["term", "value"]
        return DB.to_pandas_frame(results, header)
    
    def get_bioentry_list(self, search_term, search_on="taxid", min_bioentry_length=None):
        if search_on=="taxid":
            where = f"t1.taxon_id = {self.placeholder}"
            join = ""
        elif search_on=="seqid":
            where = f"feature.seqfeature_id = {self.placeholder} "
            join = "INNER JOIN seqfeature AS feature ON feature.bioentry_id=t1.bioentry_id "
        else:
            raise RuntimeError(f"Unsupported argument: {search_on}, must be either taxid or seqid")

        length_term = ""
        if not min_bioentry_length is None:
            length_term = f" AND length > {min_bioentry_length}"
        
        query = (
            "SELECT t2.bioentry_id, t1.accession, length, seq " 
            "FROM bioentry t1 "
            "INNER JOIN biosequence t2 ON t1.bioentry_id=t2.bioentry_id "
            f"{join}"
            f"WHERE {where} {length_term};"
        )
        
        results = self.server.adaptor.execute_and_fetchall(query, [search_term])
        header = ["bioentry_id", "accession" ,"length", "seq"]

        if search_on=="seqid":
            # assumes only one contig for a given seqid
            # bioentry_id, accession, length, seq
            return results[0]
        return DB.to_pandas_frame(results, header).set_index(["bioentry_id"])

    
    # return table of all features of a target genome with location,  (filter by term name) 
    # basically the same as get_proteins_info but with rrna and trna and location
    # possibility to merge with get_proteins_info? # locus, prot_id, gene, product
    def get_features_location(self, 
                              search_id, search_on="taxon_id",
                              seq_term_names=['CDS', 'rRNA', 'tRNA'],
                              seq_qual_term_names=['locus_tag', 'product', 'gene']):

        if search_on=="taxon_id":
            where_clause = f"t1.taxon_id = {self.placeholder}"
        elif search_on=="bioentry_id":
            where_clause = f"t1.bioentry_id = {self.placeholder}"
        else: 
            raise RuntimeError("Expected either taxon_id or bioentry_id")
        
        seq_term_query = ",".join(f"\"{name}\"" for name in seq_term_names)
        qual_term_query = ",".join(f"\"{name}\"" for name in seq_qual_term_names)
        
        query = (
            "select t1.bioentry_id, t2.seqfeature_id,t5.start_pos, t5.end_pos, t5.strand, t4.name as term_name,t6.value,t7.name as qualifier  from bioentry t1 "
            "inner join seqfeature t2 on t1.bioentry_id = t2.bioentry_id "
            "inner join term t4 on t2.type_term_id = t4.term_id "
            "inner join location t5 on t2.seqfeature_id = t5.seqfeature_id "
            "inner join seqfeature_qualifier_value t6 on t2.seqfeature_id =t6.seqfeature_id "
            "inner join term t7 on t6.term_id = t7.term_id "
            f"where {where_clause} and t4.name in ({seq_term_query}) and t7.name in ({qual_term_query}); "
        )
        
        results = self.server.adaptor.execute_and_fetchall(query, [search_id])
        
        df = DB.to_pandas_frame(results, ["bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand", "term_name", "qualifier_value", "qualifier_name"])
        df = df.set_index(["bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand", "term_name", "qualifier_name"]).unstack("qualifier_name")
        df = df.reset_index()
        df.columns = ['_'.join(col).strip('_') for col in df.columns]

        # bioentry_id  start_pos  end_pos  strand  gene locus_tag product
        df_merged = df[["seqfeature_id", "bioentry_id", "start_pos", "end_pos", "strand", "qualifier_value_gene", "term_name", "qualifier_value_locus_tag", "qualifier_value_product"]]
        
        return df_merged
        

    def get_identity_closest_homolog(self, reference_taxid, target_taxids):
        
        targets = ','.join([str(i) for i in target_taxids])
        query = (
        "select t1.id_1, t1.id_2, t1.\"identity\", t5.taxon_id from orthology_identity t1 "
        "inner join seqfeature t2 on t1.id_1=t2.seqfeature_id " 
        "inner join seqfeature t3 on t1.id_2=t3.seqfeature_id " 
        "inner join bioentry t4 on t2.bioentry_id = t4.bioentry_id " 
        "inner join bioentry t5 on t3.bioentry_id = t5.bioentry_id "
        f"where t4.taxon_id = {self.placeholder} and t5.taxon_id in ({targets}) and t4.taxon_id != t5.taxon_id " 
        "UNION "
        "select t1.id_2, t1.id_1, t1.\"identity\", t4.taxon_id from orthology_identity t1 "
        "inner join seqfeature t2 on t1.id_1=t2.seqfeature_id "
        "inner join seqfeature t3 on t1.id_2=t3.seqfeature_id "
        "inner join bioentry t4 on t2.bioentry_id = t4.bioentry_id " 
        "inner join bioentry t5 on t3.bioentry_id = t5.bioentry_id "
        f"where t4.taxon_id in ({targets}) and t5.taxon_id = {self.placeholder} and t4.taxon_id != t5.taxon_id " 
        )
        
        results = self.server.adaptor.execute_and_fetchall(query,[reference_taxid, reference_taxid])
        df = DB.to_pandas_frame(results, ["seqfeature_id_1", "seqfeature_id_2", "identity", "target_taxid"])
        # keep homolog with the highest identity
        #df_pivot = df.pivot_table(index=["seqfeature_id_1"], columns="target_taxid",values="identity", aggfunc=lambda x: max(x))
        return df
    

    def get_pfam_def(self, pfam_ids, add_ttl_count=False):
        ttl_cnt, ttl_join, ttl_grp = "", "", ""
        ttl_join_template = (
            " INNER JOIN pfam_hits AS hits ON hits.pfam_id=pfam_defs.pfam_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON hsh.hsh=hits.hsh "
        )
        ttl_grp_template = "GROUP BY pfam_defs.pfam_id"

        if pfam_ids is None:
            where = ""
            ttl_join = ttl_join_template
            ttl_grp = ttl_grp_template
        else:
            plcd = self.gen_placeholder_string(pfam_ids)
            where = f"WHERE pfam_defs.pfam_id IN ({plcd}) "

        if add_ttl_count:
            ttl_cnt = ", COUNT(*) "
            ttl_grp = ttl_grp_template
            ttl_join = ttl_join_template

        query = (
            f"SELECT pfam_defs.pfam_id, pfam_defs.definition {ttl_cnt} "
            f"FROM pfam_table AS pfam_defs {ttl_join} "
            f"{where}"
            f"{ttl_grp};"
        )
        results = self.server.adaptor.execute_and_fetchall(query, pfam_ids)
        cols = ["pfam", "def"]
        if add_ttl_count:
            cols.append("ttl_cnt")
        return DB.to_pandas_frame(results, cols).set_index(["pfam"])


    def gen_pfam_where_clause(self, search_on, entries):
        entries = self.gen_placeholder_string(entries)
        if search_on=="seqid":
            where_clause = f" hsh.seqid IN ({entries}) "
        elif search_on=="pfam":
            where_clause = f" pfam.pfam_id IN ({entries}) "
        elif search_on=="taxid":
            where_clause = f" entry.taxon_id IN ({entries}) "
        else:
            raise RuntimeError(f"Searching on {search_on} is not supported")
        return where_clause


    def get_pfam_hits_info(self, seqids):
        plcd = self.gen_placeholder_string(seqids)
        query = (
            "SELECT hsh.seqid, pfam.pfam_id, pfam.start, pfam.end "
            "FROM sequence_hash_dictionnary AS hsh "
            "INNER JOIN pfam_hits AS pfam ON hsh.hsh=pfam.hsh "
            f"WHERE hsh.seqid IN ({plcd});"
        )
        results = self.server.adaptor.execute_and_fetchall(query, seqids)
        return DB.to_pandas_frame(results, ["seqid", "pfam", "start", "end"])


    def get_pfam_hits(self, ids, indexing="taxid", search_on="taxid", plasmids=None, keep_taxid=False):

        where_clause = self.gen_pfam_where_clause(search_on, ids)
        header = []
        if indexing=="seqid":
            index = "seqid.seqfeature_id"
            header.append("seqid")
            if keep_taxid:
                index += ", entry.taxon_id "
                header.append("taxid")
        elif indexing=="taxid":
            header.append("taxid")
            index = "entry.taxon_id"
        else:
            raise RuntimeError(f"Indexing method not supported: {indexing}")

        plasmid_join = ""
        if not plasmids is None:
            index += ", CAST(is_plasmid.value AS int) "
            subclause = self.gen_pfam_where_clause(search_on, plasmids)
            where_clause = (
                    f"({where_clause} AND is_plasmid.value=0) "
                    f" OR ({subclause} AND is_plasmid.value=1)"
            )
            plasmid_join = (
                "INNER JOIN bioentry_qualifier_value AS is_plasmid ON "
                "  is_plasmid.bioentry_id=entry.bioentry_id "
                "INNER JOIN term AS plasmid_term ON plasmid_term.term_id=is_plasmid.term_id "
                "  AND plasmid_term.name=\"plasmid\""
            )
            header.append("plasmid")

        header.extend(["pfam", "count"])
        query = (
            f"SELECT {index}, pfam.pfam_id, COUNT(*) "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id = entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON seqid.seqfeature_id = hsh.seqid "
            "INNER JOIN pfam_hits AS pfam ON pfam.hsh = hsh.hsh "
            f"{plasmid_join}"
            f"WHERE {where_clause} "
            f"GROUP BY {index}, pfam.pfam_id;"
        )
        all_ids = ids
        if not plasmids is None:
            all_ids += plasmids
        results = self.server.adaptor.execute_and_fetchall(query, all_ids)

        if indexing=="taxid":
            df_index = [indexing, "pfam"]
            if not plasmids is None:
                df_index.insert(1, "plasmid")
            df = DB.to_pandas_frame(results, header)
            df = df.set_index(df_index).unstack(level=0, fill_value=0)

            if not plasmids is None:
                return df.unstack(level=0, fill_value=0)
            elif not df.empty:
                df.columns = [col for col in df["count"].columns.values]
            else:
                return df
        elif indexing=="seqid":
            if not plasmids is None:
                raise RuntimeError("Not implemented for now")
            header = ["seqid", "pfam"]
            if keep_taxid:
                header.append("taxid")
                results = ((seqid, pfam, taxid) for seqid, taxid, pfam, count in results)
            else:
                results = ((seqid, pfam) for seqid, pfam, count in results)
            df = DB.to_pandas_frame(results, header)
            # NOTE: do not index on seqid a protein can have several domains!
        return df


    def gen_cog_where_clause(self, search_on, entries):
            entries = self.gen_placeholder_string(entries)

            if search_on=="bioentry":
                where_clause = f" entry.bioentry_id IN ({entries}) "
            elif search_on=="seqid":
                where_clause = f" hsh.seqid IN ({entries}) "
            elif search_on=="cog":
                where_clause = f" cogs.cog_id IN ({entries}) "
            elif search_on=="taxid":
                where_clause = f" entry.taxon_id IN ({entries}) "
            else:
                raise RuntimeError(f"Searching on {search_on} is not supported")
            return where_clause


    # Get all cog hits for a given list of bioentries
    # The results are either indexed by the bioentry or by the seqid
    # NOTE: if indexing as bioentry or taxon_id, will return a dataframe of the form
    #           bioentry_1/taxon_1  bioentry_2/taxon_2  ...
    #   cog_1     cnt                   cnt
    #   cog_2     cnt                   cnt
    #   ...
    #
    # If the indexing is seqid, will return a dataframe with the format
    #   seqid1 cog1
    #   seqid2 cog2
    #   seqid3 cog3
    def get_cog_hits(self, ids, indexing="bioentry", search_on="bioentry", keep_taxid=False, plasmids=None):

        where_clause = self.gen_cog_where_clause(search_on, ids)
        if indexing=="seqid":
            index = "seqid.seqfeature_id"
            if keep_taxid:
                index += ", entry.taxon_id "
        elif indexing=="bioentry":
            index = "entry.bioentry_id"
        elif indexing=="taxid":
            index = "entry.taxon_id"
        else:
            raise RuntimeError(f"Indexing method not supported: {indexing}")

        plasmid_join = ""
        if not plasmids is None:
            index += ", CAST(is_plasmid.value AS int) "
            subclause = self.gen_cog_where_clause(search_on, plasmids)
            where_clause = (
                    f"({where_clause} AND is_plasmid.value=0) "
                    f" OR ({subclause} AND is_plasmid.value=1)"
            )
            plasmid_join = (
                "INNER JOIN bioentry_qualifier_value AS is_plasmid ON "
                "  is_plasmid.bioentry_id=entry.bioentry_id "
                "INNER JOIN term AS plasmid_term ON plasmid_term.term_id=is_plasmid.term_id "
                "  AND plasmid_term.name=\"plasmid\""
            )

        query = (
            f"SELECT {index}, cogs.cog_id, COUNT(*) "
            "FROM bioentry AS entry "
            "INNER JOIN seqfeature AS seqid ON seqid.bioentry_id = entry.bioentry_id "
            "INNER JOIN sequence_hash_dictionnary AS hsh ON seqid.seqfeature_id = hsh.seqid "
            "INNER JOIN cog_hits AS cogs ON cogs.hsh = hsh.hsh "
            f"{plasmid_join}"
            f"WHERE {where_clause} "
            f"GROUP BY {index}, cogs.cog_id;"
        )
        all_ids = ids
        if not plasmids is None:
            all_ids += plasmids
        results = self.server.adaptor.execute_and_fetchall(query, all_ids)

        # ugly code, open for improvements
        if indexing=="taxid" or indexing=="bioentry":
            column_names = [indexing, "cog", "count"]
            index = [indexing, "cog"]
            if not plasmids is None:
                column_names.insert(1, "plasmid")
                index.insert(1, "plasmid")
            df = DB.to_pandas_frame(results, column_names)
            df = df.set_index(index).unstack(level=0, fill_value=0)

            if not plasmids is None:
                return df.unstack(level=0, fill_value=0)
            else:
                df.columns = [col for col in df["count"].columns.values]
        elif indexing=="seqid":
            if not plasmids is None:
                raise RuntimeError("Not implemented for now")
            header = ["seqid", "cog"]
            if keep_taxid:
                header.append("taxid")
                results = ((seqid, cog, taxid) for seqid, taxid, cog, count in results)
            else:
                results = ((seqid, cog) for seqid, cog, count in results)

            df = DB.to_pandas_frame(results, header)
            df = df.set_index(["seqid"])
        return df


    def get_filenames_to_taxon_id(self):
        sql = (
            f"SELECT filename, taxon_id FROM filenames;"
        )
        hsh_filenames_to_entry = {}
        results = self.server.adaptor.execute_and_fetchall(sql)
        return {filename: entry_id for filename, entry_id in results}


    def count_orthogroups(self):
        sql = ("SELECT COUNT(*)"
                "FROM filenames;"
        )
        number = str(self.server.adaptor.execute_and_fetchall(sql))
        return number [2]


    def get_taxon_id_to_filenames(self):
        sql = (
            f"SELECT filename, taxon_id FROM filenames;"
        )
        hsh_entry_to_filenames = {}
        results = self.server.adaptor.execute_and_fetchall(sql)
        return {entry_id: filename for filename, entry_id in results}


    def get_bioentry_id_to_plasmid_contigs(self):
        sql = (
            "SELECT bioentry_id, accession "
            "FROM bioentry "
            "WHERE bioentry_id IN (SELECT bioentry_qualifier_value.bioentry_id "
            "FROM bioentry_qualifier_value "
            "WHERE bioentry_qualifier_value.term_id=47 AND bioentry_qualifier_value.value=1) "
        )
        bioentry_id_to_plasmid_contigs = {}
        results = self.server.adaptor.execute_and_fetchall(sql)
        return {bioentry_id: accession for bioentry_id, accession in results}
    

    def load_chlamdb_config_tables(self, entries):
        sql = (
            "CREATE TABLE biodb_config"
            "(name varchar(200), type varchar(200), status BOOLEAN);"
        )
        self.server.adaptor.execute(sql)
        self.load_data_into_table("biodb_config", entries)


    def update_taxon_ids(self, update_lst):
        query = (
            "UPDATE bioentry SET taxon_id=? WHERE bioentry_id=?;"
        )
        for bioentry_id, taxon_id in update_lst:
            self.server.adaptor.execute(query, (taxon_id, bioentry_id))


    def update_plasmid_status(self, plasmid_bioentries):
        plasmid_term_id = self.get_term_id("plasmid", create_if_absent=True)
        sql = (
            "INSERT INTO bioentry_qualifier_value VALUES  "
            " (?, ?, ?, 0);"
        )
        data = [(bioentry_id, plasmid_term_id, is_plasmid)
                for bioentry_id, is_plasmid in plasmid_bioentries]
        self.server.adaptor.executemany(sql, data)

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

    
    def gen_placeholder_string(self, args):
        return ",".join(self.placeholder for _ in args)


    # wrapper methods
    def commit(self):
        self.server.commit()


    def load_gbk_wrapper(self, records):
        self.server[self.db_name].load(records)


    # Maybe return different instance of a subclass depending on the type
    # of database? Would allow to avoid code duplication if several database
    # types are to be included.
    def load_db(db_file, params):
        sqlpsw = params.get("zdb.psswd", "")
        db_type = params["zdb.db_type"]
        db_name = params["zdb.db_name"]

        if db_type != "sqlite":
            server = BioSeqDatabase.open_database(driver="MySQLdb", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1", 
                                                  db=db_file, 
                                                  charset='utf8',
                                                  use_unicode=True)
        else:
            server = BioSeqDatabase.open_database(driver="sqlite3", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1",
                                                  db=db_file)
        return DB(server, db_name)


    def load_db_from_name(db_name, db_type = "sqlite"):
        params = {"zdb.db_type" : db_type, "zdb.db_name" : db_name}
        return DB.load_db(db_name, params)


    def location2sequence(self, accession, start, end):
        sql = (
            "select substr(seq, %s, %s) from biosequence "
            "inner join bioentry on bioentry.bioentry_id=biosequence.bioentry_id "
            "inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id "
            "where accession='%s' " % (start, end, accession)
            )
        
        #sequence = str(self.server.adaptor.execute_and_fetchall(sql))
        #return sequence 
        sequence = str(self.server.adaptor.execute_and_fetchall(sql))
        return sequence [0]

