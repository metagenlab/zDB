#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-

from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
import sys
import mysqldb_load_mcl_output

# A faire
# add arthogroup to db_xref, chose an orthogroup prefix (chlam_x)
# get_orthogroup(orthogroup_id), return list of protein ids
# get_genome_id_from_locus_tag


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


def load_db(db_name=False):
    
    server = BioSeqDatabase.open_database(driver="MySQLdb", user="tpillone",
                       passwd = "agnathe3", host = "localhost", db="biosqldb")
    if db_name:
        try:
            db = server[db_name]
            return server, db
        except:
            if query_yes_no("Database %s do not exist, create it?" % db_name):
                print "creating new db...", db_name
                db = server.new_database(db_name, description="db %s" % db_name)
                return server, db
            else:
                sys.exit("Stopping execution")
    else:
        return server


def _to_dict(tuple_list):
    temp_dict = {}
    for i in tuple_list:
        if len(i[1:]) >1:
            temp_dict[i[0]] = i[1:]
        else:
            temp_dict[i[0]] = i[1]
    return temp_dict


def get_orthogroup(orthogroup_id, biodatabase_name):
    sql = ""



# get all seqfeature qualifiers from seqfeature id
# SELECT name, value FROM seqfeature_qualifier_value join term using (term_id) WHERE seqfeature_id = 5498;
def seqfeature_id2seqfeature_object_dict(*DBSeqRecord_objects):
    seqfeature_id2seqfeature_object = {}
    for DBSeqRecord_object in DBSeqRecord_objects:
        temp_dict = {}
        for seqfeature in DBSeqRecord_object.features:
            temp_dict[seqfeature._seqfeature_id] = seqfeature
        seqfeature_id2seqfeature_object.update(temp_dict)
    return seqfeature_id2seqfeature_object


# get seqfeature_if from protein_name
# Attention: should filter by biodatabase
#SELECT seqfeature_id FROM seqfeature_qualifier_value join term using (term_id) WHERE value = "CAQ51124.1";
def get_bioentry_and_seqfeature_id_from_locus_tag_list(server, locus_tag_list, biodatabase_name):


    query_locus_tag = '('
    for i in range(0, len(locus_tag_list)-1):
        query_locus_tag+='value = "%s" or ' % str(locus_tag_list[i])
    query_locus_tag+='value = "%s")' % locus_tag_list[-1]
    
    sql = 'select t6.accession, t1.seqfeature_id from biosqldb.seqfeature as t1' \
          ' inner join seqfeature_qualifier_value as t2 on t1.seqfeature_id = t2.seqfeature_id and %s' \
          ' inner join term as t3 on t2.term_id = t3.term_id and t3.name = "locus_tag"' \
          ' inner join term as t4 on t1.type_term_id = t4.term_id and t4.name = "CDS"' \
          ' inner join term as t5 on t1.source_term_id = t5.term_id' \
          ' inner join bioentry as t6 on t1.bioentry_id = t6.bioentry_id' \
          ' inner join biodatabase on t6.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"'

    
    sql = sql % (query_locus_tag, biodatabase_name)
    print sql
    b = server.adaptor.execute_and_fetchall(sql, )
    bioentry_id_list = [i[0] for i in b]
    seqfeature_id_list = [i[1] for i in b]
    return (bioentry_id_list, seqfeature_id_list)


def get_bioentry_id_from_locus_tag(server, locus_tag, biodatabase_name):
    
    sql_specific_locus_tag = "select bioentry.name from seqfeature_qualifier_value" \
          " inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = %s and value = %s" \
          " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
          " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
          " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s"

    b = server.adaptor.execute_and_fetchall(sql_specific_locus_tag, ("locus_tag", locus_tag, biodatabase_name))
    return b[0][0]


def taxon_id2genome_description(server, biodatabase_name):
    print "bonjour"
    sql = 'select taxon_id, bioentry.description from bioentry ' \
          ' inner join biodatabase on biodatabase.biodatabase_id = bioentry.biodatabase_id' \
          ' where biodatabase.name = "%s" and bioentry.description not like "%%%%plasmid%%%%" ' % biodatabase_name
    print sql
    result = server.adaptor.execute_and_fetchall(sql, )
    return _to_dict(result)


def accession2description(server, biodatabase_name):
    sql = 'select bioentry.accession, bioentry.description from bioentry' \
          ' inner join biodatabase on biodatabase.biodatabase_id = bioentry.biodatabase_id' \
          ' where biodatabase.name = "%s"' % biodatabase_name
    result = server.adaptor.execute_and_fetchall(sql, )
    return _to_dict(result)


def locus_tag2bioentry_id_dict(server, biodatabase_name):
    
    sql_locus_tag_bioentry_id_table = "select value, bioentry.name from seqfeature_qualifier_value" \
          " inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = %s" \
          " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
          " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
          " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s"

    result = server.adaptor.execute_and_fetchall(sql_locus_tag_bioentry_id_table, ("locus_tag", biodatabase_name))

    
    
    locus_tag2bioentry_id = {}
    for protein in result:
       locus_tag2bioentry_id[protein[0]] = protein[1] 
    return locus_tag2bioentry_id


def locus_tag2seqfeature_id_dict(server, biodatabase_name):



    sql_locus_tag_seqfeature_id_table = 'select t2.value, t1.seqfeature_id from biosqldb.seqfeature as t1' \
                                        ' inner join seqfeature_qualifier_value as t2 on t1.seqfeature_id = t2.seqfeature_id' \
                                        ' inner join term as t3 on t2.term_id = t3.term_id and t3.name = "locus_tag"' \
                                        ' inner join term as t4 on t1.type_term_id = t4.term_id and t4.name = "CDS"' \
                                        ' inner join term as t5 on t1.source_term_id = t5.term_id' 


    result = server.adaptor.execute_and_fetchall(sql_locus_tag_seqfeature_id_table, )

    locus_tag2seqfeature_id = {}
    for locus_tag in result:
        locus_tag2seqfeature_id[locus_tag[0]] = locus_tag[1]
    return locus_tag2seqfeature_id


def protein_id2seqfeature_id_dict(server, biodatabase_name):

    sql_protein_id2_seqfeature_id_table = "select value, seqfeature.seqfeature_id from seqfeature_qualifier_value" \
    " inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = %s" \
    " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
    " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s"
    
    result = server.adaptor.execute_and_fetchall(sql_protein_id2_seqfeature_id_table, ("protein_id", biodatabase_name))
    protein_id2seqfeature_id = {}
    for protein in result:
       print "protein", protein
       #seqfeature_id = locus_tag2CDS_seqfeature_id(server, protein[0], biodatabase_name) 
       protein_id2seqfeature_id[protein[0]] = protein[1]
    return protein_id2seqfeature_id

#def locus_tag2seqfeature_id(server, locus_tag, biodatabase_name):
#    sql_locus_tag_seqfeature_id_table = "select value, seqfeature.seqfeature_id from seqfeature_qualifier_value" \
#    " inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = %s and value = %s" \
#    " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
#    " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
#    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s"
    
#    result = server.adaptor.execute_and_fetchall(sql_locus_tag_seqfeature_id_table, ("locus_tag", locus_tag, biodatabase_name))
#    return result[0][0]


def bioentry_id2genome_name_dict(server, biodatabase_name):

    sql_bioentry_id2genome_name = "select bioentry.name, bioentry.description from bioentry" \
                                     " join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name =%s;"
    
    result = server.adaptor.execute_and_fetchall(sql_bioentry_id2genome_name, (biodatabase_name))
    bioentry_id2genome_name = {}
    for protein in result:
       bioentry_id2genome_name[protein[0]] = protein[1] 
    return bioentry_id2genome_name


def locus_tag2bioentry_accession_and_description(server, locus_tag, biodatabase_name):
    
    sql ="select bioentry.accession, bioentry.description from seqfeature_qualifier_value" \
    " inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = %s and value = %s" \
    " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
    " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s" \
    
    result = server.adaptor.execute_and_fetchall(sql, ("locus_tag", locus_tag, biodatabase_name))

    return [i for i in result[0]]


def locus_tag2orthogroup_id(server, locus_tag, biodatabase_name):
    sql_locus_tag2 = "select value, seqfeature.seqfeature_id from seqfeature_qualifier_value" \
    " inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = %s" \
    " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
    " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s"
    
    result = server.adaptor.execute_and_fetchall(sql_locus_tag_seqfeature_id_table, ("locus_tag", biodatabase_name))


def seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodatabase_name):
    
    sql ="select term.name,  seqfeature_qualifier_value.value  from seqfeature_qualifier_value" \
    " inner join term on seqfeature_qualifier_value.term_id = term.term_id and seqfeature_id = %s" \
    " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
    " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s" 
    
    result = server.adaptor.execute_and_fetchall(sql, (seqfeature_id, biodatabase_name))
    return _to_dict(result)


def locus_tag2seqfeature_qualifier_values(server, locus_tag, biodatabase_name):

    seqfeature_id = locus_tag2CDS_seqfeature_id(server, locus_tag, biodatabase_name)
    
    return seqfeature_id2seqfeature_qualifier_values(server, seqfeature_id, biodatabase_name)


def orthogroup_id2seqfeature_id_list(server, orthogroup_id, biodatabase_name):
    sql ="select seqfeature.seqfeature_id   from seqfeature_qualifier_value" \
    " inner join term on seqfeature_qualifier_value.term_id = term.term_id and value = %s" \
    " inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id" \
    " inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s"

    result = server.adaptor.execute_and_fetchall(sql, (orthogroup_id, biodatabase_name))
    return [i[0] for i in result]


def get_biodatabase_list(server):
    sql = "select name from biodatabase;"
    return server.adaptor.execute_and_fetchall(sql,)


def get_genome_description_list(server, biodatabase_name):
    """
    return the list of genome description for a given database
    """

    sql = "select bioentry.name from bioentry" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s" 
    result = server.adaptor.execute_and_fetchall(sql, (biodatabase_name))
    return [i[0] for i in result]


def get_taxon_id_list(server, biodatabase_name):


    sql = "select taxon_id from bioentry" \
    " inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id" \
    " where biodatabase.name = %s" \
    " group by taxon_id"
    result = server.adaptor.execute_and_fetchall(sql, (biodatabase_name))
    return [i[0] for i in result]
    
    
def orthogroup_id2locus_tag_list(server, orthogroup_id, biodatabase_name):
    
    seqfeaure_id_list =  orthogroup_id2seqfeature_id_list(server, orthogroup_id, biodatabase_name)

    query_seqfeature_id = "("
    for i in range(0, len(seqfeaure_id_list)-1):
        query_seqfeature_id+="seqfeature_id = %s or " % str(seqfeaure_id_list[i])
    query_seqfeature_id+="seqfeature_id = %s)" % seqfeaure_id_list[-1]
    #print query_seqfeature_id


    sql_seqfeature_id2locus_tag = 'select bioentry.accession, seqfeature_qualifier_value.seqfeature_id, seqfeature_qualifier_value.value, bioentry.description from seqfeature_qualifier_value' \
    ' inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "locus_tag" and %s' \
    ' inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id' \
    ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
    ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' \
    ' group by seqfeature_qualifier_value.seqfeature_id' % (query_seqfeature_id, biodatabase_name)

    #print sql_seqfeature_id2locus_tag
    
    result = server.adaptor.execute_and_fetchall(sql_seqfeature_id2locus_tag,)
    return result


def protein_id2seqfeature_id(server, protein_id, biodatabase_name):
    
    sql = 'select seqfeature_qualifier_value.seqfeature_id from seqfeature_qualifier_value' \
          ' inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id and value = "%s"' \
          ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
          ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' % (protein_id, biodatabase_name)
    result = server.adaptor.execute_and_fetchall(sql,)
    return result[0][0]


def protein_id2locus_tag(server, protein_id, biodatabase_name):

    seqfeature_id = protein_id2seqfeature_id(server, protein_id, biodatabase_name)

    sql = 'select seqfeature_qualifier_value.value from seqfeature_qualifier_value' \
          ' inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "locus_tag"' \
          ' inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id and seqfeature_qualifier_value.seqfeature_id = %s' \
          ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
          ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' % (seqfeature_id, biodatabase_name)
    result = server.adaptor.execute_and_fetchall(sql,)
    return result[0][0]  

   
def locus_tag2seqfeature_ids(server, locus_tag, biodatabase_name): 
    sql_locus_tag_seqfeature_id = 'select seqfeature_qualifier_value.seqfeature_id from seqfeature_qualifier_value' \
    ' inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "locus_tag" and value = %s' \
    ' inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id' \
    ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
    ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s'
    result = server.adaptor.execute_and_fetchall(sql_locus_tag_seqfeature_id, (locus_tag, biodatabase_name))
    return [i[0] for i in result]


def locus_tag2CDS_seqfeature_id(server, locus_tag, biodatabase_name):
    seqfeature_ids = locus_tag2seqfeature_ids(server, locus_tag, biodatabase_name)

    sql_seqids = "seqfeature.seqfeature_id = %s" % seqfeature_ids[0]

    if len(seqfeature_ids) >1:
        for i in range(1,len(seqfeature_ids)):
            sql_seqids+=" or seqfeature.seqfeature_id = %s" % seqfeature_ids[i]

    sql = 'select seqfeature_id  from seqfeature' \
          ' inner join term on seqfeature.type_term_id = term.term_id' \
          ' where name = "CDS" and (%s)' % (sql_seqids)
        
    result = server.adaptor.execute_and_fetchall(sql, )
    return result[0][0]

def locus_or_protein_id2taxon_id(server, db_name):
    sql = 'select value, taxon_id, term.name from seqfeature_qualifier_value' \
          'inner join term on term.term_id = seqfeature_qualifier_value.term_id' \
          ' inner join seqfeature as t2 on t2.seqfeature_id = seqfeature_qualifier_value.seqfeature_id' \
          ' inner join bioentry on t2.bioentry_id = bioentry.bioentry_id' \
          ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "%s"' \
          ' where term.name = "locus_tag" or term.name = "protein_id"'' % db_name

    result = server.adaptor.execute_and_fetchall(sql, )
    return _to_dict(result)




def bioentry_name2orthoup_size(server, biodatabase_name, bioentry_name):

    # get number of proteins for each orthogroup of taxon_id 4
    sql = 'select value, count(value) from seqfeature_qualifier_value as ortho_table' \
    ' inner join term on ortho_table.term_id = term.term_id and name = "orthogroup"' \
    ' inner join seqfeature on ortho_table.seqfeature_id = seqfeature.seqfeature_id' \
    ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
    ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s' \
    ' where bioentry.name = %s' \
    ' group by value'

    result = server.adaptor.execute_and_fetchall(sql, (biodatabase_name, bioentry_name))
    result = _to_dict(result)
    return result


def get_column_names(server, table):
  sql= 'show columns from %s' % table
  result = server.adaptor.execute_and_fetchall(sql, )
  return [i[0] for i in result]
  

def taxon_id2orthogroup_size(server, biodatabase_name, taxon_id):

    # get number of proteins for each orthogroup of taxon_id 4
    sql = 'select value, count(value) from seqfeature_qualifier_value as ortho_table' \
    ' inner join term on ortho_table.term_id = term.term_id and name = "orthogroup"' \
    ' inner join seqfeature on ortho_table.seqfeature_id = seqfeature.seqfeature_id' \
    ' inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id' \
    ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = %s' \
    ' where bioentry.taxon_id = %s' \
    ' group by value'

    result = server.adaptor.execute_and_fetchall(sql, (biodatabase_name, taxon_id))
    result = _to_dict(result)
    return result


def seqfeature_id2feature_location(server, seqfeature_id):
    sql ='select start_pos, end_pos, strand from location where seqfeature_id= %s and rank = 1'
    result = server.adaptor.execute_and_fetchall(sql, seqfeature_id)
    return result[0]
    

def get_orthology_table(server, biodb_name):
    sql ='select * from orthology_%s' % biodb_name
    result = server.adaptor.execute_and_fetchall(sql,)
    return result


def get_orthology_table_subset(server, biodb_name, reference_taxon_id):
    sql = 'select * from orthology_%s where `%s` > 0' % reference_taxon_id
    result = server.adaptor.execute_and_fetchall(sql,)
    return result

    
if __name__ == '__main__':
    pass
    #server, db = load_db("Chlamydiales_subgroup")

    #print locus_tag2bioentry_id_dict(server, "test")
    #locus_tag2seqfeature_id = locus_tag2seqfeature_id_dict(server, "chlamydiales")
    #print add_orthogroup_term(server)
    #print bioentry_id2genome_name_dict(server, "test")
    #del server["test"]
    #server.commit()
    #print locus_tag2seqfeature_id(server, "CAQ51124.1", "chlamydiales")
    #print orthogroup_id2seqfeature_id_list(server, "group_1", "chlamydiales")
    #print orthogroup_id2locus_tag_list(server, "group_1", "chlamydiales")
    #print locus_tag2bioentry_accession_and_description(server, "chlamydiales", "CAQ51124.1")
    #locus_tag2orthogroup_id, orthomcl_groups2proteins, genome_orthomcl_code2proteins, locus_tag2genome_ortho_mcl_code = parse_orthomcl_output("chlam_subgroup/merged_mcl.txt")
    #locus_tag2seqfeature_id = locus_tag2seqfeature_id_dict(server, "Chlamydiales_subgroup")
    #add_orthogroup_to_seq(server, locus_tag2orthogroup_id, locus_tag2seqfeature_id)

    

# seqfeature
#| seqfeature_id | bioentry_id | type_term_id | source_term_id | display_name | rank |

# seqfeature_qualifier_value

# | term_id | seqfeature_id | rank | value      | name       | definition | identifier | is_obsolete | ontology_id |




"""
template
SELECT ...
FROM a
LEFT JOIN b on a.somefield=b.somefield
LEFT JOIN c on b.otherfield=c.otherfield
"""

# get accession, no bioentry filtering
#select accession from bioentry
#LEFT JOIN seqfeature on bioentry.bioentry_id = seqfeature.bioentry_id
#LEFT JOIN seqfeature_qualifier_value on seqfeature.seqfeature_id = seqfeature_qualifier_value.seqfeature_id
#where seqfeature_qualifier_value.value = "CAQ51124.1";


# get accession, no bioentry filtering
#select accession from bioentry
#LEFT JOIN biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id where biodatabase.name="tata"
#LEFT JOIN seqfeature on bioentry.bioentry_id = seqfeature.bioentry_id
#LEFT JOIN seqfeature_qualifier_value on seqfeature.seqfeature_id = seqfeature_qualifier_value.seqfeature_id
#where seqfeature_qualifier_value.value = "CAQ51124.1";

# obtenir table locus_tag | sequfeature_id
#select name,value,seqfeature_id from seqfeature_qualifier_value
#join term on seqfeature_qualifier_value.term_id = term.term_id
#where name = "locus_tag";


# obtenir table locus_tag | accession | biodatabase
#select value, bioentry_id  from seqfeature_qualifier_value, seqfeature
#inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "locus_tag";



# table avec protein id | bioentry id
#select name, value, seqfeature.bioentry_id from seqfeature_qualifier_value
#inner join term on seqfeature_qualifier_value.term_id = term.term_id
#inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
#where name = "locus_tag";

# table avec protein id | bioentry name
#select value, bioentry.name, biodatabase.name  from seqfeature_qualifier_value
#inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "locus_tag" and value = "CAQ51124.1"
#inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
#inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id
#inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "tata"


"""
# return group_1 accession and description
select * from seqfeature_qualifier_value
inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "uniprot_test"
"""

"""
select * from bioentry_qualifier_value
inner join bioentry on bioentry_qualifier_value.bioentry_id = bioentry.bioentry_id
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "uniprot_test"
inner join term on term.term_id = bioentry_qualifier_value.term_id
where accession = "P0A5B8";
"""

"""
# return group size
select COUNT(*) from seqfeature_qualifier_value
inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "orthogroup" and value = "group_222"
inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "Chlamydiales_subgroup_no_estrella"
"""

"""

"""


"""
# retourne le nombre de contigs par génome
select COUNT(*), description from bioentry
GROUP BY taxon_id;
"""

"""
# retourne le nombre de contigs par génome
select * from bioentry where name = "FR872592"
GROUP BY description;
"""




"""
# return group_1 accession and description 
select bioentry.accession, bioentry.description, seqfeature_qualifier_value.value, bioentry.taxon_id  from seqfeature_qualifier_value
inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "orthogroup" and value = "group_1777"
inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id 
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "chlamydiales"
"""





"""
# return group_1 accession and description
select * from seqfeature_qualifier_value as ortho_table 
inner join term on ortho_table.term_id = term.term_id and name = "orthogroup" 
inner join seqfeature on ortho_table.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id 
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "Chlamydiales_subgroup_no_estrella"
where taxon_id = 64
group by ortho_table.value
"""


"""
# get group size
select count(*) from seqfeature_qualifier_value as ortho_table 
inner join term on ortho_table.term_id = term.term_id and name = "orthogroup" and value = "group_76"
inner join seqfeature on ortho_table.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id 
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "chlamydiales"
"""


"""
# get number of proteins for each orthogroup of taxon_id 4
select value, count(value) from seqfeature_qualifier_value as ortho_table 
inner join term on ortho_table.term_id = term.term_id and name = "orthogroup" 
inner join seqfeature on ortho_table.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id 
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "Chlamydiales_subgroup_no_estrella"
where bioentry.name = "WchW"
group by value
"""

"""
# get number of proteins for each orthogroup => pour faire distrib taille des groupes
select value, count(value) from seqfeature_qualifier_value as ortho_table 
inner join term on ortho_table.term_id = term.term_id and name = "orthogroup" 
inner join seqfeature on ortho_table.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id 
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "Chlamydiales_subgroup_no_estrella"
group by value
"""


"""
select bioentry.accession, seqfeature_qualifier_value.seqfeature_id, seqfeature_qualifier_value.value from seqfeature_qualifier_value
inner join term on seqfeature_qualifier_value.term_id = term.term_id and name = "locus_tag" and (seqfeature_id = 26919 or seqfeature_id = 29555 or seqfeature_id = 31677 or seqfeature_id = 32569 or seqfeature_id = 35258 or seqfeature_id = 37374 or seqfeature_id = 38015 or seqfeature_id = 40665 or seqfeature_id = 42767 or seqfeature_id = 43518 or seqfeature_id = 46206 or seqfeature_id = 48320 or seqfeature_id = 49065 or seqfeature_id = 51534 or seqfeature_id = 53646 or seqfeature_id = 54586 or seqfeature_id = 59421 or seqfeature_id = 63250 or seqfeature_id = 64498 or seqfeature_id = 69670 or seqfeature_id = 73420 or seqfeature_id = 74864 or seqfeature_id = 77515 or seqfeature_id = 79511 or seqfeature_id = 80113 or seqfeature_id = 82919 or seqfeature_id = 84966 or seqfeature_id = 85659 or seqfeature_id = 88313 or seqfeature_id = 90885 or seqfeature_id = 91695 or seqfeature_id = 94136 or seqfeature_id = 96252 or seqfeature_id = 97046 or seqfeature_id = 100070 or seqfeature_id = 102225 or seqfeature_id = 103083 or seqfeature_id = 105619 or seqfeature_id = 107850 or seqfeature_id = 108669 or seqfeature_id = 111747 or seqfeature_id = 113931 or seqfeature_id = 114745 or seqfeature_id = 117413 or seqfeature_id = 119579 or seqfeature_id = 120398 or seqfeature_id = 123021 or seqfeature_id = 125112 or seqfeature_id = 125762 or seqfeature_id = 128535 or seqfeature_id = 130625 or seqfeature_id = 131969 or seqfeature_id = 132681 or seqfeature_id = 134865 or seqfeature_id = 136925 or seqfeature_id = 139364 or seqfeature_id = 141515 or seqfeature_id = 142247 or seqfeature_id = 144871 or seqfeature_id = 147141 or seqfeature_id = 147810 or seqfeature_id = 150229 or seqfeature_id = 152390 or seqfeature_id = 153124 or seqfeature_id = 155794 or seqfeature_id = 157928 or seqfeature_id = 158649 or seqfeature_id = 161147 or seqfeature_id = 163095 or seqfeature_id = 163758 or seqfeature_id = 166302 or seqfeature_id = 168468 or seqfeature_id = 169210 or seqfeature_id = 171763 or seqfeature_id = 173798 or seqfeature_id = 174619 or seqfeature_id = 177203 or seqfeature_id = 179332 or seqfeature_id = 180024 or seqfeature_id = 182543 or seqfeature_id = 184551 or seqfeature_id = 185251 or seqfeature_id = 187769 or seqfeature_id = 189641 or seqfeature_id = 190384 or seqfeature_id = 193146 or seqfeature_id = 195184 or seqfeature_id = 195913 or seqfeature_id = 198559 or seqfeature_id = 200755 or seqfeature_id = 201555 or seqfeature_id = 203670 or seqfeature_id = 206266 or seqfeature_id = 207056 or seqfeature_id = 209562 or seqfeature_id = 211848 or seqfeature_id = 212543 or seqfeature_id = 216940 or seqfeature_id = 217729 or seqfeature_id = 220435 or seqfeature_id = 222742 or seqfeature_id = 223873 or seqfeature_id = 229383 or seqfeature_id = 233819 or seqfeature_id = 235219 or seqfeature_id = 237679 or seqfeature_id = 239703 or seqfeature_id = 240509 or seqfeature_id = 243450 or seqfeature_id = 245850 or seqfeature_id = 246731 or seqfeature_id = 249502 or seqfeature_id = 251461 or seqfeature_id = 252167 or seqfeature_id = 254685 or seqfeature_id = 256915 or seqfeature_id = 257658 or seqfeature_id = 260188 or seqfeature_id = 262420 or seqfeature_id = 263251 or seqfeature_id = 265711 or seqfeature_id = 267737 or seqfeature_id = 268434 or seqfeature_id = 270898 or seqfeature_id = 272922 or seqfeature_id = 273692 or seqfeature_id = 276162 or seqfeature_id = 278188 or seqfeature_id = 279011 or seqfeature_id = 281475 or seqfeature_id = 283503 or seqfeature_id = 284281 or seqfeature_id = 286745 or seqfeature_id = 288775 or seqfeature_id = 289545 or seqfeature_id = 292005 or seqfeature_id = 294031 or seqfeature_id = 294801 or seqfeature_id = 297267 or seqfeature_id = 299293)
inner join seqfeature on seqfeature_qualifier_value.seqfeature_id = seqfeature.seqfeature_id
inner join bioentry on seqfeature.bioentry_id = bioentry.bioentry_id
inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id and biodatabase.name = "chlamydiales"
group by seqfeature_qualifier_value.seqfeature_id;
"""

