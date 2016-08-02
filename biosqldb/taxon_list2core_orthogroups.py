#!/usr/bin/python

import biosql_own_sql_tables
import manipulate_biosqldb



accession_list = ["NZ_CP012300",
"NC_022566",
"NZ_CP006738",
"NZ_LN824139",
"NZ_CP011315",
"NZ_CP009208",
"NZ_CP014005",
"NZ_CP006799",
"NZ_CP012754",
"NC_018522",
"NZ_CP003996",
"KpGe",
"NZ_CP012743",
"NZ_CP012744",
"NZ_CP012993",
"NZ_CP012989",
"NZ_CP009461",
"NZ_CP012884",
"NZ_FO834904",
"NZ_CP006722",
"NZ_CP014009",
"NZ_CP013712",
"NC_017541",
"NC_006625",
"NZ_CP014011",
"NZ_CP009864",
"NZ_CP007732",
"NZ_CP008930",
"NZ_CP009877",
"NC_009649",
"NZ_CP011623",
"NZ_CP013324",
"NZ_CP006663",
"NC_022078",
"NZ_CP011577",
"NZ_CP010394",
"NZ_CP008798",
"NZ_CP010574",
"NZ_CP009776",
"NZ_CP011979",
"NZ_CP006927",
"NZ_CP007728",
"NZ_CP008828",
"NZ_CP009116",
"NC_016838",
"NZ_CP011643",
"NZ_CP009772",
"NZ_CP006919",
"NZ_CP011984",
"NZ_CP011986",
"NZ_CP011991",
"NZ_CP009873",
"NZ_CP008832"]

print len(accession_list)

taxon_list_chlamydiae = [
"291",
"67",
"1279767",
"1279774",
"1279496",
"48",
"46",
"55",
"87925",
"1279815",
"62",
"1279822",
"66",
"52",
"49",
"64",
"60",
"804807",
"886707",
"283",
"314",
"1069693",
"1069694",
"1137444",
"1143376",
"313",
"1172027",
"1172028",
"1035343",
"307",
"293",
"1279839",
"1279497"
]

'''
server, db = manipulate_biosqldb.load_db("chlamydia_04_16")

accession2taxon = manipulate_biosqldb.accession2taxon_id(server, "chlamydia_04_16")

taxon_list = []
for i in accession_list:
    taxon_list.append(accession2taxon[i])

cosson_taxons = list(set(accession2taxon.values()))
'''
biosql_own_sql_tables.taxon_subset2core_orthogroups("chlamydia_04_16", taxon_list_chlamydiae, "protein")
