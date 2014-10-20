#!/usr/bin/python

import heatmap
import numpy as np
import manipulate_biosqldb


# heatmap Chlamydiales pan-genome

server, db = manipulate_biosqldb.load_db("saureus1")
matrix = np.array(manipulate_biosqldb.get_orthology_table(server, "saureus1"))

taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, "saureus1")
print taxon_id2description

group_names = matrix[:,0]
taxons_ids = manipulate_biosqldb.get_taxon_id_list(server, "saureus1")
taxons_ids = [taxon_id2description[i] for i in taxons_ids]
print taxons_ids


M = matrix[:, 1:].astype(float)
M = heatmap.randomize_table(M)


heatmap.heatmap(M, output="heat_saureus.pdf", breaks="-0.5, 0.5, 1.5, 2.5", rows=False, columns = taxons_ids)
