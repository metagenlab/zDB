#!/usr/bin/env python

import pygenprop
from pygenprop.database_file_parser import parse_genome_properties_flat_file


groups = ['METAPATH', 'CATEGORY', 'GUILD']

geneprop2category = {}

tree = parse_genome_properties_flat_file(open("/home/tpillone/tmp/genomeProperties.txt", "r"))
print("genprop.id\thierarchy\tdescription\treferences\tgenprop.type\tgenprop.name\tn_step\tn_element\tn_evidence\tname\treq_step\treq_functional_element\tid_evidence\tip_id\tsuff")
for genprop in tree:
    for n_step, step in enumerate(genprop.steps):
        for n_element, functional_element in enumerate(step.functional_elements):
            for n_evidence, evidence in enumerate(functional_element.evidence):
                ip_id_list = evidence.interpro_identifiers
                if len(ip_id_list)> 1:
                    print("More than one!!!!!",ip_id_list)
                    import sys 
                    sys.exit()
                elif len(ip_id_list) == 0:
                    ip_id = None
                else:
                    ip_id = ip_id_list[0]
                suff_evidence = evidence.sufficient
                if len(evidence.property_identifiers) > 0:
                    id_evidence = evidence.property_identifiers[0]
                else:
                    id_evidence = None

                if genprop.type in groups and id_evidence is not None:
                    if id_evidence not in geneprop2category:
                        geneprop2category[id_evidence] = [genprop.id]
                    else:
                        geneprop2category[id_evidence].append(genprop.id)
                id_functional_element = functional_element.id
                name = functional_element.name 
                req_functional_element = functional_element.required
                req_step = step.required
                if genprop.id in geneprop2category:
                    hier = geneprop2category[genprop.id]
                else:
                    hier = ''
                print(f"{genprop.id}\t{hier}\t{genprop.description}\t{genprop.references}\t{genprop.type}\t{genprop.name}\t{n_step}\t{n_element}\t{n_evidence}\t{id_functional_element}\t{name}\t{req_step}\t{req_functional_element}\t{id_evidence}\t{ip_id}\t{suff_evidence}")

#tree['GenProp0002'].steps[1].functional_elements[1].evidence[2].interpro_identifierss