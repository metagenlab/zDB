#!/usr/bin/python

def intero_data2features_string(interpro_data):

    template = '''

    ft2.addFeature({
        data: [{x:%s,y:%s,description:"%s: %s",id:"a%s"}],
        name: "%s",
        className: "%s",
        color: "#81BEAA",
        type: "rect",
        filter: "type2"
    });

        '''

    merged_features = ''

    for interpro_entry in interpro_data:
        # signature_accession, signature_description,start,stop
        start = interpro_entry[2]
        stop = interpro_entry[3]
        accession = interpro_entry[0]
        description = interpro_entry[1]
        merged_features+=template % (start,
                                           stop,
                                           accession,
                                           description,
                                           accession,
                                           accession,
                                           accession)
    return merged_features

import manipulate_biosqldb
server, db = manipulate_biosqldb.load_db("chlamydia_04_16")
sql_pfam = 'select signature_accession, signature_description,start,stop' \
           ' from interpro_%s where locus_tag="%s" ' \
           ' and analysis="Pfam";' % ("chlamydia_04_16", "WCW_RS07525")
pfam_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pfam, )]

print intero_data2features_string(pfam_data)
