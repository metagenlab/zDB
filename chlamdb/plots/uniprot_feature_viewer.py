#!/usr/bin/python

def phobius_data2features_string(phobius_data):

    template = '''

    ft2.addFeature({
        data: [{x:%s,y:%s,description:"%s",id:"a%s"}],
        name: "%s",
        className: "%s",
        color: "%s",
        type: "rect",
        filter: "type2"
    });

        '''
    template_TM = '''

    ft2.addFeature({
        data: [%s],
        name: "%s",
        className: "%s",
        color: "%s",
        type: "rect",
        filter: "type2"
    });

        '''
    one_TM = '{x:%s,y:%s,description:"%s",id:"a%s"},'
    all_TM = ''
    merged_features = ''

    for phobius_entry in phobius_data:
        # signature_accession, signature_description,start,stop
        start = phobius_entry[1]
        stop = phobius_entry[2]
        accession = phobius_entry[0]

        if accession == 'TRANSMEMBRANE':
            accession = 'Transmembrane'
            all_TM+=one_TM % (start, stop, accession, accession)
        else:
            col = '#FE2E2E'
            accession = 'S_peptide'

            merged_features+=template % (start,
                                               stop,
                                               accession,
                                               accession,
                                               accession,
                                               accession,
                                               col)

    if len(all_TM)>0:
        merged_features+=template_TM % (all_TM[0:-1],
                                     'TM',
                                     'TM',
                                     '#FE9A2E')
    return merged_features
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


def superfamily_data2features_string(superfamily_data):

    template = '''

    ft2.addFeature({
        data: [{x:%s,y:%s,description:"%s: %s",id:"a%s"}],
        name: "%s",
        className: "%s",
        color: "#5882FA",
        type: "rect",
        filter: "type2"
    });

        '''

    merged_features = ''

    for superfamily_entry in superfamily_data:
        # signature_accession, signature_description,start,stop
        start = superfamily_entry[2]
        stop = superfamily_entry[3]
        accession = superfamily_entry[0]
        description = superfamily_entry[1]
        merged_features+=template % (start,
                                           stop,
                                           accession,
                                           description,
                                           accession,
                                           accession,
                                           accession)
    return merged_features







'''
from chlamdb.biosqldb import manipulate_biosqldb
server, db = manipulate_biosqldb.load_db("chlamydia_04_16")
sql_pfam = 'select signature_accession, signature_description,start,stop' \
           ' from interpro where locus_tag="%s" ' \
           ' and analysis="Pfam";' % ("chlamydia_04_16", "WCW_RS07525")
#pfam_data = [list(i) for i in server.adaptor.execute_and_fetchall(sql_pfam, )]

#print intero_data2features_string(pfam_data)
sql18 = 'select signature_accession,start,stop from interpro where analysis="Phobius" and locus_tag="%s" ' \
        ' and signature_accession in ("TRANSMEMBRANE",' \
        ' "SIGNAL_PEPTIDE_C_REGION","SIGNAL_PEPTIDE_N_REGION");' % ("chlamydia_04_16", "BN1013_01691")


#phobius_data = server.adaptor.execute_and_fetchall(sql18, )
#print phobius_data2features_string(phobius_data)
'''
