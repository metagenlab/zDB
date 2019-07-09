#! /usr/bin/env python


def median_RBBH2species(biodb, cutoff=97):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table if not exists species_%s (taxon_id INT, species_id INT)' % biodb
    server.adaptor.execute(sql,)
    server.commit()

    sql_taxon = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % biodb
    taxon_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_taxon,)]

    sql2 = 'select taxon_1,taxon_2,median_identity from comparative_tables.shared_og_av_id_%s;' % biodb

    taxon2taxon2identity = {}
    for row in server.adaptor.execute_and_fetchall(sql2,):
        if row[0] not in taxon2taxon2identity:
            taxon2taxon2identity[row[0]] = {}
            taxon2taxon2identity[row[0]][row[1]] = row[2]
        else:
            taxon2taxon2identity[row[0]][row[1]] = row[2]
        if row[1] not in taxon2taxon2identity:
            taxon2taxon2identity[row[1]] = {}
            taxon2taxon2identity[row[1]][row[0]] = row[2]
        else:
            taxon2taxon2identity[row[1]][row[0]] = row[2]
    #print taxon2taxon2identity
    species_index = 0
    taxons_classified = []
    for taxon_1 in taxon_id_list:
        if taxon_1 in taxons_classified:
            continue
        species_list = [taxon_1]
        for taxon_2 in taxon_id_list:
            if taxon_1 == taxon_2:
                continue
            else:
                if float(taxon2taxon2identity[taxon_1][taxon_2]) >= cutoff:
                    species_list.append(taxon_2)
        for taxon in species_list:
            taxons_classified.append(taxon)
            sql = 'insert into species_%s values(%s,%s)' % (biodb, taxon, species_index)
            server.adaptor.execute(sql,)
        server.commit()
        species_index+=1
        

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")


    args = parser.parse_args()
    median_RBBH2species(args.biodb)