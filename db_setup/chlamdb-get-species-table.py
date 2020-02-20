#! /usr/bin/env python


def median_RBBH2species(biodb, cutoff=97):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table if not exists taxid2species (taxon_id INT, species_id INT)'
    server.adaptor.execute(sql,)
    server.commit()

    sql_taxon = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % biodb
    taxon_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_taxon,)]

    sql2 = 'select taxon_1,taxon_2,median_identity from comparative_tables_shared_og_av_id;'

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
            sql = 'insert into taxid2species values(%s,%s)' % (taxon, species_index)
            server.adaptor.execute(sql,)
        server.commit()
        species_index+=1


def bioentry_metadata(biodb):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio import Entrez
    
    Entrez.email = "trestan.pillonel@chuv.ch"
    
    server, db = manipulate_biosqldb.load_db(biodb)

    #sql = 'create table if not exists bioentry_metadata%s (taxon_id INT, species_id INT)' % biodb
    #server.adaptor.execute(sql,)
    #server.commit()

    sql_accessions = 'select accession, bioentry_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                ' where t1.name="%s"' % biodb
                
    accession2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_accessions,))

    sql = 'create table assembly_metadata (assembly_id INT AUTO_INCREMENT PRIMARY KEY, AssemblyAccession varchar(200), ReleaseLevel TEXT, PartialGenomeRepresentation varchar(20), LastUpdateDate varchar(200), ' \
          ' RefSeq_category varchar(40), SpeciesTaxid INTEGER, BioSampleAccn varchar(200), ' \
          ' SubmitterOrganization TEXT, Taxid INTEGER, FtpPath_GenBank TEXT, ExclFromRefSeq TEXT, AsmReleaseDate_GenBank varchar(400), ' \
          ' FtpPath_RefSeq TEXT, AssemblyName varchar(400), SpeciesName TEXT, AnomalousList TEXT, AssemblyStatus varchar(200), Coverage INTEGER, Organism TEXT)'
    print(sql)
    server.adaptor.execute(sql,)
    
    sql2 = 'create table bioentry2assembly (bioentry_id INTEGER, assembly_id INTEGER);'
    server.adaptor.execute(sql2,)
    
    assembly2data = {}
    bioentry_id2assembly_accession = {}
    for n, accession in enumerate(accession2bioentry_id):
        print("%s / %s -- %s" % (n, len(accession2bioentry_id), accession))
        try:
            handle1 = Entrez.esearch(db="nuccore", term="%s[Accession]" % accession)
            record1 = Entrez.read(handle1)
            ncbi_id = record1['IdList'][0]
        except IndexError:
            #accession2assembly_accession[accession2bioentry_id[accession]] = None
            continue

        handle2 = Entrez.elink(dbfrom="nuccore", db="assembly", id=ncbi_id)
        record2 = Entrez.read(handle2)
        try:
            id = record2[0]['LinkSetDb'][0]['Link'][0]['Id']
        except:
            print("No assembly link for %s --> try from biosample" % accession)
            handle_biosample = Entrez.elink(dbfrom="nuccore", db="biosample", id=ncbi_id)
            record_biosample = Entrez.read(handle_biosample)
            id_sample = record_biosample[0]['LinkSetDb'][0]['Link'][0]['Id']
            handle2 = Entrez.elink(dbfrom="biosample", db="assembly", id=id_sample)
            record2 = Entrez.read(handle2)
            id = record2[0]['LinkSetDb'][0]['Link'][0]['Id']

        handle3 = Entrez.esummary(db='assembly',id=id, retmode='xml')
        record = Entrez.read(handle3, validate=False)
        # ['GB_BioProjects', 
        # 'AsmReleaseDate_RefSeq', 'ReleaseLevel', 'PartialGenomeRepresentation', 'LatestAccession', 'SortOrder', 'Primary', 
        # 'FtpPath_Stats_rpt', 'UCSCName', 'Synonym', 'RsUid', 'WGS', 'LastMajorReleaseAccession', 'LastUpdateDate', 
        # 'AssemblyDescription', 'Biosource', 'AssemblyAccession', 'RefSeq_category', 'ReleaseType', 'SpeciesTaxid', 
        # 'BioSampleAccn', 'AssemblyType', 'AsmUpdateDate', 'AsmReleaseDate_GenBank', 'SubmitterOrganization', 'AssemblyClass',
        # 'ChainId', 'Taxid', 'FtpPath_GenBank', 'SubmissionDate', 'RS_BioProjects', 'FtpPath_Regions_rpt', 'PropertyList', 'SeqReleaseDate', 
        # 'ExclFromRefSeq', 'GbUid', 'FtpPath_RefSeq', 'AssemblyName', 'ContigN50', 'SpeciesName', 'AnomalousList', 'FromType', 
        # 'BioSampleId', 'AssemblyStatus', 'Coverage', 'Organism', 'FtpPath_Assembly_rpt', 'GB_Projects', 'Meta', 'EnsemblName', 'ScaffoldN50', 'RS_Projects']

        assembly_accession =  record['DocumentSummarySet']['DocumentSummary'][0]["AssemblyAccession"]
        bioentry_id2assembly_accession[accession2bioentry_id[accession]] = assembly_accession
        
        if assembly_accession not in assembly2data:
            assembly2data[assembly_accession] = {}
            assembly2data[assembly_accession]["ReleaseLevel"] = record['DocumentSummarySet']['DocumentSummary'][0]["ReleaseLevel"]
            assembly2data[assembly_accession]["PartialGenomeRepresentation"] = record['DocumentSummarySet']['DocumentSummary'][0]["PartialGenomeRepresentation"]
            assembly2data[assembly_accession]["LastUpdateDate"] = record['DocumentSummarySet']['DocumentSummary'][0]["LastUpdateDate"]
            
            assembly2data[assembly_accession]["RefSeq_category"] = record['DocumentSummarySet']['DocumentSummary'][0]["RefSeq_category"]
            #assembly2data[assembly_accession]["ReleaseType"] = record['DocumentSummarySet']['DocumentSummary'][0]["ReleaseType"]
            assembly2data[assembly_accession]["SpeciesTaxid"] = record['DocumentSummarySet']['DocumentSummary'][0]["SpeciesTaxid"]
            assembly2data[assembly_accession]["BioSampleAccn"] = record['DocumentSummarySet']['DocumentSummary'][0]["BioSampleAccn"]
            assembly2data[assembly_accession]["SubmitterOrganization"] = record['DocumentSummarySet']['DocumentSummary'][0]["SubmitterOrganization"]
            assembly2data[assembly_accession]["Taxid"] = record['DocumentSummarySet']['DocumentSummary'][0]["Taxid"]
            assembly2data[assembly_accession]["FtpPath_GenBank"] = record['DocumentSummarySet']['DocumentSummary'][0]["FtpPath_GenBank"]
            assembly2data[assembly_accession]["ExclFromRefSeq"] = record['DocumentSummarySet']['DocumentSummary'][0]["ExclFromRefSeq"]
            assembly2data[assembly_accession]["AsmReleaseDate_GenBank"] = record['DocumentSummarySet']['DocumentSummary'][0]["AsmReleaseDate_GenBank"]
            assembly2data[assembly_accession]["FtpPath_RefSeq"] = record['DocumentSummarySet']['DocumentSummary'][0]["FtpPath_RefSeq"]
            assembly2data[assembly_accession]["AssemblyName"] = record['DocumentSummarySet']['DocumentSummary'][0]["AssemblyName"]
            assembly2data[assembly_accession]["SpeciesName"] = record['DocumentSummarySet']['DocumentSummary'][0]["SpeciesName"]
            assembly2data[assembly_accession]["AnomalousList"] = record['DocumentSummarySet']['DocumentSummary'][0]["AnomalousList"]
            #assembly2data[assembly_accession]["BioSampleId"] = record['DocumentSummarySet']['DocumentSummary'][0]["BioSampleId"])
            assembly2data[assembly_accession]["AssemblyStatus"] = record['DocumentSummarySet']['DocumentSummary'][0]["AssemblyStatus"]
            assembly2data[assembly_accession]["Coverage"] = record['DocumentSummarySet']['DocumentSummary'][0]["Coverage"]
            assembly2data[assembly_accession]["Organism"] = record['DocumentSummarySet']['DocumentSummary'][0]["Organism"]
                
    for assembly_accession in assembly2data:
        print(assembly_accession)
        sql = 'insert into assembly_metadata (AssemblyAccession, ReleaseLevel, PartialGenomeRepresentation, LastUpdateDate, RefSeq_category, ' \
              ' SpeciesTaxid, BioSampleAccn, SubmitterOrganization, Taxid, FtpPath_GenBank, ExclFromRefSeq, AsmReleaseDate_GenBank, FtpPath_RefSeq,' \
              ' AssemblyName, SpeciesName, AnomalousList, AssemblyStatus,Coverage, Organism) ' \
              ' values (%s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)' 
              
        values = [assembly_accession,
                  assembly2data[assembly_accession]["ReleaseLevel"] ,
                  assembly2data[assembly_accession]["PartialGenomeRepresentation"] ,
                  assembly2data[assembly_accession]["LastUpdateDate"],
                  assembly2data[assembly_accession]["RefSeq_category"] ,
                  assembly2data[assembly_accession]["SpeciesTaxid"],
                  assembly2data[assembly_accession]["BioSampleAccn"] ,
                  assembly2data[assembly_accession]["SubmitterOrganization"],
                  assembly2data[assembly_accession]["Taxid"] ,
                  assembly2data[assembly_accession]["FtpPath_GenBank"],
                  assembly2data[assembly_accession]["ExclFromRefSeq"],
                  assembly2data[assembly_accession]["AsmReleaseDate_GenBank"] ,
                  assembly2data[assembly_accession]["FtpPath_RefSeq"],
                  assembly2data[assembly_accession]["AssemblyName"],
                  assembly2data[assembly_accession]["SpeciesName"] ,
                  assembly2data[assembly_accession]["AnomalousList"] ,
                  assembly2data[assembly_accession]["AssemblyStatus"],
                  assembly2data[assembly_accession]["Coverage"],
                  assembly2data[assembly_accession]["Organism"]]
        server.adaptor.execute(sql, values)
    server.commit()
    
    sql = 'select AssemblyAccession, assembly_id from assembly_metadata'
    assembly_accession2assembly_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    for bioentry_id in bioentry_id2assembly_accession:
        assembly_accession = bioentry_id2assembly_accession[bioentry_id]
        assembly_id = assembly_accession2assembly_id[assembly_accession]
        sql = 'insert into bioentry2assembly values (%s, %s)' % (bioentry_id, assembly_id)
        server.adaptor.execute(sql,)
    server.commit()


def create_species_curated_taxonomy(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, fb = manipulate_biosqldb.load_db(biodb)
    
    sql1 = 'select distinct species_id from taxid2species'
    species_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1,)]
    # get mot frequent classification for each species defined at 97% median identity
    sql = 'select t5.species_id,t3.phylum,t3.order,t3.family,t3.genus,t3.species, count(*) as n from bioentry2assembly t1 ' \
          ' inner join assembly_metadata t2 on t1.assembly_id=t2.assembly_id ' \
          ' inner join blastnr_blastnr_taxonomy t3 on t2.taxid=t3.taxon_id ' \
          ' inner join bioentry t4 on t1.bioentry_id=t4.bioentry_id ' \
          ' inner join taxid2species t5 on t4.taxon_id=t5.taxon_id ' \
          ' group by t5.species_id,t3.phylum,t3.order,t3.family,t3.genus,t3.species order by n DESC;'
          
    data = server.adaptor.execute_and_fetchall(sql,)
    species_id2most_frequent_taxonomy = {}
    for row in data:
        if row[0] not in species_id2most_frequent_taxonomy:
            species_id2most_frequent_taxonomy[row[0]] = row
    
    sql = 'create table species_curated_taxonomy (species_id INTEGER, phylum varchar(200), `order` varchar(200), family varchar(200), genus varchar(200), species varchar(600))'
    server.adaptor.execute(sql)
    for species_id in species_id_list:
        if species_id not in species_id2most_frequent_taxonomy:
            sql = 'select description from bioentry t1 inner join taxid2species t2 on t1.taxon_id=t2.taxon_id where species_id=%s and description not like "%%%%plasmid%%%%" limit 1;' % (species_id)
            description = server.adaptor.execute_and_fetchall(sql,)[0][0]
            sql = 'insert into species_curated_taxonomy values (%s,%s,%s,%s,%s,%s)'
            values = [species_id, "-","-","-","-",description]
            server.adaptor.execute(sql, values)
        else:
            sql = 'insert into species_curated_taxonomy values (%s,%s,%s,%s,%s,%s)'
            values = [species_id,
                    species_id2most_frequent_taxonomy[species_id][1],
                    species_id2most_frequent_taxonomy[species_id][2],
                    species_id2most_frequent_taxonomy[species_id][3],
                    species_id2most_frequent_taxonomy[species_id][4],
                    species_id2most_frequent_taxonomy[species_id][5]]
            server.adaptor.execute(sql, values)
    server.commit()
    

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")


    args = parser.parse_args()
    median_RBBH2species(args.biodb)
    bioentry_metadata(args.biodb)
    create_species_curated_taxonomy(args.biodb)