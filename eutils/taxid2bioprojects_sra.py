#!/usr/bin/env python

from Bio import Entrez, SeqIO
import os
import re
from ftplib import FTP
import download_from_ftp

Entrez.email = "trestan.pillonel@unil.ch"

def gi(ncbi_term, database="nuccore", retmax=20):
    handle = Entrez.esearch(db=database, term=ncbi_term, retmax=retmax)
    record = Entrez.read(handle)
    return record["IdList"]

def bioproject2data(project_uid):

    # get project summary data
    handle_project = Entrez.esummary(db="bioproject", id=project_uid)
    project_summary_record = Entrez.read(handle_project, validate=False)

    description = project_summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Description']
    project_data_type = project_summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Data_Type']
    poject_name = project_summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Name']
    project_title = project_summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Title']
    registration_date = project_summary_record['DocumentSummarySet']['DocumentSummary'][0]['Registration_Date']
    project_id = project_summary_record['DocumentSummarySet']['DocumentSummary'][0]['Project_Id']

    return project_id, poject_name, project_title, description, registration_date

def chunks(l, n):
    "return sublists of l of minimum length n (work subdivision for the subprocesing module"
    return [l[i:i+n] for i in range(0, len(l), n)]

def merge_dictionnaries(dico1, dico2):
    '''

    merge two dictionnaries with lists as value
    if key already present in the fist dictionnary, add values to the existing list

    :param dico1:
    :param dico2:
    :return:
    '''
    for key in dico2:
        if key in dico1:
            dico1[key]+=dico2[key]
        else:
            dico1[key] = dico2[key]
    return dico1
def get_genomic_projects(taxon_id):

    filter = 'txid%s[Organism:exp] AND "strategy wgs"[Properties]' % taxon_id
    handle = Entrez.esearch(db="sra", term=filter, retmax=1000000)
    record = Entrez.read(handle)
    sra_ids = record['IdList']
    chunk_lists = chunks(sra_ids, 1000)
    project_id2sra_list = {}
    for n, sra_list in enumerate(chunk_lists):
        if n%100 == 0:
            print ('%s/%s' % (n, len(chunk_lists)))
        temp_dico = sra2bioproject2(sra_list)
        #print 'n projects, n sra:', len(temp_dico), sum([len(i) for i in temp_dico.values()])
        project_id2sra_list = merge_dictionnaries(project_id2sra_list, temp_dico)
    return project_id2sra_list

def sra2bioproject(sra_id_list):
    # get assembly links
    filter_id = ','.join([str(i) for i in sra_id_list])
    handle = Entrez.elink(dbfrom="sra", db="bioproject", id=filter_id, retmax=100)
    links = Entrez.read(handle)
    project_ids = [i['Id'] for i in links[0]['LinkSetDb'][1]['Link']]
    sra_ids = links[0]['IdList']
    project_id2sra_list={}
    for project_id, sra_id in zip(project_ids, sra_ids):
        if project_id not in project_id2sra_list:
            project_id2sra_list[project_id] = [sra_id]
        else:
            project_id2sra_list[project_id].append(sra_id)
    return project_id2sra_list

def sra2bioproject2(sra_id_list):
    # get assembly links
    import re
    filter_id = ','.join([str(i) for i in sra_id_list])
    handle = Entrez.esummary(db="sra", id=filter_id, retmax=1000)
    links = Entrez.read(handle)
    project_id2sra_list = {}
    for i in links:
        data = i['ExpXml']
        try:
            bioproject = re.findall('<Bioproject>(.*)</Bioproject>',data)[0]
        except:
            print ('problem:', data)
            continue
        id = i['Id']
        if bioproject not in project_id2sra_list:
            project_id2sra_list[bioproject] = [id]
        else:
            project_id2sra_list[bioproject].append(id)

    return project_id2sra_list

'''
print '%s\t%s\t%s\t%s' % (Project_Id,
                          Project_Title,
                          len(id_list),
                          len(id_list_sra))
'''

test = \
    [4506179, 4506178]

#print sra2bioproject2(test)

project2sra_list = get_genomic_projects("2")

print ('n projects:', len(project2sra_list))
print ('n sra:', sum([len(i) for i in project2sra_list.values()]))
with open('sra_projects_bacteria.tab', 'w') as f:
    for n, project in enumerate(project2sra_list):
        print ("%s\t%s\t%s" % (project, n, len(project2sra_list)))
        if len(project2sra_list[project]) < 20:
            continue
        project_id, poject_name, project_title, description, registration_date = bioproject2data(gi(project, database="bioproject")[0])
        line = "%s\t%s\t%s\t%s\t%s\n" % (project, str(len(project2sra_list[project])),
                                  registration_date,
                                  poject_name,
                                  project_title)
        f.write(line.encode("utf-8"))
