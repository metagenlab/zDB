#!/usr/bin/env python

from Bio import Entrez
from bs4 import BeautifulSoup
import pandas

Entrez.email = "trestan.pillonel@unil.ch"


def check_sra(sra_run_list, sra_sample_table):

    s = open('sample_remove.tsv', 'w')
    t = open('sra_sample_check.tsv', 'w')
    u = open('sra_check_detail.tsv', 'w')

    run_accession2sample_name = sra_sample_table.set_index("Run").to_dict()[
        "SampleName"]
    run_accession2LibraryLayout = sra_sample_table.set_index("Run").to_dict()[
        "LibraryLayout"]
    run_accession2ScientificName = sra_sample_table.set_index("Run").to_dict()[
        "ScientificName"]

    t.write('Run\tSampleName\tLibraryLayout\tScientificName\n')
    u.write(
        'project\texperiment\tsample\trun_number\tlargest_run\tlargest_size_mb\treplace\n')

    for run_accession in sra_run_list:
        filter = '%s[Accession]' % run_accession
        handle = Entrez.esearch(db="sra", term=filter, retmax=1000000)
        record = Entrez.read(handle)
        sra_id = record['IdList'][0]
        handle = Entrez.esummary(db="sra", id=sra_id, retmax=1)
        data = Entrez.read(handle)
        soup_exp = BeautifulSoup(data[0]['ExpXml'], 'html.parser')

        # library_layout
        # bioproject
        # biosample
        # sample acc
        # experiment
        experiment = soup_exp.find('experiment').attrs['acc']
        bioproject = soup_exp.find('bioproject').text
        sample = soup_exp.find('sample').attrs['acc']
        soup = BeautifulSoup(data[0]['Runs'], 'html.parser')
        largest = [0, '-']
        run_acc2data = {}
        for n, run in enumerate(soup.findAll('run')):
            run_size = round(float(run.attrs['total_bases']) / 1000000, 2)
            run_acc2data[run.attrs['acc']] = [bioproject,
                                              experiment, sample, n, run.attrs['acc'], run_size]
            if run_size > largest[0]:
                largest = [run_size, run.attrs['acc']]
        if largest[1] != run_accession:
            u.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (run_acc2data[largest[1]][0],
                                                      run_acc2data[largest[1]][1],
                                                      run_acc2data[largest[1]][2],
                                                      n + 1,
                                                      run_acc2data[largest[1]][4],
                                                      run_acc2data[largest[1]][5],
                                                      'replace %s (%sMb)' % (run_accession,
                                                                             run_acc2data[run_accession][5])))
            s.write("%s\n" % (run_accession2sample_name[run_accession]))
        else:
            u.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (run_acc2data[largest[1]][0],
                                                      run_acc2data[largest[1]][1],
                                                      run_acc2data[largest[1]][2],
                                                      n + 1,
                                                      run_acc2data[largest[1]][4],
                                                      run_acc2data[largest[1]][5],
                                                      '-'))
        t.write("%s\t%s\t%s\t%s\n" % (largest[1],
                                      run_accession2sample_name[run_accession],
                                      run_accession2LibraryLayout[run_accession],
                                      run_accession2ScientificName[run_accession]))
    s.close()
    t.close()
    u.close()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--sra_sample_file',
                        default=False, type=str, help="sra sample file")

    args = parser.parse_args()

    sra_table = pandas.read_csv(args.sra_sample_file, sep="\t")
    sra_accession_list = list(sra_table['Run'])

    sra_exp = sum([1 for i in sra_accession_list if 'SRX' in i])
    if sra_exp > 0:
        raise IOError('%s SRX accession found in the Run column' % sra_exp)
    else:
        check_sra(sra_accession_list,
                  sra_table)
