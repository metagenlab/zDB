#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import os
import shutil
import generate_bsub_file
import splitFasta
from Bio import SeqIO
from shell_command import *
import re


class Orthomcl():
    def __init__(self, input_fasta,
                 index_mysql_tables,
                 working_dir,
                 BLAST_file,
                 evalue_exponent,
                 percent_match_cutoff,
                 inflation_value,
                 local=False):

        import sys
        evalue_string = "%f" % 10 ** int(evalue_exponent)

        if local and not BLAST_file:
            raise IOError('Blast file should be provided in case of local execution')

        if not BLAST_file:
            self.generate_config_file(index_mysql_tables, working_dir, str(percent_match_cutoff), str(evalue_exponent))
            sys.stdout.write("install shchema\n")
            ortho_install_shema = "orthomclInstallSchema orthomcl_config.txt install_schema.log"

            sterr, stout, process = shell_command(ortho_install_shema)
            print sterr
            sys.stdout.write(stout + "\n")
            sys.stdout.write("Adjust fasta\n")
            self.adjust_fasta(input_fasta)
            ortho_filter = "orthomclFilterFasta fasta 10 20"
            stdout_str, stderr_str, process = shell_command(ortho_filter)
            sys.stdout.write("Running BLAST\n")
            self.run_blast("goodProteins.fasta", evalue_string)
            sys.stdout.write("orthomclBlastParser\n")
            ortho_blast_parser = "orthomclBlastParser merged_blast.txt fasta >> similarSequences.txt"
            script = generate_bsub_file.BSUB_script(command=ortho_blast_parser, mem_in_GB=2, name="parsing",
                                                    log_file="out_BlastParser.txt", error_file="err_BlastParser.txt")
            generate_bsub_file.run_job_and_wait(script)


        else:
            self.generate_config_file(index_mysql_tables, working_dir, str(percent_match_cutoff), str(evalue_exponent),
                                      local)

            sys.stdout.write("install shchema\n")
            ortho_install_shema = "orthomclInstallSchema orthomcl_config.txt install_schema.log"
            sterr, stout, process = shell_command(ortho_install_shema)
            print sterr, stout, process
            sys.stdout.write("Parsing %s BLAST file\n" % BLAST_file)
            ortho_blast_parser = "orthomclBlastParser %s fasta >> similarSequences.txt" % (BLAST_file)
            if not local:
                script = generate_bsub_file.BSUB_script(command=ortho_blast_parser, mem_in_GB=2, name="parsing",
                                                        log_file="out_BlastParser.txt",
                                                        error_file="err_BlastParser.txt")
                generate_bsub_file.run_job_and_wait(script)
            else:
                sterr, stout, process = shell_command(ortho_blast_parser)
                print sterr, stout, process


        if not local:
            sys.stdout.write("orthomclLoadBlast \n")
            ortho_load_blast = "orthomclLoadBlast orthomcl_config.txt similarSequences.txt"
            script = generate_bsub_file.BSUB_script(command=ortho_load_blast, mem_in_GB=2, name="load blast",
                                                    log_file="out_LoadBlast.txt", error_file="err_LoadBlast.txt")
            generate_bsub_file.run_job_and_wait(script)

            sys.stdout.write("orthomclPairs\n")
            ortho_pairs = "orthomclPairs orthomcl_config.txt orthomcl_pairs.log cleanup=yes"
            script = generate_bsub_file.BSUB_script(command=ortho_pairs, mem_in_GB=2, name="pairs",
                                                    log_file="out_Pairs.txt", error_file="err_Pairs.txt")
            generate_bsub_file.run_job_and_wait(script)

            sys.stdout.write("orthomclDumpPairsFiles\n")
            ortho_dump = "orthomclDumpPairsFiles orthomcl_config.txt"

            script = generate_bsub_file.BSUB_script(command=ortho_dump, mem_in_GB=2, name="dump",
                                                    log_file="out_DumpPairsFiles.txt",
                                                    error_file="err_DumpPairsFiles.txt")
            generate_bsub_file.run_job_and_wait(script)

            sys.stdout.write("mcl\n")
            ortho_mcl = "mcl mclInput --abc -I %s -o mclOutput" % inflation_value
            script = generate_bsub_file.BSUB_script(command=ortho_mcl, mem_in_GB=2, name="mcl", log_file="out_mcl.txt",
                                                    error_file="err_mcl.txt")
            generate_bsub_file.run_job_and_wait(script)

            sys.stdout.write("drop database\n")
            cmd = "drop.sh orthomcl orthomcl orthomcl"
            script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=2, name="drop", log_file="out_drop.txt",
                                                    error_file="err_drop.txt")
            generate_bsub_file.run_job_and_wait(script)

            sys.stdout.write("groups\n")
            stdout_str, stderr_str, process = shell_command("orthomclMclToGroups group_ 0 < mclOutput > mclOutput.grp")

            sys.stdout.write("singletons\n")
            stdout_str, stderr_str, process = shell_command(
                "orthomclSingletons goodProteins.fasta mclOutput.grp > singletons.txt")

        else:


            sys.stdout.write("orthomclLoadBlast \n")
            ortho_load_blast = "orthomclLoadBlast orthomcl_config.txt similarSequences.txt"
            stdout_str, stderr_str, process = shell_command(ortho_load_blast)
            sys.stdout.write(stdout_str)
            sys.stderr.write(stderr_str)

            sys.stdout.write("orthomclPairs\n")
            ortho_pairs = "orthomclPairs orthomcl_config.txt orthomcl_pairs.log cleanup=yes"
            stdout_str, stderr_str, process = shell_command(ortho_pairs)
            sys.stdout.write(stdout_str)
            sys.stderr.write(stderr_str)

            sys.stdout.write("orthomclDumpPairsFiles\n")
            ortho_dump = "orthomclDumpPairsFiles orthomcl_config.txt"
            stdout_str, stderr_str, process = shell_command(ortho_dump)
            sys.stdout.write(stdout_str)
            sys.stderr.write(stderr_str)

            sys.stdout.write("mcl\n")
            ortho_mcl = "mcl mclInput --abc -I %s -o mclOutput" % inflation_value
            stdout_str, stderr_str, process = shell_command(ortho_mcl)
            sys.stdout.write(stdout_str)
            sys.stderr.write(stderr_str)

            # sys.stdout.write("drop database\n")
            # cmd = "drop.sh orthomcl orthomcl orthomcl"

            sys.stdout.write("groups\n")
            stdout_str, stderr_str, process = shell_command("orthomclMclToGroups group_ 0 < mclOutput > mclOutput.grp")
            sys.stdout.write(stdout_str)
            sys.stderr.write(stderr_str)

            sys.stdout.write("singletons\n")
            stdout_str, stderr_str, process = shell_command(
                "orthomclSingletons goodProteins.fasta mclOutput.grp > singletons.txt")
            sys.stdout.write(stdout_str)
            sys.stderr.write(stderr_str)

        if not os.path.exists("log_files"):
            os.mkdir("log_files")
        stdout_str, stderr_str, process = shell_command("mv out* log_files")
        stdout_str, stderr_str, process = shell_command("mv err* log_files")
        stdout_str, stderr_str, process = shell_command("mv *log log_files")

        if not os.path.exists("input_fasta"):
            os.mkdir("input_fasta")
        for one_fasta in input_fasta:
            stdout_str, stderr_str, process = shell_command("mv %s input_fasta" % one_fasta)
        stdout_str, stderr_str, process = shell_command("rm goodProteins.fasta.p*")

    def adjust_fasta(self, fasta_files):
        count = 0
        job_ids = []
        for one_file in fasta_files:
            count += 1
            # outname = re.sub(".faa", "", one_file)
            outname = "genome_%s" % str(count)
            fasta = [record for record in SeqIO.parse(open(one_file), "fasta")]
            fields = fasta[0].id.split("|")
            if len(fields) == 1:
                cmd = "orthomclAdjustFasta %s %s 1" % (outname, one_file)
            elif len(fields) == 5:
                cmd = "orthomclAdjustFasta %s %s 4" % (outname, one_file)
            elif len(fields) == 2:
                cmd = "orthomclAdjustFasta %s %s 2" % (outname, one_file)
            else:
                print fields
                raise NameError('Unknown fasta header format')
            script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=2, name="format",
                                                    log_file="out_AdjustFasta.txt", error_file="err_AdjustFasta.txt")
            job_ids.append(generate_bsub_file.run_job(script))
        sys.stdout.write(job_ids + "\n")
        generate_bsub_file.wait_multi_jobs(job_ids)
        if not os.path.exists("fasta"):
            os.mkdir("fasta")
        shell_command("mv *fasta fasta")

    def run_blast(self, fasta_file, evalue_string):
        n_files = splitFasta.split_fasta(fasta_file, 50)
        # formatdb = "formatdb -i %s  -p T" % fasta_file
        formatdb = 'makeblastdb -in %s -dbtype "prot"' % fasta_file
        script = generate_bsub_file.BSUB_script(command=formatdb, mem_in_GB=2, name="format",
                                                log_file="out_formatdb.txt", error_file="err_formatdb.txt")
        generate_bsub_file.run_job_and_wait(script)
        job_ids = []
        for i in range(1,
                       n_files + 1):  # old cutoff: 0.00001 | giant virus  0.01 | kleb 0.00001 | giant virus mai 15 : 0.1
            cmd = "blastp -db %s -query group_%s.fasta -out blast_result_%s.txt -outfmt 6 -num_threads 8 -evalue %s" % (
            fasta_file, i, i, evalue_string)
            script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=4, name="blast", log_file="out_blast.txt",
                                                    error_file="err_blast.txt")
            job_ids.append(generate_bsub_file.run_job(script))
        generate_bsub_file.wait_multi_jobs(job_ids)
        cmd = "for i in {1..%s}; do echo $i; cat blast_result_$i.txt>>merged_blast.txt;done" % (i)
        script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=4, name="concat", log_file="out_cat_fasta.txt",
                                                error_file="err_cat_fasta.txt")
        generate_bsub_file.run_job_and_wait(script)
        if not os.path.exists("blast_groups"):
            os.mkdir("blast_groups")
        shell_command("mv blast_result_* blast_groups")
        shell_command("mv group_* blast_groups")

    # original cutoffs: 50%, evolue -5 giant virus 50/-2 ; klebsiella 50/-5 | giant virus mai 15 50/-1
    def get_orthomcl_config(self, n, percent_match_cutoff, evalue_exponent, local=False):

        orthomcl_config = """dbVendor=mysql
dbConnectString=dbi:mysql:orthomcl;host=sib-sql04.vital-it.ch
dbLogin=orthomcl
dbPassword=orthomcl
similarSequencesTable=similar_table_%s
orthologTable=ortho_table_test_%s
orthologTable=ortholog_%s
inParalogTable=in_paralog_%s
coOrthologTable=co_ortholog_%s
interTaxonMatchView=inter_taxon_match_%s
percentMatchCutoff=%s
evalueExponentCutoff=%s
oracleIndexTblSpc=NONE""" % (n, n, n, n, n, n, percent_match_cutoff, evalue_exponent)

        orthomcl_config_local = """dbVendor=mysql
dbConnectString=dbi:mysql:orthomcl:mysql_local_infile=1
dbLogin=root
dbPassword=estrella3
similarSequencesTable=similar_table_%s
orthologTable=ortho_table_test_%s
orthologTable=ortholog_%s
inParalogTable=in_paralog_%s
coOrthologTable=co_ortholog_%s
interTaxonMatchView=inter_taxon_match_%s
percentMatchCutoff=%s
evalueExponentCutoff=%s
oracleIndexTblSpc=NONE""" % (n, n, n, n, n, n, percent_match_cutoff, evalue_exponent)

        if local:
            return orthomcl_config_local
        else:
            return orthomcl_config


            # drop_db = """
            # mysqladmin -u orthomcl -p orthomcl drop orthomcl
            # """

    def generate_config_file(self, tables_index, working_dir, percent_match_cutoff, evalue_exponent, local=False):
        config = self.get_orthomcl_config(tables_index, percent_match_cutoff, evalue_exponent, local)
        f = open(os.path.join(working_dir, "orthomcl_config.txt"), "w")
        f.write(config)
        f.close()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file", nargs='+')
    parser.add_argument("-n", '--number', type=int, help="number of sequences per file", default=1000)
    parser.add_argument("-x", '--index', type=int, help="index for orthmcl tables", default=1)
    parser.add_argument("-b", '--BLAST_file', type=str, help="BLAST file", default=False)
    parser.add_argument("-e", '--e_value_exp', type=str,
                        help="e-value exponent cutoff for BLAST and orthomcl clustering (default=-5)", default="-5")
    parser.add_argument("-m", '--percent_match_cutoff', type=str,
                        help="percent_match_cutoff for orthomcl clustering (default=50)", default="50")
    parser.add_argument("-f", '--inflation_value', type=str, help="mcl inflation value (default=1.5)", default="1.5")
    parser.add_argument("-l", '--local', action='store_true', help="run orthoml locally")

    working_dir = os.getcwd()

    args = parser.parse_args()

    test = Orthomcl(args.input_fasta, args.index, working_dir, args.BLAST_file, args.e_value_exp,
                    args.percent_match_cutoff, args.inflation_value, args.local)
# test.adjust_fasta(args.input_fasta)
