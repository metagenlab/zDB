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
	def __init__(self, input_fasta, index_mysql_tables, working_dir, BLAST_file):


		if not BLAST_file:
			self.generate_config_file(index_mysql_tables, working_dir)
			print "install shchema"
			ortho_install_shema = "orthomclInstallSchema orthomcl_config.txt install_schema.log"
	       		
			sterr, stout, process = shell_command(ortho_install_shema)
			print sterr
			print stout
			print "Adjust fasta"
			self.adjust_fasta(input_fasta)
			ortho_filter = "orthomclFilterFasta fasta 10 20"
			stdout_str, stderr_str, process = shell_command(ortho_filter)
			print "Running BLAST"
			self.run_blast("goodProteins.fasta")
			print "orthomclBlastParser"
			ortho_blast_parser = "orthomclBlastParser merged_blast.txt fasta >> similarSequences.txt"
			script = generate_bsub_file.BSUB_script(command=ortho_blast_parser, mem_in_GB=2, name="format", log_file="out_BlastParser.txt", error_file="err_BlastParser.txt")
			generate_bsub_file.run_job_and_wait(script)


		else:
			print "Parsing %s BLAST file" % BLAST_file
			ortho_blast_parser = "orthomclBlastParser %s fasta >> similarSequences.txt" % (BLAST_file)
			script = generate_bsub_file.BSUB_script(command=ortho_blast_parser, mem_in_GB=2, name="format", log_file="out_BlastParser.txt", error_file="err_BlastParser.txt")

		print "orthomclLoadBlast "
		ortho_load_blast = "orthomclLoadBlast orthomcl_config.txt similarSequences.txt"
		script = generate_bsub_file.BSUB_script(command=ortho_load_blast, mem_in_GB=2, name="format", log_file="out_LoadBlast.txt", error_file="err_LoadBlast.txt")
		generate_bsub_file.run_job_and_wait(script)

		print "orthomclPairs"
		ortho_pairs = "orthomclPairs orthomcl_config.txt orthomcl_pairs.log cleanup=yes"
		script = generate_bsub_file.BSUB_script(command=ortho_pairs, mem_in_GB=2, name="format", log_file="out_Pairs.txt", error_file="err_Pairs.txt")
		generate_bsub_file.run_job_and_wait(script)
		

		print "orthomclDumpPairsFiles"
		ortho_dump = "orthomclDumpPairsFiles orthomcl_config.txt"

		script = generate_bsub_file.BSUB_script(command=ortho_dump, mem_in_GB=2, name="format", log_file="out_DumpPairsFiles.txt", error_file="err_DumpPairsFiles.txt")
		generate_bsub_file.run_job_and_wait(script)

		print "mcl"
		ortho_mcl = "mcl mclInput --abc -I 1.5 -o mclOutput"
		script = generate_bsub_file.BSUB_script(command=ortho_mcl, mem_in_GB=2, name="format", log_file="out_mcl.txt", error_file="err_mcl.txt")
		generate_bsub_file.run_job_and_wait(script)

		print "drop database"
		cmd = "drop.sh orthomcl orthomcl orthomcl"
		script = generate_bsub_file.BSUB_script(command=ortho_mcl, mem_in_GB=2, name="format", log_file="out_drop.txt", error_file="err_drop.txt")
		generate_bsub_file.run_job_and_wait(script)

		print "groups"
		stdout_str, stderr_str, process = shell_command("orthomclMclToGroups group_ 0 < mclOutput > mclOutput.grp")

		print "singletons"
		stdout_str, stderr_str, process = shell_command("orthomclSingletons goodProteins.fasta mclOutput.grp > singletons.txt")
		 
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
			count+=1
			#outname = re.sub(".faa", "", one_file)
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
			script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=2, name="format",  log_file="out_AdjustFasta.txt", error_file="err_AdjustFasta.txt")
			job_ids.append(generate_bsub_file.run_job(script))
		print job_ids
		generate_bsub_file.wait_multi_jobs(job_ids)
		if not os.path.exists("fasta"):
			os.mkdir("fasta")
		shell_command("mv *fasta fasta")
				
				
	def run_blast(self, fasta_file):
		n_files = splitFasta.split_fasta(fasta_file, 100)
		formatdb = "formatdb -i %s  -p T" % fasta_file
		script = generate_bsub_file.BSUB_script(command=formatdb, mem_in_GB=2, name="format", log_file="out_formatdb.txt", error_file="err_formatdb.txt")
		generate_bsub_file.run_job_and_wait(script)
		job_ids = []
		for i in range(1, n_files+1):
			cmd = "blastp -db %s -query group_%s.fasta -out blast_result_%s.txt -outfmt 6 -num_threads 8 -evalue 0.00001" % (fasta_file, i, i)
			script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=4, name="blast", log_file="out_blast.txt", error_file="err_blast.txt")
			job_ids.append(generate_bsub_file.run_job(script))
		generate_bsub_file.wait_multi_jobs(job_ids)
		cmd = "for i in {1..%s}; do echo $i; cat blast_result_$i.txt>>merged_blast.txt;done" % (i)
		script = generate_bsub_file.BSUB_script(command=cmd, mem_in_GB=4, name="concat",  log_file="out_cat_fasta.txt", error_file="err_cat_fasta.txt")
		generate_bsub_file.run_job_and_wait(script)
		if not os.path.exists("blast_groups"):
			os.mkdir("blast_groups")
		shell_command("mv blast_result_* blast_groups")
		shell_command("mv group_* blast_groups")

	def get_orthomcl_config(self, n):
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
percentMatchCutoff=50
evalueExponentCutoff=-5
oracleIndexTblSpc=NONE""" % (n, n, n, n, n, n)

	    return orthomcl_config


    #drop_db = """
    #mysqladmin -u orthomcl -p orthomcl drop orthomcl
    #"""

	def generate_config_file(self, tables_index, working_dir):
		config = self.get_orthomcl_config(tables_index)
		f = open(os.path.join(working_dir, "orthomcl_config.txt"), "w")
		f.write(config)
		f.close()






if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta',type=str,help="input fasta file", nargs='+')
    parser.add_argument("-n", '--number',type=int,help="number of sequences per file", default = 1000)
    parser.add_argument("-x", '--index',type=int,help="index for orthmcl tables", default=1)
    parser.add_argument("-b", '--BLAST_file',type=int,help="BLAST file", default=False)

    working_dir = os.getcwd()
    
    args = parser.parse_args()
    
    test = Orthomcl(args.input_fasta, args.index, working_dir, args.BLAST_file)
    #test.adjust_fasta(args.input_fasta)
