#!/usr/bin/env nextflow
/*
 * Author:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 *
 */


// Each Sample
if (params.ncbi_sample_sheet != false){
  Channel.fromPath( file(params.ncbi_sample_sheet) )
                      .splitCsv(header: true, sep: '\t')
                      .map{ row -> "$row.Genbank" } 
                      .into{ assembly_accession_list }

  Channel.fromPath( file(params.ncbi_sample_sheet) )
                      .splitCsv(header: true, sep: '\t')
                      .map{ row -> "$row.RefSeq" }
                      .set{ assembly_accession_list_refseq }

}

process prokka {
	publishDir 'data/prokka_output', mode: 'copy', overwrite: true

	// NOTE : according to Prokka documentation, a linear acceleration
	// is obtained up to 8 processes and after that, the overhead becomes
	// more important
	cpus 8

	when:
	params.prokka

	input:
	file genome_fasta from Channel.fromPath(params.fna_dir + "/*.fna")

	output:
	file "prokka_results/*.gbk" into gbk_from_local_assembly

	script:
	"""
	prokka $genome_fasta --outdir prokka_results --prefix ${genome_fasta.baseName} \\
		 --centre X --compliant --cpus ${task.cpus}
	"""
}

// leave all the contigs with no coding region (check_gbk doesn't like them)
process prokka_filter_CDS {
	publishDir 'data/prokka_output_filtered', mode: 'copy', overwrite: true
	input:
	file prokka_file from gbk_from_local_assembly.collect()

	when:
	params.prokka

	output:
	file "*_filtered.gbk" into gbk_prokka_filtered

	script:
	"""
	#!/usr/bin/env python
	
	import annotations

	gbk_files = "${prokka_file}".split()
	for i in gbk_files:
		annotations.filter_out_unannotated(i)
	"""
}

if(!params.prokka) {
	gbk_prokka_filtered = Channel.empty()
}


if(params.local_sample_sheet){
	local_genomes = Channel.fromPath(params.local_sample_sheet)
				  .splitCsv(header: true, sep: '\t')
				  .map{row -> "$row.gbk_path" }
				  .map { file(it) }
} else {
	local_genomes = Channel.empty()
} 

process copy_local_assemblies {
    publishDir 'data/gbk_local', mode: 'copy', overwrite: true

    cpus 1

    when:
    params.local_sample_sheet || params.prokka

	// Mixing inputs from local assemblies annotated with Prokka
	// and local annotated genomes in .gbk format
    input:
    file local_gbk from local_genomes.mix(gbk_prokka_filtered.flatMap())

    output:
    file "${local_gbk}.gz" into raw_local_gbffs

    script:
    """
    gzip -f ${local_gbk}
    """
}

if(!params.prokka && !params.local_sample_sheet) {
	raw_local_gbffs = Channel.empty()	
}

if(params.ncbi_sample_sheet == false) {
	raw_ncbi_gbffs = Channel.empty()
}else {
  process download_assembly {

    // conda 'biopython=1.73'

    publishDir 'data/gbk_ncbi', mode: 'copy', overwrite: true

    maxForks 2
    maxRetries 3
    //errorStrategy 'ignore'

    echo false

    when:
    params.ncbi_sample_sheet != false

    input:
    each accession from assembly_accession_list

    cpus 1

    output:
    file '*.gbff.gz' into raw_ncbi_gbffs

    script:
    //accession = assembly_accession[0]
    """
  #!/usr/bin/env python

  import re
  import time
  from ftplib import FTP
  from Bio import Entrez, SeqIO
  Entrez.email = "trestan.pillonel@chuv.ch"
  Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"

  print("${accession}")

  handle1 = Entrez.esearch(db="assembly", term="${accession}")
  record1 = Entrez.read(handle1)

  ncbi_id = record1['IdList'][-1]
  print(ncbi_id)
  handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
  assembly_record = Entrez.read(handle_assembly, validate=False)

  if 'genbank_has_annotation' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]["PropertyList"]:
      ftp_path = re.findall('<FtpPath type="GenBank">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
  elif 'refseq_has_annotation' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]["PropertyList"]:
      ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
  else:
    raise("${accession} assembly not annotated! --- exit ---")
  print(ftp_path)
  ftp=FTP('ftp.ncbi.nih.gov')
  ftp.login("anonymous","trestan.pillonel@unil.ch")
  ftp.cwd(ftp_path)
  filelist=ftp.nlst()
  filelist = [i for i in filelist if 'genomic.gbff.gz' in i]
  print(filelist)
  for file in filelist:
    ftp.retrbinary("RETR "+file, open(file, "wb").write)
    """
  }

  process download_assembly_refseq {

    conda 'biopython=1.73'

    publishDir 'data/refseq_corresp', mode: 'copy', overwrite: true

    maxForks 2
    maxRetries 3
    //errorStrategy 'ignore'

    echo false

    when:
    params.ncbi_sample_sheet != false && params.get_refseq_locus_corresp == true

    input:
    each accession from assembly_accession_list_refseq

    cpus 1

    output:
    file '*.gbff.gz' optional true into raw_ncbi_gbffs_refseq

    script:
    //accession = assembly_accession[0]
    """
  #!/usr/bin/env python

  import re
  import time
  from ftplib import FTP
  from Bio import Entrez, SeqIO
  Entrez.email = "trestan.pillonel@chuv.ch"
  Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"

  if "${accession}" != "":
      handle1 = Entrez.esearch(db="assembly", term="${accession}")
      record1 = Entrez.read(handle1)

      ncbi_id = record1['IdList'][-1]
      print(ncbi_id)
      handle_assembly = Entrez.esummary(db="assembly", id=ncbi_id)
      assembly_record = Entrez.read(handle_assembly, validate=False)

      # only consider annotated genbank => get corresponding refseq assembly
      if 'genbank_has_annotation' in assembly_record['DocumentSummarySet']['DocumentSummary'][0]["PropertyList"]:
          ftp_path = re.findall('<FtpPath type="RefSeq">ftp[^<]*<', assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Meta'])[0][50:-1]
          print("ftp_path", ftp_path)
          ftp=FTP('ftp.ncbi.nih.gov')
          ftp.login("anonymous","trestan.pillonel@unil.ch")
          ftp.cwd(ftp_path)
          filelist=ftp.nlst()
          filelist = [i for i in filelist if 'genomic.gbff.gz' in i]
          print(filelist)
          for file in filelist:
            ftp.retrbinary("RETR "+file, open(file, "wb").write)
      else:
        print("no genbank annotation for ${accession}")
    """
  }

  process refseq_locus_mapping {

    conda 'biopython=1.73'

    publishDir 'data/refseq_corresp', mode: 'copy', overwrite: true

    maxForks 2
    maxRetries 3
    //errorStrategy 'ignore'

    echo false

    when:
    params.get_refseq_locus_corresp == true

    input:
    file accession_list from raw_ncbi_gbffs_refseq.collect()

    cpus 1

    output:
    file 'refseq_corresp.tab'

    script:
    """
  #!/usr/bin/env python

  import re
  import time
  from Bio import SeqIO
  import gzip

  o = open("refseq_corresp.tab", "w")

  for accession in "${accession_list}".split(" "):
      print(accession)
      with gzip.open(accession, "rt") as handle:
        records = SeqIO.parse(handle, "genbank")
        for record in records:
            for feature in record.features:
                if 'protein_id' in feature.qualifiers and 'old_locus_tag' in feature.qualifiers:
                    refseq_locus_tag = feature.qualifiers["locus_tag"][0]
                    protein_id = feature.qualifiers["protein_id"][0]
                    old_locus_tag = feature.qualifiers["old_locus_tag"][0]
                    o.write(f"{old_locus_tag}\\t{refseq_locus_tag}\\t{protein_id}\\n")
    """
  }
}

all_raw_gbff = raw_ncbi_gbffs.mix(raw_local_gbffs)

process gbk_check {
  publishDir 'data/gbk_edited', mode: 'copy', overwrite: true

  // conda 'bioconda::biopython=1.68'

  cpus 2

  input:
  file all_gbff from all_raw_gbff.collect()

  output:
  file "*merged.gbk" into edited_gbks

  script:
  """
  #!/usr/bin/env python
  
  import annotations
  annotations.check_gbk("$all_gbff".split())
  """
}

process convert_gbk_to_faa {

  publishDir 'data/faa_locus', mode: 'copy', overwrite: true

  // conda 'bioconda::biopython=1.68'

  echo false

  cpus 1

  input:
  each file(edited_gbk) from edited_gbks

  output:
  file "*.faa" into faa_files

  script:
  """
	#!/usr/bin/env python
	import annotations

	annotations.convert_gbk_to_faa("${edited_gbk}", "${edited_gbk.baseName}.faa")
  """
}

faa_files.into{ faa_locus1; faa_locus2 }


faa_locus1.into { faa_genomes1
                  faa_genomes2
                  faa_genomes3
                  faa_genomes4
                  faa_genomes5
                  faa_genomes6 }

faa_locus2.collectFile(name: 'merged.faa', newLine: true)
    .set { merged_faa0 }


process get_nr_sequences {

  // conda 'bioconda::biopython=1.68'

  publishDir 'data/', mode: 'copy', overwrite: true

  input:
  file(seq) from merged_faa0
  file genome_list from faa_genomes4.collect()

  output:

  file 'nr.faa' into nr_seqs
  file 'nr_mapping.tab' into nr_mapping

  script:
  fasta_file = seq.name
  """
	#!/usr/bin/env python
	
	import annotations
	annotations.get_nr_sequences("${fasta_file}", "${genome_list}".split())
  """
}

nr_seqs.collectFile(name: 'merged_nr.faa', newLine: true)
.into { merged_faa_chunks
        merged_faa1
        merged_faa2
        merged_faa3
        merged_faa4
        merged_faa5
        merged_faa6
        merged_faa7 }

merged_faa_chunks.splitFasta( by: 1000, file: "chunk_" )
.into { faa_chunks1
        faa_chunks2
        faa_chunks3
        faa_chunks4
        faa_chunks5
        faa_chunks6
        faa_chunks7
        faa_chunks8
        faa_chunks9 }

process prepare_orthofinder {

  // conda 'bioconda::orthofinder=2.2.7'

  when:
  params.orthofinder

  input:
    file genome_list from faa_genomes1.collect()

  output:
    file "OrthoFinder/Results_$params.orthofinder_output_dir/WorkingDirectory/Species*.fa" into species_fasta
    file "OrthoFinder/Results_$params.orthofinder_output_dir/WorkingDirectory/BlastDBSpecies*.phr" into species_blastdb
    path "OrthoFinder/Results_$params.orthofinder_output_dir/" into result_dir

  script:
  """
  orthofinder -op -a 8 -n "$params.orthofinder_output_dir" -S blast -f . > of_prep.tab
  """
}

process blast_orthofinder {

  // conda 'bioconda::blast=2.7.1'

  cpus 2

  input:
  file complete_dir from result_dir
  each seq from species_fasta
  each blastdb from species_blastdb


  when:
  params.orthofinder

  output:
  file "${complete_dir.baseName}/WorkingDirectory/Blast${species_1}_${species_2}.txt" into blast_results
  

  script:
  blastdb_name = blastdb.getBaseName()
  blastdb_path = blastdb.getParent()
  seq_name = seq.getBaseName()
  species_1 =  (seq_name =~ /Species(\d+)/)[0][1]
  species_2 =  (blastdb_name =~ /BlastDBSpecies(\d+)/)[0][1]

  """
  blastp -outfmt 6 -evalue 0.001 -query $seq -db $blastdb_path/$blastdb_name -num_threads ${task.cpus} > ${complete_dir}/WorkingDirectory/Blast${species_1}_${species_2}.txt
  """
}

process orthofinder_main {

  // conda 'bioconda::orthofinder=2.2.7'

  publishDir 'orthology', mode: 'copy', overwrite: true
  echo true

  cpus 8

  input:
  file complete_dir from result_dir
  val blast_results from blast_results.collect()

  output:
  // for some reason, executing orthofinder on previous blast results makes it change directory
  // and output its results in the new directory... TODO : fixit!
  file "Results_$params.orthofinder_output_dir/WorkingDirectory/OrthoFinder/Results_$params.orthofinder_output_dir/Orthogroups/Orthogroups.txt" into orthogroups
  file "Results_$params.orthofinder_output_dir/WorkingDirectory/OrthoFinder/Results_$params.orthofinder_output_dir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" into singletons

  script:
  """
  echo "${complete_dir.baseName}"
  ls
  orthofinder -og -t ${task.cpus} -a ${task.cpus} \
	-b ./Results_$params.orthofinder_output_dir/WorkingDirectory/ > of_grouping.txt \
	-n "$params.orthofinder_output_dir"
  """
}

orthogroups
.into { orthogroups_1
        orthogroups_2}

process orthogroups2fasta {
  // conda 'bioconda::biopython=1.70'

  publishDir 'orthology/orthogroups_fasta', mode: 'copy', overwrite: true

  input:
  file 'Orthogroups.txt' from orthogroups_1
  file genome_list from faa_genomes2.collect()

  output:
  file "*faa" into orthogroups_fasta

  """
  #!/usr/bin/env python

  from Bio import SeqIO
  import os

  fasta_list = "${genome_list}".split(' ')

  sequence_data = {}
  for fasta_file in fasta_list:
      sequence_data.update(SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")))

  # write fasta
  with open("Orthogroups.txt") as f:
      all_grp = [i for i in f]
      for n, line in enumerate(all_grp):
          groups = line.rstrip().split(' ')
          group_name = groups[0][0:-1]
          groups = groups[1:len(groups)]
          if len(groups)>1:
              new_fasta = [sequence_data[i] for i in groups]
              out_path = "%s.faa" % group_name
              out_handle = open(out_path, "w")
              SeqIO.write(new_fasta, out_handle, "fasta")
  """
}


process align_with_mafft {

  // conda 'bioconda::mafft=7.407'

  publishDir 'orthology/orthogroups_alignments', mode: 'copy', overwrite: true

  input:
  file og from orthogroups_fasta.flatten().collate( 20 )

  output:
  file "*_mafft.faa" into mafft_alignments

  script:
  """
  unset MAFFT_BINARIES
  for faa in ${og}; do
  mafft \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}

/* Get only alignments with more than than two sequences */
mafft_alignments.collect().into {all_alignments_1
                                 all_alignments_2
                                 all_alignments_3
                                 all_alignments_4}

all_alignments_1.flatten().map { it }.filter { (it.text =~ /(>)/).size() > 3 }.set { alignments_larget_tah_3_seqs }
all_alignments_2.flatten().map { it }.filter { (it.text =~ /(>)/).size() == 3 }.set { alignments_3_seqs }
all_alignments_4.flatten().map { it }.filter { (it.text =~ /(>)/).size() > 2 }.set { alignement_larger_than_2_seqs }

/*
process orthogroups_phylogeny_with_raxml {

  echo false
  conda 'bioconda::raxml=8.2.9'
  cpus 4
  publishDir 'orthology/orthogroups_phylogenies', mode: 'copy', overwrite: false

  input:
  each file(og) from large_alignments

  output:
    file "${og}.nwk"

  script:
  """
  raxmlHPC -m PROTGAMMALG -p 12345 -s ${og} -n ${og.getBaseName()} -c 4 -T 4;
  raxmlHPC -f J -m PROTGAMMALG -s ${og} -p 12345 -t RAxML_result.${og.getBaseName()} -n sh_test_${og.getBaseName()} -T 4
  """
}
*/

process orthogroups_phylogeny_with_fasttree3 {

  conda 'bioconda::fasttree=2.1.10'
  cpus 4
  publishDir 'orthology/orthogroups_phylogenies_fasttree', mode: 'copy', overwrite: true

  when:
  params.orthogroups_phylogeny_with_fasttree == true

  input:
  each file(og) from alignement_larger_than_2_seqs

  output:
    file "${og.baseName}.nwk"

  script:
  """
  FastTree ${og} > ${og.baseName}.nwk
  """
}


process orthogroups_phylogeny_with_iqtree {

  conda 'bioconda::iqtree=1.6.8'
  cpus 2
  publishDir 'orthology/orthogroups_phylogenies_iqtree', mode: 'copy', overwrite: true

  when:
  params.orthogroups_phylogeny_with_iqtree == true

  input:
  each file(og) from alignments_larget_tah_3_seqs

  output:
    file "${og.getBaseName()}.iqtree"
    file "${og.getBaseName()}.treefile"
    file "${og.getBaseName()}.log"
    file "${og.getBaseName()}.bionj"
    file "${og.getBaseName()}.ckp.gz"
    file "${og.getBaseName()}.mldist"
    file "${og.getBaseName()}.model.gz"
    file "${og.getBaseName()}.splits.nex"

  script:
  """
  iqtree -nt 2 -s ${og} -alrt 1000 -bb 1000 -pre ${og.getBaseName()}
  """
}

process orthogroups_phylogeny_with_iqtree_no_boostrap {

  conda 'bioconda::iqtree=1.6.8'
  cpus 2
  publishDir 'orthology/orthogroups_phylogenies_iqtree', mode: 'copy', overwrite: true

  when:
  params.orthogroups_phylogeny_with_iqtree == true

  input:
  each file(og) from alignments_3_seqs

  output:
    file "${og.getBaseName()}.iqtree"
    file "${og.getBaseName()}.treefile"
    file "${og.getBaseName()}.log"
    file "${og.getBaseName()}.bionj"
    file "${og.getBaseName()}.ckp.gz"
    file "${og.getBaseName()}.mldist"
    file "${og.getBaseName()}.model.gz"

  script:
  """
  iqtree -nt 2 -s ${og} -pre ${og.getBaseName()}
  """
}

process get_core_orthogroups {

  // conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'
  cpus 1
  memory '16 GB'
  echo false
  publishDir 'orthology/core_groups', mode: 'copy', overwrite: true

  input:
  file 'Orthogroups.txt' from orthogroups
  file genome_list from faa_genomes3.collect()
  file fasta_files from all_alignments_3.collect()

  output:
  file '*_taxon_ids.faa' into core_orthogroups

  script:

  """
  #!/usr/bin/env python
  import annotations
  annotations.get_core_orthogroups(genome_list.split(), int("${params.core_missing}"))
  """
}

process concatenate_core_orthogroups {

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: true

  input:
  file core_groups from core_orthogroups.collect()

  output:
  file 'msa.faa' into core_msa

  script:

  """
  #!/usr/bin/env python

  fasta_files = "${core_groups}".split(" ")
  print(fasta_files)
  out_name = 'msa.faa'

  from Bio import AlignIO
  from Bio.SeqRecord import SeqRecord
  from Bio.Seq import Seq
  from Bio.Align import MultipleSeqAlignment
  # identification of all distinct fasta headers id (all unique taxons ids) in all fasta
  # storing records in all_seq_data (dico)
  taxons = []
  all_seq_data = {}
  for one_fasta in fasta_files:
      all_seq_data[one_fasta] = {}
      with open(one_fasta) as f:
          alignment = AlignIO.read(f, "fasta")
      for record in alignment:
          if record.id not in taxons:
              taxons.append(record.id)
          all_seq_data[one_fasta][record.id] = record


  # building dictionnary of the form: dico[one_fasta][one_taxon] = sequence
  concat_data = {}

  start_stop_list = []
  start = 0
  stop = 0

  for one_fasta in fasta_files:
      start = stop + 1
      stop = start + len(all_seq_data[one_fasta][list(all_seq_data[one_fasta].keys())[0]]) - 1
      start_stop_list.append([start, stop])
      print(len(taxons))
      for taxon in taxons:
          # check if the considered taxon is present in the record
          if taxon not in all_seq_data[one_fasta]:
              # if taxon absent, create SeqRecord object "-"*len(alignments): gap of the size of the alignment
              seq = Seq("-"*len(all_seq_data[one_fasta][list(all_seq_data[one_fasta].keys())[0]]))
              all_seq_data[one_fasta][taxon] = SeqRecord(seq, id=taxon)
          if taxon not in concat_data:
              concat_data[taxon] = all_seq_data[one_fasta][taxon]
          else:
              concat_data[taxon] += all_seq_data[one_fasta][taxon]

  # concatenating the alignments, writing to fasta file
  MSA = MultipleSeqAlignment([concat_data[i] for i in concat_data])
  with open(out_name, "w") as handle:
      AlignIO.write(MSA, handle, "fasta")
  """
}

process build_core_phylogeny_with_fasttree {

  conda 'bioconda::fasttree=2.1.10'

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: true

  when:
  params.core_genome_phylogeny_with_fasttree == true

  input:
  file 'msa.faa' from core_msa

  output:
  file 'core_genome_phylogeny.nwk'

  script:
  '''
  FastTree -gamma -spr 4 -mlacc 2 -slownni msa.faa > core_genome_phylogeny.nwk
  '''
}


/* process checkm_analyse {
  '''
  Get fasta file of each orthogroup
  '''

  // conda 'bioconda::checkm-genome=1.0.12'

  publishDir 'data/checkm/analysis', mode: 'copy', overwrite: true

  when:
  params.checkm == true

  input:
  file genome_list from faa_genomes5.collect()

  output:
  file "checkm_results/*" into checkm_analysis

  """
  checkm analyze --genes -x faa \$CHECKM_SET_PATH . checkm_results -t 8 --nt
  """
}


process checkm_qa {
  '''
  Get fasta file of each orthogroup
  '''

  conda 'bioconda::checkm-genome=1.0.12'

  publishDir 'data/checkm/analysis', mode: 'copy', overwrite: true

  when:
  params.checkm == true

  input:
  file checkm_analysis_results from checkm_analysis

  output:
  file "checkm_results.tab" into checkm_table

  """
  checkm qa \$CHECKM_SET_PATH . -o 2 --tab_table > checkm_results.tab
  """
} */


process rpsblast_COG {

  // conda 'bioconda::blast=2.7.1'

  cpus 4

  when:
  params.cog == true

  input:
  file 'seq' from faa_chunks1

  output:
  file 'blast_result' into blast_result

  script:
  n = seq.name
  """
  rpsblast -db $params.databases_dir/cdd/profiles/Cog -query seq -outfmt 6 -evalue 0.001 -num_threads ${task.cpus} > blast_result
  """
}

blast_result.collectFile(name: 'annotation/COG/blast_COG.tab')

process blast_swissprot {

  // conda 'bioconda::blast=2.7.1'

  publishDir 'annotation/blast_swissprot', mode: 'copy', overwrite: true

  cpus 4

  when:
  params.blast_swissprot

  input:
  file(seq) from faa_chunks3

  output:
  file '*tab' into swissprot_blast

  script:

  n = seq.name
  """
  blastp -db $params.databases_dir/uniprot/swissprot/uniprot_sprot.fasta -query ${n} -outfmt 6 -evalue 0.001  -num_threads ${task.cpus} > ${n}.tab
  """
}


process plast_refseq {

  publishDir 'annotation/plast_refseq', mode: 'copy', overwrite: true

  cpus 12

  when:
  params.plast_refseq == true

  input:
  file(seq) from faa_chunks5

  output:
  file '*tab' into refseq_plast
  file '*log' into refseq_plast_log

  script:

  n = seq.name
  """
  # 15'000'000 vs 10'000'000
  # 100'000'000 max
  # -s 45
  export PATH="\$PATH:\$PLAST_HOME"
  plast -p plastp -a ${task.cpus} -d $params.databases_dir/refseq/merged.faa.pal -i ${n} -M BLOSUM62 -s 60 -seeds-use-ratio 45 -max-database-size 50000000 -e 1e-5 -G 11 -E 1 -o ${n}.tab -F F -bargraph -verbose -force-query-order 1000 -max-hit-per-query 200 -max-hsp-per-hit 1 > ${n}.log;
  """
}

process diamond_refseq {

  publishDir 'annotation/diamond_refseq', mode: 'copy', overwrite: true

  cpus 4
  // conda 'bioconda::diamond=0.9.24'

  when:
  params.diamond_refseq

  input:
  file(seq) from faa_chunks6

  output:
  file '*tab' into refseq_diamond

  script:

  n = seq.name
  """
  diamond blastp -p ${task.cpus} -d $params.databases_dir/refseq/merged_refseq.dmnd -q ${n} -o ${n}.tab --max-target-seqs 200 -e 0.01 --max-hsps 1
  """
}

refseq_diamond.collectFile().set { refseq_diamond_results_sqlitedb }


//refseq_diamond_results_taxid_mapping
//.splitCsv(header: false, sep: '\t')
//.map{row ->
//    def protein_accession = row[1]
//    return "${protein_accession}"
//}
//.unique()
//.collectFile(name: 'nr_refseq_hits.tab', newLine: true)
//.set {refseq_diamond_nr}

//.collate( 300 )
//.set {
//    nr_refseq_hits_chunks
//}

process get_uniparc_mapping {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true

  when:
  params.uniparc == true

  input:
  file(seq) from merged_faa1

  output:
  file 'uniparc_mapping.tab' into uniparc_mapping_tab
  file 'uniprot_mapping.tab' into uniprot_mapping_tab
  file 'no_uniprot_mapping.faa' into no_uniprot_mapping_faa
  file 'no_uniparc_mapping.faa' into no_uniparc_mapping_faa
  file 'uniparc_mapping.faa' into uniparc_mapping_faa

  script:
  fasta_file = seq.name
  """
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("${params.databases_dir}/uniprot/uniparc/uniparc.db")
cursor = conn.cursor()

fasta_file = "${fasta_file}"

uniparc_map = open('uniparc_mapping.tab', 'w')
uniprot_map = open('uniprot_mapping.tab', 'w')
no_uniprot_mapping = open('no_uniprot_mapping.faa', 'w')
no_uniparc_mapping = open('no_uniparc_mapping.faa', 'w')
uniparc_mapping_faa = open('uniparc_mapping.faa', 'w')

uniparc_map.write("locus_tag\\tuniparc_id\\tuniparc_accession\\tstatus\\n")
uniprot_map.write("locus_tag\\tuniprot_accession\\ttaxon_id\\tdescription\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_mapping_uniprot_records = []
no_mapping_uniparc_records = []
mapping_uniparc_records = []

for record in records:
    match = False
    sql = 'select t1.uniparc_id,uniparc_accession,accession,taxon_id,description, db_name, status from uniparc_accession t1 inner join uniparc_cross_references t2 on t1.uniparc_id=t2.uniparc_id inner join crossref_databases t3 on t2.db_id=t3.db_id where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_mapping_uniparc_records.append(record)
        no_mapping_uniprot_records.append(record)
    else:
        mapping_uniparc_records.append(record)
        all_status = [i[6] for i in hits]
        if 1 in all_status:
            status = 'active'
        else:
            status = 'dead'
        uniparc_map.write("%s\\t%s\\t%s\\t%s\\n" % (record.id,
                                               hits[0][0],
                                               hits[0][1],
                                               status))
        for uniprot_hit in hits:
            if uniprot_hit[5] in ["UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL"] and uniprot_hit[6] == 1:
                match = True
                uniprot_map.write("%s\\t%s\\t%s\\t%s\\t%s\\n" % (record.id,
                                                                 uniprot_hit[2],
                                                                 uniprot_hit[3],
                                                                 uniprot_hit[4],
                                                                 uniprot_hit[5]))
        if not match:
            no_mapping_uniprot_records.append(record)

SeqIO.write(no_mapping_uniprot_records, no_uniprot_mapping, "fasta")
SeqIO.write(no_mapping_uniparc_records, no_uniparc_mapping, "fasta")
SeqIO.write(mapping_uniparc_records, uniparc_mapping_faa, "fasta")
  """
}

uniparc_mapping_tab.into{
  uniparc_mapping_tab1
  uniparc_mapping_tab2
}

uniprot_mapping_tab.into{
  uniprot_mapping_tab1
  uniprot_mapping_tab2
}


uniprot_mapping_tab2.splitCsv(header: false, sep: '\t')
.map{row ->
    def uniprot_accession = row[1]
    return "${uniprot_accession}"
}
.collect()
.unique()
.into {uniprot_nr_accessions}


process get_uniprot_data {

  conda 'biopython=1.73=py36h7b6447c_0'

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true
  echo true

  when:
  params.uniprot_data == true

  input:
  file(table) from uniprot_mapping_tab1

  output:
  file 'uniprot_data.tab' into uniprot_data

  script:

  """
#!/usr/bin/env python3.6

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum
import datetime
import logging

logging.basicConfig(filename='uniprot_data.log',level=logging.DEBUG)

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def uniprot_accession2score(uniprot_accession_list):

    import urllib.request
    from urllib.error import URLError

    # https://www.uniprot.org/uniprot/?query=id:V8TQN7+OR+id:V8TR74&format=xml
    link = "http://www.uniprot.org/uniprot/?query=id:%s&columns=id,annotation%%20score&format=tab" % ("+OR+id:".join(uniprot_accession_list))

    try:
        page = urllib.request.urlopen(link)
        data = page.read().decode('utf-8').split('\\n')
        rows = [i.rstrip().split('\\t') for i in data]
    except URLError:
        success = False
        while not success:
            import time
            logging.debug('connection problem, trying again...\\n')
            logging.debug('%s\\n' % link)
            time.sleep(10)
            try:
                page = urllib.request.urlopen(req)
                data = page.read().decode('utf-8').split('\\n')
                rows = [i.rstrip().split('\\t') for i in data]
                success=True
            except:
                success = False
    unirpot2score = {}
    for row in rows:
        if len(row) > 0:
            if row[0] == 'Entry':
                continue
            elif len(row)<2:
                continue
            else:
                unirpot2score[row[0]] = row[1]

    return unirpot2score


conn = sqlite3.connect("${params.databases_dir}/uniprot/idmapping/uniprot_sprot_trembl.db")
cursor = conn.cursor()

uniprot_table = open("${table}", 'r')
uniprot_data = open('uniprot_data.tab', 'w')

uniprot_data.write("uniprot_accession\\tuniprot_score\\tuniprot_status\\tproteome\\tcomment_function\\tec_number\\tcomment_subunit\\tgene\\trecommendedName_fullName\\tproteinExistence\\tdevelopmentalstage\\tcomment_similarity\\tcomment_catalyticactivity\\tcomment_pathway\\tkeywords\\n")

uniprot_accession_list = [row.rstrip().split("\\t")[1].split(".")[0] for row in uniprot_table if row.rstrip().split("\\t")[1] != 'uniprot_accession']
uniprot_accession_chunks = chunks(uniprot_accession_list, 300)

sql_uniprot_annot = 'select * from uniprot_annotation where uniprot_accession in ("%s");'

for n, one_chunk in enumerate(uniprot_accession_chunks):
    logging.debug("chunk: %s -- %s\\n" % (n, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    uniprot2score = uniprot_accession2score(one_chunk)
    filter = '","'.join(one_chunk)
    uniprot_annot_data = cursor.execute(sql_uniprot_annot % filter,).fetchall()

    for m, uniprot_annotation in enumerate(uniprot_annot_data):
        uniprot_accession = uniprot_annotation[0]
        #print("chunk:", n, "entry", m , uniprot_accession)
        try:
            uniprot_score = uniprot2score[uniprot_accession]
        # deal with eventual removed entries
        except KeyError:
            uniprot_score = 0
        comment_function = uniprot_annotation[1]
        ec_number = uniprot_annotation[2]
        comment_similarity = uniprot_annotation[3]
        comment_catalyticactivity = uniprot_annotation[4]
        comment_pathway = uniprot_annotation[5]
        keywords = uniprot_annotation[6]
        comment_subunit = uniprot_annotation[7]
        gene = uniprot_annotation[8]
        recommendedName_fullName = uniprot_annotation[9]
        proteinExistence = uniprot_annotation[10]
        developmentalstage = uniprot_annotation[11]
        proteome = uniprot_annotation[12]
        reviewed = uniprot_annotation[14]

        uniprot_data.write("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" % ( uniprot_accession,
                                                                                                             uniprot_score,
                                                                                                             reviewed,
                                                                                                             proteome,
                                                                                                             comment_function,
                                                                                                             ec_number,
                                                                                                             comment_subunit,
                                                                                                             gene,
                                                                                                             recommendedName_fullName,
                                                                                                             proteinExistence,
                                                                                                             developmentalstage,
                                                                                                             comment_similarity,
                                                                                                             comment_catalyticactivity,
                                                                                                             comment_pathway,
                                                                                                             keywords)
                                                                                                            )
    """
}


process get_uniprot_proteome_data {

  conda 'biopython=1.73=py36h7b6447c_0'

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true
  echo true

  when:
  params.uniprot_data == true

  input:
  file "uniprot_data.tab" from uniprot_data

  output:
  file 'uniprot_match_annotations.db' into uniprot_db

  script:
  """
  uniprot_annotations.py
  """
}


process get_string_mapping {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

  when:
  params.string == true

  input:
  file(seq) from merged_faa2


  output:
  file 'string_mapping.tab' into string_mapping

  script:
  fasta_file = seq.name
  """
	#!/usr/bin/env python
	import annotations

	annotations.get_string_mapping("${fasta_file}", "${params.databases_dir}")
  """
}

process get_string_PMID_mapping {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

  when:
  params.string == true

  input:
  file(string_map) from string_mapping


  output:
  file 'string_mapping_PMID.tab' into string_mapping_BMID

  script:

  """
	#!/usr/bin/env python
	import annotations
	annotations.get_string_PMID_mapping("${string_map}")
  """
}


process get_PMID_data {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

  when:
  params.string == true

  input:
  file 'string_mapping_PMID.tab' from string_mapping_BMID


  output:
  file 'string_mapping_PMID.db' into PMID_db

  script:

  """
#!/usr/bin/env python

from Bio import Entrez
import sqlite3

Entrez.email = "trestan.pillonel@chuv.ch"


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def pmid2abstract_info(pmid_list):
    from Bio import Medline

    # make sure that pmid are strings
    pmid_list = [str(i) for i in pmid_list]

    try:
        handle = Entrez.efetch(db="pubmed", id=','.join(pmid_list), rettype="medline", retmode="text")
        records = Medline.parse(handle)
    except:
        print("FAIL:", pmid)
        return None

    pmid2data = {}
    for record in records:
      try:
          pmid = record["PMID"]
      except:
          print(record)
          #{'id:': ['696885 Error occurred: PMID 28696885 is a duplicate of PMID 17633143']}
          if 'duplicate' in record['id:']:
              duplicate = record['id:'].split(' ')[0]
              correct = record['id:'].split(' ')[-1]
              print("removing duplicated PMID... %s --> %s" % (duplicate, correct))
              # remove duplicate from list
              pmid_list.remove(duplicate)
              return pmid2abstract_info(pmid_list)

      pmid2data[pmid] = {}
      pmid2data[pmid]["title"] = record.get("TI", "?")
      pmid2data[pmid]["authors"] = record.get("AU", "?")
      pmid2data[pmid]["source"] = record.get("SO", "?")
      pmid2data[pmid]["abstract"] = record.get("AB", "?")
      pmid2data[pmid]["pmid"] = pmid

    return pmid2data

conn = sqlite3.connect("string_mapping_PMID.db")
cursor = conn.cursor()

sql = 'create table if not exists hash2pmid (hash binary, ' \
      ' pmid INTEGER)'

sql2 = 'create table if not exists pmid2data (pmid INTEGER, title TEXT, authors TEXT, source TEXT, abstract TEXT)'
cursor.execute(sql,)
cursor.execute(sql2,)

# get PMID nr list and load PMID data into sqldb
pmid_nr_list = []
sql_template = 'insert into hash2pmid values (?, ?)'
with open("string_mapping_PMID.tab", "r") as f:
    n = 0
    for row in f:
        data = row.rstrip().split("\t")
        if data[1] != 'None':
            n+=1
            cursor.execute(sql_template, data)
            if data[1] not in pmid_nr_list:
                pmid_nr_list.append(data[1])
        if n % 1000 == 0:
            print(n, 'hash2pmid ---- insert ----')
            conn.commit()

pmid_chunks = [i for i in chunks(pmid_nr_list, 50)]

# get PMID data and load into sqldb
sql_template = 'insert into pmid2data values (?, ?, ?, ?, ?)'
for n, chunk in enumerate(pmid_chunks):
    print("pmid2data -- chunk %s / %s" % (n, len(pmid_chunks)))
    pmid2data = pmid2abstract_info(chunk)
    for pmid in pmid2data:
        cursor.execute(sql_template, (pmid, pmid2data[pmid]["title"], str(pmid2data[pmid]["authors"]), pmid2data[pmid]["source"], pmid2data[pmid]["abstract"]))
    if n % 10 == 0:
        conn.commit()
conn.commit()
  """
}


process get_tcdb_mapping {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/tcdb_mapping/', mode: 'copy', overwrite: true

  when:
  params.tcdb == true

  input:
  file(seq) from merged_faa3


  output:
  file 'tcdb_mapping.tab' into tcdb_mapping
  file 'no_tcdb_mapping.faa' into no_tcdb_mapping

  script:
  fasta_file = seq.name
  """
	#!/usr/bin/env python
	
	import annotations
	annotations.get_tcdb_mapping("${fasta_file}", "${params.databases_dir}")
  """
}

no_tcdb_mapping.splitFasta( by: 1000, file: "chunk" )
.set { faa_tcdb_chunks }

process tcdb_gblast3 {

  container 'metagenlab/chlamdb_annotation:1.0.3'

  publishDir 'annotation/tcdb_mapping', mode: 'copy', overwrite: true

  cpus 1
  maxForks 20

  echo false

  when:
  params.tcdb_gblast == true

  input:
  file(seq) from faa_tcdb_chunks

  output:
  file 'TCDB_RESULTS_*' into tcdb_results

  script:

  n = seq.name
  """
  # activate conda env
  source activate gblast3
  /usr/local/bin/BioVx/scripts/gblast3.py -i ${seq} -o TCDB_RESULTS_${seq} --db_path /usr/local/bin/tcdb_db/
  """
}

process get_pdb_mapping {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/pdb_mapping/', mode: 'copy', overwrite: true

  when:
  params.pdb == true

  input:
  file(seq) from merged_faa4


  output:
  file 'pdb_mapping.tab' into pdb_mapping
  file 'no_pdb_mapping.faa' into no_pdb_mapping

  script:
  fasta_file = seq.name
  """
	#!/usr/bin/env python
	import annotations
	annotations.get_pdb_mapping("${fasta_file}", "${params.databases_dir}")
  """
}

process get_oma_mapping {

  // conda 'bioconda::biopython=1.68'

  publishDir 'annotation/oma_mapping/', mode: 'copy', overwrite: true

  when:
  params.oma == true

  input:
  file(seq) from merged_faa5


  output:
  file 'oma_mapping.tab' into oma_mapping
  file 'no_oma_mapping.faa' into no_oma_mapping

  script:
  fasta_file = seq.name
  """
#!/usr/bin/env python


from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("${params.databases_dir}/oma/oma.db")
cursor = conn.cursor()

fasta_file = "${fasta_file}"

oma_map = open('oma_mapping.tab', 'w')
no_oma_mapping = open('no_oma_mapping.faa', 'w')

oma_map.write("locus_tag\\toma_id\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_oma_mapping_records = []
for record in records:
    sql = 'select accession from hash_table where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_oma_mapping_records.append(record)
    else:
        for hit in hits:
          oma_map.write("%s\\t%s\\n" % (record.id,
                                              hit[0]))


SeqIO.write(no_oma_mapping_records, no_oma_mapping, "fasta")

  """
}

process execute_interproscan_no_uniparc_matches {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: true

  cpus 8
  memory '16 GB'
  // conda 'anaconda::openjdk=8.0.152'

  when:
  params.interproscan

  input:
  file(seq) from no_uniparc_mapping_faa.splitFasta( by: 300, file: "no_uniparc_match_chunk_" )

  output:
  file '*gff3' into interpro_gff3_no_uniparc
  file '*html.tar.gz' into interpro_html_no_uniparc
  file '*svg.tar.gz' into interpro_svg_no_uniparc
  file '*tsv' into interpro_tsv_no_uniparc
  file '*xml' into interpro_xml_no_uniparc
  file '*log' into interpro_log_no_uniparc

  script:
  n = seq.name
  """
  echo $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus}
  bash $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus} >> ${n}.log
  """
}


process execute_interproscan_uniparc_matches {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: true

  cpus 8
  memory '8 GB'
  conda 'anaconda::openjdk=8.0.152'

  when:
  params.interproscan == true

  input:
  file(seq) from uniparc_mapping_faa.splitFasta( by: 1000, file: "uniparc_match_chunk_" )

  output:
  file '*gff3' into interpro_gff3_uniparc
  file '*html.tar.gz' into interpro_html_uniparc
  file '*svg.tar.gz' into interpro_svg_uniparc
  file '*tsv' into interpro_tsv_uniparc
  file '*xml' into interpro_xml_uniparc
  file '*log' into interpro_log_uniparc

  script:
  n = seq.name
  """
  echo $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . -iprlookup -cpu ${task.cpus}
  bash $INTERPRO_HOME/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . -iprlookup -cpu ${task.cpus} >> ${n}.log
  """
}


process execute_kofamscan {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  // conda 'hmmer=3.2.1 parallel ruby=2.4.5'

  cpus 4
  memory '8 GB'

  when:
  params.ko == true

  input:
  file(seq) from faa_chunks7

  output:
  file '*tab'

  script:
  n = seq.name

  // Hack for now
  """
  ${params.databases_dir}/kegg/exec_annotation ${n} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o ${n}.tab
  """
}


process execute_PRIAM {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  container 'metagenlab/chlamdb_annotation:1.0.3'

  cpus 2
  memory '4 GB'

  when:
  params.PRIAM == true

  input:
  file(seq) from faa_chunks8

  output:
  file 'results/PRIAM_*/ANNOTATION/sequenceECs.txt' into PRIAM_results

  script:
  n = seq.name
  """
  conda activate /opt/conda/envs/priam
  java -jar  /usr/local/bin/PRIAM/PRIAM_search.jar -i ${n} -o results -p $params.databases_dir/PRIAM/PRIAM_JAN18 --num_proc ${task.cpus}
  """
}


PRIAM_results.collectFile(name: 'annotation/PRIAM/sequenceECs.txt')


process setup_orthology_db {
  /*
  - [ ] load locus_tag2hash table
  - [ ] load locus2orthogroup data into sqlite db (locus2orthogroup)
  */

  // conda 'biopython=1.73=py36h7b6447c_0'
  publishDir 'orthology/', mode: 'link', overwrite: true

  cpus 4
  memory '8 GB'

  when:
  params.refseq_diamond_BBH_phylogeny

  input:
  file nr_mapping_file from nr_mapping
  file orthogroup from orthogroups_2
  file nr_fasta from merged_faa6

  output:
  file 'orthology.db' into orthology_db

  script:
  """
  #!/usr/bin/env python

  import annotations
  annotations.setup_orthology_db("${nr_fasta}", "${nr_mapping_file}", "${orthogroup}")
  """
}


process setup_diamond_refseq_db {
  /*
  - [ ] load blast results into sqlite db
  */

  // conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'
  publishDir 'annotation/diamond_refseq', mode: 'copy', overwrite: true
  echo false
  cpus 4
  memory '8 GB'

  when:
  params.refseq_diamond_BBH_phylogeny == true

  input:
  file diamond_tsv_list from refseq_diamond_results_sqlitedb.collect()

  output:
  file 'diamond_refseq.db' into diamond_refseq_db
  file 'nr_refseq_hits.tab' into refseq_diamond_nr

  script:
  """
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum
import pandas

conn = sqlite3.connect("diamond_refseq.db")
cursor = conn.cursor()

sql1 = 'create table diamond_refseq(hit_count INTEGER, qseqid varchar(200), sseqid varchar(200), pident FLOAT, length INTEGER, mismatch INTEGER, gapopen INTEGER, qstart INTEGER, qend INTEGER, sstart INTEGER, send INTEGER, evalue FLOAT, bitscore FLOAT)'
cursor.execute(sql1,)
conn.commit()

sql = 'insert into diamond_refseq values (?,?, ?,?,?,?,?,?,?,?,?,?,?)'

diamond_file_list = "${diamond_tsv_list}".split(' ')
for one_file in diamond_file_list:
    diamond_table = pandas.read_csv(one_file, sep="\\t")
    accession = ''
    count = ''
    # add hit count as first column
    for index, row in diamond_table.iterrows():
        # remove version number from accession
        row[1] = row[1].split(".")[0]
        # if new protein, reinitialise the count
        if row[0] != accession:
            accession = row[0]
            count = 1
        else:
            count+=1
        cursor.execute(sql, [count] + row.tolist())
    conn.commit()

# index query accession (hash) + hit number
sql_index_1 = 'create index hitc on diamond_refseq (hit_count);'
sql_index_2 = 'create index qacc on diamond_refseq (qseqid);'
sql_index_3 = 'create index sacc on diamond_refseq (sseqid);'

cursor.execute(sql_index_1)
cursor.execute(sql_index_2)
cursor.execute(sql_index_3)
conn.commit()
sql = 'select distinct sseqid from diamond_refseq'
with open("nr_refseq_hits.tab", 'w') as f:
    for acc in cursor.execute(sql,):
        f.write("%s\\n" % acc[0])
  """
}


process get_refseq_hits_taxonomy {

  // conda 'biopython=1.73=py36h7b6447c_0'

  publishDir 'annotation/diamond_refseq/', mode: 'copy', overwrite: true

  echo true

  when:
  params.diamond_refseq_taxonomy == true

  input:
  file refseq_hit_table from refseq_diamond_nr

  output:
  file 'refseq_taxonomy.db' into refseq_hit_taxid_mapping_db

  script:

  """
#!/usr/bin/env python

from Bio import SeqIO
import http.client
import re
import sqlite3
import datetime
import logging

logging.basicConfig(filename='refseq_taxonomy.log',level=logging.DEBUG)

conn = sqlite3.connect("refseq_taxonomy.db")
cursor = conn.cursor()

cursor.execute("PRAGMA synchronous = OFF")
cursor.execute("BEGIN TRANSACTION")


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]


def accession2taxid_entrez(accession):
    from Bio import Entrez
    Entrez.email = "trestan.pillonel@chuv.ch"
    Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"
    from socket import error as SocketError
    import errno
    try:
        handle = Entrez.esummary(db="protein", id="%s" % accession, retmax=1)
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise('error connexion with %s' % accession)
        else:
            import time
            print ('connexion error, trying again...')
            time.sleep(60)
            accession2taxid_entrez(accession)
    record = Entrez.read(handle, validate=False)[0]
    # record['AccessionVersion'],
    # record['TaxId'],
    # record['Title'],
    # record['Length']
    return int(record['TaxId'])

conn_refseq = sqlite3.connect("$params.databases_dir/refseq/merged_refseq.db")
cursor_refseq = conn_refseq.cursor()

conn_taxid = sqlite3.connect("$params.databases_dir/ncbi-taxonomy/prot_accession2taxid.db")
cursor_taxid = conn_taxid.cursor()

sql = 'create table refseq_hits (accession varchar(200), taxid INTEGER, description TEXT, length INTEGER)'
cursor.execute(sql,)
conn.commit()

f = open("${refseq_hit_table}", 'r')
# get list of accessions, remove version number
hit_list = [i.rstrip().split(".")[0] for i in f]

chunk_lists = chunks(hit_list, 5000)
sql_template = 'insert into refseq_hits values (?,?,?,?)'
list_of_lists = []
template_annotation = 'select accession, description, sequence_length from refseq where accession in ("%s")'
template_taxid = 'select accession, taxid from accession2taxid where accession in ("%s")'
for n, acc_list in enumerate(chunk_lists):
    logging.debug("chunk: %s -- %s" % (n, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    filter = '","'.join(acc_list)
    #logging.debug('SQL ------- %s' % (n, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    data_annotation = cursor_refseq.execute(template_annotation % filter,).fetchall()
    data_taxid = cursor_taxid.execute(template_taxid % filter,).fetchall()
    #logging.debug('SQL-done -- %s' % (n, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    # create dictionnary

    accession2taxid = {}
    for row in data_taxid:
        accession2taxid[row[0]] = row[1]
    for annot in data_annotation:
        annot = list(annot)
        accession = annot[0]
        seq_length = annot[2]
        # AccessionVersion, TaxId, Title, Length
        description = re.sub("%s.[0-9]+ " % annot[0], "", annot[1])
        try:
            taxid = accession2taxid[annot[0]]
        except:
            logging.debug("Missing taxid:\t%s" % accession)
            taxid = accession2taxid_entrez(accession)
        list_of_lists.append([accession, taxid, description, seq_length])
    if n % 10 == 0:
        logging.debug('insert: %s' % n)
        cursor.executemany(sql_template, list_of_lists)
        conn.commit()
        list_of_lists = []
logging.debug('LAST insert: %s' % n)
cursor.executemany(sql_template, list_of_lists)
conn.commit()
# index accession and taxid columns
sql_index_1 = 'create index acc on refseq_hits (accession);'
sql_index_2 = 'create index taxid on refseq_hits (taxid);'

cursor.execute(sql_index_1)
cursor.execute(sql_index_2)
conn.commit()

  """
}


process get_diamond_refseq_top_hits {
  /*
  - [ ] load blast results into sqlite db
  */

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'
  publishDir 'annotation/diamond_refseq_BBH_phylogenies', mode: 'copy', overwrite: true
  echo false
  cpus 4
  memory '8 GB'

  when:
  params.refseq_diamond_BBH_phylogeny == true

  input:
  file 'orthology.db' from orthology_db
  file 'refseq_taxonomy.db' from refseq_hit_taxid_mapping_db
  file 'diamond_refseq.db' from diamond_refseq_db

  output:
  file '*_nr_hits.faa' into diamond_refseq_hits_fasta

  script:
  """
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum
import pandas
from collections import Iterable

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

def refseq_accession2fasta(accession_list):
    from Bio import Entrez
    from Bio import SeqIO
    import urllib2

    Entrez.email = "trestan.pillonel@unil.ch"
    try:
        handle = Entrez.efetch(db='protein', id=','.join(accession_list), rettype="fasta", retmode="text")
    except urllib2.HTTPError as e:
        print (e)
        print (accession_list)
        print ('connexion problem, waiting 60s. and trying again...')
        # run again if connexion error
        import time
        time.sleep(60)
        return refseq_accession2fasta(accession_list)
    except urllib2.URLError:
        print (e)
        print (accession_list)
        print ('connexion problem, waiting 60s. and trying again...')
        # run again if connexion error
        import time
        time.sleep(60)
        return refseq_accession2fasta(accession_list)

    records = [i for i in SeqIO.parse(handle, "fasta")]

    return records

    from collections import Iterable


def flatten(items):
    # flatten list of lists
    merged_list = []
    for x in items:
        merged_list+=x
    return merged_list

conn = sqlite3.connect("orthology.db")
cursor = conn.cursor()
#sql1 = 'create table diamond_top_%s_nr_hits(orthogroup varchar, hit_accession, hit_sequence)' % ${params.refseq_diamond_BBH_phylogeny_top_n_hits}

sql2 = 'attach "${params.databases_dir}/ncbi-taxonomy/linear_taxonomy.db" as linear_taxonomy'
sql3 = 'attach "refseq_taxonomy.db" as refseq_taxonomy'
sql4 = 'attach "diamond_refseq.db" as diamond_refseq'

#cursor.execute(sql1,)
cursor.execute(sql2,)
cursor.execute(sql3,)
cursor.execute(sql4,)
conn.commit()

#sql_template = 'insert into diamond_top_%s_nr_hits values (?,?,?)' % ${params.refseq_diamond_BBH_phylogeny_top_n_hits}

sql_filtered_hits = 'select t1.orthogroup,t1.locus_tag,t3.sseqid,t3.hit_count from locus_tag2orthogroup t1 inner join locus_tag2sequence_hash t2 on t1.locus_tag=t2.locus_tag inner join diamond_refseq.diamond_refseq t3 on t2.sequence_hash=t3.qseqid inner join refseq_taxonomy.refseq_hits t4 on t3.sseqid=t4.accession inner join linear_taxonomy.ncbi_taxonomy t5 on t4.taxid=t5.tax_id where t5.phylum not in ("%s") order by t1.orthogroup,t1.locus_tag,t3.hit_count;'  % '","'.join(${params.refseq_diamond_BBH_phylogeny_phylum_filter})
print(len(sql_filtered_hits))
filtered_hits = cursor.execute(sql_filtered_hits,).fetchall()

# retrieve top hits excluding phylum of choice
orthogroup2locus2top_hits = {}
for row in filtered_hits:
    orthogroup = row[0]
    locus_tag = row[1]
    hit_id = row[2]
    if orthogroup not in orthogroup2locus2top_hits:
        orthogroup2locus2top_hits[orthogroup] = {}
    if locus_tag not in orthogroup2locus2top_hits[orthogroup]:
        orthogroup2locus2top_hits[orthogroup][locus_tag] = []
    if len(orthogroup2locus2top_hits[orthogroup][locus_tag]) < int(${params.refseq_diamond_BBH_phylogeny_top_n_hits}):
        orthogroup2locus2top_hits[orthogroup][locus_tag].append(hit_id)

# retrieve aa sequences
sql = 'select orthogroup,t1.locus_tag,sequence from locus_tag2orthogroup t1 inner join locus_tag2sequence_hash t2 on t1.locus_tag=t2.locus_tag inner join sequence_hash2aa_sequence t3 on t2.sequence_hash=t3.sequence_hash '
orthogroup_sequences = cursor.execute(sql,).fetchall()
orthogroup2locus_and_sequence = {}
for row in orthogroup_sequences:
    orthogroup = row[0]
    locus_tag = row[1]
    sequence = row[2]
    if orthogroup not in orthogroup2locus_and_sequence:
        orthogroup2locus_and_sequence[orthogroup] = []
    orthogroup2locus_and_sequence[orthogroup].append([locus_tag, sequence])


# for each group, retrieve aa sequence from NCBI
# write it to fasta file with orthogroup sequences
for group in orthogroup2locus2top_hits:
    ortho_sequences = orthogroup2locus_and_sequence[group]
    refseq_sequence_records = []
    top_hits = [orthogroup2locus2top_hits[group][i] for i in orthogroup2locus2top_hits[group]]
    nr_top_hit = list(set(flatten(top_hits)))
    split_lists = chunks(nr_top_hit, 50)
    for one_list in split_lists:
        refseq_sequence_records += refseq_accession2fasta(one_list)

    if len(refseq_sequence_records) > 0:
        with open("%s_nr_hits.faa" % group, 'w') as f:
            for one_ortholog in ortho_sequences:
                f.write(">%s\\n%s\\n" % (one_ortholog[0], one_ortholog[1]))
            for record in refseq_sequence_records:
                name = record.name
                f.write(">%s\\n%s\\n" % (name, str(record.seq)))
  """
}

process align_refseq_BBH_with_mafft {

  conda 'bioconda::mafft=7.407'

  publishDir 'orthology/orthogroups_refseq_diamond_BBH_alignments', mode: 'copy', overwrite: true

  when:
  params.refseq_diamond_BBH_phylogeny == true

  input:
  file og from diamond_refseq_hits_fasta.flatten().collate( 20 )

  output:
  file "*_mafft.faa" into mafft_alignments_refseq_BBH

  script:
  """
  unset MAFFT_BINARIES
  for faa in ${og}; do
  mafft --anysymbol \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}

mafft_alignments_refseq_BBH.flatten().set{ diamond_BBH_alignments }

process orthogroup_refseq_BBH_phylogeny_with_fasttree {

  conda 'bioconda::fasttree=2.1.10'
  cpus 4
  publishDir 'orthology/orthogroups_refseq_diamond_BBH_phylogenies', mode: 'copy', overwrite: true

  when:
  params.refseq_diamond_BBH_phylogeny == true

  input:
  each file(og) from diamond_BBH_alignments

  output:
    file "${og.baseName}.nwk"

  script:
  """
  FastTree ${og} > ${og.baseName}.nwk
  """
}


// Filter out small sequences and ambiguous AA
process filter_sequences {
  // conda 'biopython=1.73'

  publishDir 'data/', mode: 'copy', overwrite: true

  when:
  params.effector_prediction

  input: 
    file nr_fasta from merged_faa7

  output:
    file "filtered_sequences.faa" into filtered_sequences
  
  script:
	"""
	#!/usr/bin/env python
	import annotations
	annotations.filter_sequences("${nr_fasta}")
	"""
}

filtered_sequences.splitFasta( by: 5000, file: "chunk_" ).into{ nr_faa_large_sequences_chunks1
                                                                    nr_faa_large_sequences_chunks2
                                                                    nr_faa_large_sequences_chunks3
                                                                    nr_faa_large_sequences_chunks4 }

process execute_BPBAac {

  container 'metagenlab/chlamdb_annotation:1.0.3'

  when:
  params.effector_prediction == true

  input:
  file "nr_fasta.faa" from nr_faa_large_sequences_chunks1

  output:
  file 'BPBAac_results.csv' into BPBAac_results

  script:
  """
  ln -s /usr/local/bin/BPBAac/BPBAacPre.R
  ln -s /usr/local/bin/BPBAac/BPBAac.Rdata
  ln -s /usr/local/bin/BPBAac/Classify.pl
  ln -s /usr/local/bin/BPBAac/DataPrepare.pl
  ln -s /usr/local/bin/BPBAac/Feature100.pl
  ln -s /usr/local/bin/BPBAac/Integ.pl
  ln -s /usr/local/bin/BPBAac/Neg100AacFrequency
  ln -s /usr/local/bin/BPBAac/Pos100AacFrequency
  sed 's/CRC-/CRC/g' nr_fasta.faa > edited_fasta.faa
  perl DataPrepare.pl edited_fasta.faa;
  Rscript BPBAacPre.R;
  sed -i 's/CRC/CRC-/g' FinalResult.csv;
  mv FinalResult.csv BPBAac_results.csv;
  """
}

BPBAac_results.collectFile(name: 'annotation/T3SS_effectors/BPBAac_results.tab')

process execute_effectiveT3 {

  container 'metagenlab/chlamdb_annotation:1.0.3'

  when:
  params.effector_prediction == true

  input:
  file "nr_fasta.faa" from nr_faa_large_sequences_chunks2

  output:
  file 'effective_t3.out' into effectiveT3_results

  script:
  """
  ln -s /usr/local/bin/effective/module
  java -jar /usr/local/bin/effective/TTSS_GUI-1.0.1.jar -f nr_fasta.faa -m TTSS_STD-2.0.2.jar -t cutoff=0.9999 -o effective_t3.out -q
  """
}

effectiveT3_results.collectFile(name: 'annotation/T3SS_effectors/effectiveT3_results.tab')

process execute_DeepT3 {

  container 'metagenlab/chlamdb_annotation:1.0.3'

  when:
  params.effector_prediction == true && false

  input:
  file "nr_fasta.faa" from nr_faa_large_sequences_chunks3

  output:
  file 'DeepT3_results.tab' into DeepT3_results

  script:
  """
  DeepT3_scores.py -f nr_fasta.faa -o DeepT3_results.tab -d /usr/local/bin/DeepT3/DeepT3/DeepT3-Keras/
  """
}

DeepT3_results.collectFile(name: 'annotation/T3SS_effectors/DeepT3_results.tab')

process execute_T3_MM {

  container 'metagenlab/chlamdb_annotation:1.0.3'

  when:
  params.effector_prediction == true

  input:
  file "nr_fasta.faa" from nr_faa_large_sequences_chunks4

  output:
  file 'T3_MM_results.csv' into T3_MM_results

  script:
  """
  ln -s /usr/local/bin/T3_MM/log.freq.ratio.txt
  perl /usr/local/bin/T3_MM/T3_MM.pl nr_fasta.faa
  mv final.result.csv T3_MM_results.csv
  """
}


T3_MM_results.collectFile(name: 'annotation/T3SS_effectors/T3_MM_results.tab')


process execute_macsyfinder_T3SS {

  container 'metagenlab/chlamdb_annotation:1.0.3'

  publishDir 'annotation/macsyfinder/T3SS/', mode: 'copy', overwrite: true

  when:
  params.macsyfinder == true

  input:
  file(genome) from faa_genomes6

  output:
  file "macsyfinder-${genome.baseName}/macsyfinder.conf"
  file "macsyfinder-${genome.baseName}/macsyfinder.log"
  file "macsyfinder-${genome.baseName}/macsyfinder.out"
  file "macsyfinder-${genome.baseName}/macsyfinder.report"
  file "macsyfinder-${genome.baseName}/macsyfinder.summary"

  script:
  """
  python /usr/local/bin/macsyfinder/bin/macsyfinder --db-type unordered_replicon --sequence-db ${genome} -d /usr/local/bin/annotation_pipeline_nextflow/data/macsyfinder/secretion_systems/TXSS/definitions -p /usr/local/bin/annotation_pipeline_nextflow/data/macsyfinder/secretion_systems/TXSS/profiles all --coverage-profile 0.2 --i-evalue-select 1 --e-value-search 1
  mv macsyfinder-* macsyfinder-${genome.baseName}
  """
}


process execute_psortb {

  container 'metagenlab/psort:3.0.6'

  publishDir 'annotation/psortb/', mode: 'copy', overwrite: true

  when:
  params.psortb == true

  input:
  file(chunk) from faa_chunks9

  output:
  file "psortb_${chunk}.txt" into psortb_results

  script:
  """
  /usr/local/psortb/bin/psort --negative ${chunk} > psortb_${chunk}.txt
  """
}


process blast_pdb {

  // conda 'bioconda::blast=2.7.1'

  publishDir 'annotation/pdb_mapping', mode: 'copy', overwrite: true

  cpus 4

  when:
  params.pdb == true

  input:
  file(seq) from no_pdb_mapping.splitFasta( by: 1000, file: "chunk_" )

  output:
  file '*tab' into pdb_blast

  script:

  n = seq.name
  """
  blastp -db $params.databases_dir/pdb/pdb_seqres.faa -query ${n} -outfmt 6 -evalue 0.001  -num_threads ${task.cpus} > ${n}.tab
  """
}


process get_uniparc_crossreferences {

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true

  when:
  params.uniparc == true

  input:
  file(table) from uniparc_mapping_tab1

  output:
  file 'uniparc_crossreferences.tab'

  script:

  """
#!/usr/bin/env python3.6
import sqlite3
import datetime
import logging

logging.basicConfig(filename='uniparc_crossrefs.log',level=logging.DEBUG)


conn = sqlite3.connect("${params.databases_dir}/uniprot/uniparc/uniparc.db")
cursor = conn.cursor()

o = open("uniparc_crossreferences.tab", "w")

sql = 'select db_name,accession, status from uniparc_cross_references t1 inner join crossref_databases t2 on t1.db_id=t2.db_id where uniparc_id=? and db_name not in ("SEED", "PATRIC", "EPO", "JPO", "KIPO", "USPTO");'

with open("${table}", 'r') as f:
    for n, row in enumerate(f):
        if n == 0:
            continue
        data = row.rstrip().split("\t")
        uniparc_id = str(data[1])
        checksum = data[0]
        cursor.execute(sql, [uniparc_id])
        crossref_list = cursor.fetchall()
        for crossref in crossref_list:
            db_name = crossref[0]
            db_accession = crossref[1]
            entry_status = crossref[2]
            o.write("%s\\t%s\\t%s\\t%s\\n" % ( checksum,
                                          db_name,
                                          db_accession,
                                          entry_status))

    """
}


process get_idmapping_crossreferences {

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true
  echo true

  when:
  params.uniprot_idmapping

  input:
  file(table) from uniparc_mapping_tab2

  output:
  file 'idmapping_crossreferences.tab'

  script:

  """
#!/usr/bin/env python3.6
import sqlite3
import datetime
import logging

logging.basicConfig(filename='uniparc_crossrefs.log',level=logging.DEBUG)


conn = sqlite3.connect("${params.databases_dir}/uniprot/idmapping/idmapping.db")
cursor = conn.cursor()

o = open("idmapping_crossreferences.tab", "w")

sql = 'select uniprokb_accession, db_name,accession from uniparc2uniprotkb t1 inner join uniprotkb_cross_references t2 on t1.uniprotkb_id=t2.uniprotkb_id inner join database t3 on t2.db_id=t3.db_id where uniparc_accession=?;'

with open("${table}", 'r') as f:
    for n, row in enumerate(f):
        if n == 0:
            continue
        data = row.rstrip().split("\t")
        uniparc_accession = data[2]
        checksum = data[0]
        cursor.execute(sql, [uniparc_accession])
        crossref_list = cursor.fetchall()
        for crossref in crossref_list:
            uniprot_accession = crossref[0]
            db_name = crossref[1]
            db_accession = crossref[2]
            o.write("%s\\t%s\\t%s\\t%s\\n" % ( checksum,
                                               uniprot_accession,
                                               db_name,
                                               db_accession))

    """
}


process get_uniprot_goa_mapping {

  publishDir 'annotation/goa/', mode: 'copy', overwrite: true
  echo true

  when:
  params.uniprot_goa == true

  input:
  val (uniprot_acc_list) from uniprot_nr_accessions

  output:
  file 'goa_uniprot_exact_annotations.tab'

  script:

  """
#!/usr/bin/env python3.6
import sqlite3
import datetime

conn = sqlite3.connect("${params.databases_dir}/goa/goa_uniprot.db")
cursor = conn.cursor()

o = open("goa_uniprot_exact_annotations.tab", "w")

sql = 'select GO_id, reference,evidence_code,category from goa_table where uniprotkb_accession=?;'

for accession in "${uniprot_acc_list}".split(","):
    accession_base = accession.split(".")[0].strip()
    cursor.execute(sql, [accession_base])
    #print(f'select GO_id, reference,evidence_code,category from goa_table where uniprotkb_accession="{accession_base}";')
    go_list = cursor.fetchall()
    for go in go_list:
        GO_id = go[0]
        reference = go[1]
        evidence_code = go[2]
        category = go[3]
        o.write(f"{accession_base}\\t{GO_id}\\t{reference}\\t{evidence_code}\\t{category}\\n")
    """
}


workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  mail = [ to: 'trestan.pillonel@gmail.com',
           subject: 'Annotation Pipeline - DONE',
           body: 'SUCCESS!' ]
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
