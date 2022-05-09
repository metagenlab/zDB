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

	container "$params.annotation_container"

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
	container "$params.annotation_container"

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
  
  // TODO : create a single query to the database
  process download_assembly {

	container "$params.annotation_container"

    publishDir 'data/gbk_ncbi', mode: 'copy', overwrite: true

    maxForks 2
    maxRetries 3
    //errorStrategy 'ignore'

    when:
    params.ncbi_sample_sheet

    input:
    each accession from assembly_accession_list

    cpus 1

    output:
    file '*.gbff.gz' into raw_ncbi_gbffs

    script:
    //accession = assembly_accession[0]
    """
	#!/usr/bin/env python
	import annotations
	annotations.download_assembly("$accession")
    """
  }

  // TODO: create a single query to the database
  process download_assembly_refseq {

    container "$params.annotation_container"

    publishDir 'data/refseq_corresp', mode: 'copy', overwrite: true

    maxForks 1
    maxRetries 3
    //errorStrategy 'ignore'

    echo false

    when:
    params.ncbi_sample_sheet && params.get_refseq_locus_corresp

    input:
    each accession from assembly_accession_list_refseq

    cpus 1

    output:
    file '*.gbff.gz' optional true into raw_ncbi_gbffs_refseq

    script:
    //accession = assembly_accession[0]
    """
	#!/usr/bin/env python
	import annotations
	annotations.download_assembly_refseq("$accession")
    """
  }

  // TODO : test with containers
  process refseq_locus_mapping {

	container "$params.annotation_container"

    publishDir 'data/refseq_corresp', mode: 'copy', overwrite: true

    maxForks 1
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
	import annotations
	accessions_list = "$accession_list".split(' ')
	annotations.refseq_locus_mapping(accessions_list)
    """
  }
}

all_raw_gbff = raw_ncbi_gbffs.mix(raw_local_gbffs)

process gbk_check {
  publishDir 'data/gbk_edited', mode: 'copy', overwrite: true

  container "$params.annotation_container"

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

  container "$params.annotation_container"

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
                  faa_genomes6 
                  inc_effectors_prediction }

faa_locus2.collectFile(name: 'merged.faa', newLine: true)
    .set { merged_faa0 }


process get_nr_sequences {

  container "$params.annotation_container"

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

  container "$params.orthofinder_container"

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

  container "$params.blast_container"

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

  container "$params.annotation_container"

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
  container "$params.annotation_container"

  publishDir 'orthology/orthogroups_fasta', mode: 'copy', overwrite: true

  input:
  file 'Orthogroups.txt' from orthogroups_1
  file genome_list from faa_genomes2.collect()

  output:
  file "*faa" into orthogroups_fasta

  """
	#!/usr/bin/env python
	import annotations
	annotations.orthogroups_to_fasta("$genome_list")
  """
}


process align_with_mafft {
  container "$params.mafft_container"

  publishDir 'orthology/orthogroups_alignments', mode: 'copy', overwrite: true

  input:
  file og from orthogroups_fasta.flatten().collate(20)

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
  container "$params.annotation_container"
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

  container "$params.fasttree_container"
  cpus 4
  publishDir 'orthology/orthogroups_phylogenies_fasttree', mode: 'copy', overwrite: true

  when:
  params.orthogroups_phylogeny_with_fasttree

  input:
  each file(og) from alignement_larger_than_2_seqs

  output:
    file "${og.baseName}.nwk"

  script:
  """
  FastTree ${og} > ${og.baseName}.nwk
  """
}


process get_core_orthogroups {

  container "$params.annotation_container"
  cpus 1
  memory '16 GB'
  echo false
  publishDir 'orthology/core_groups', mode: 'copy', overwrite: true

  input:
  file 'Orthogroups.txt' from orthogroups
  file genomes_list from faa_genomes3.collect()
  file fasta_files from all_alignments_3.collect()

  output:
  file '*_taxon_ids.faa' into core_orthogroups

  script:

  """
  #!/usr/bin/env python
  import annotations
  genomes_list = "$genomes_list".split()
  annotations.get_core_orthogroups(genomes_list, int("${params.core_missing}"))
  """
}

process concatenate_core_orthogroups {

  container "$params.annotation_container"

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: true

  input:
  file core_groups from core_orthogroups.collect()

  output:
  file 'msa.faa' into core_msa

  script:

  """
  #!/usr/bin/env python
  import annotations
  fasta_files = "${core_groups}".split(" ")
  annotations.concatenate_core_orthogroups(fasta_files)

  """
}

process build_core_phylogeny_with_fasttree {

  container "$params.fasttree_container"

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: true

  when:
  params.core_genome_phylogeny_with_fasttree

  input:
  file 'msa.faa' from core_msa

  output:
  file 'core_genome_phylogeny.nwk'

  script:
  '''
  FastTree -gamma -spr 4 -mlacc 2 -slownni msa.faa > core_genome_phylogeny.nwk
  '''
}


process checkm_analyse {
  container "$params.checkm_container"

  publishDir 'data/checkm/analysis', mode: 'copy', overwrite: true

  when:
  params.checkm

  input:
  file genome_list from faa_genomes5.collect()

  output:
  file "checkm_results/*" into checkm_analysis
  file "checkm_results.tab" into checkm_table

  """
  checkm analyze --genes -x faa $params.databases_dir/checkm/bacteria.ms . checkm_results -t 8 --nt
  checkm qa $params.databases_dir/checkm/bacteria.ms checkm_results -o 2 --tab_table > checkm_results.tab
  """
}


process rpsblast_COG {
  
  container "$params.annotation_container"

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

  container "$params.blast_container"

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
  ls $params.databases_dir/uniprot/swissprot/
  blastp -db $params.databases_dir/uniprot/swissprot/uniprot_sprot.fasta -query ${n} -outfmt 6 -evalue 0.001  -num_threads ${task.cpus} > ${n}.tab
  """
}

// cut sequences in smaller chunks to speed up execution?
// use a caching method to speed up queries if two identical sequences come up?
process diamond_uniref {

  publishDir 'annotation/diamond_uniref', mode: 'copy', overwrite: true

  cpus 4
  container "$params.diamond_container"

  when:
  params.diamond_uniref

  input:
  file(seq) from faa_chunks6

  output:
  file '*tab' into uniref_diamond

  script:

  n = seq.name
  """
  diamond blastp -p ${task.cpus} -d $params.databases_dir/uniref/uniref100.dmnd -q ${n} -o ${n}.tab --max-target-seqs 200 -e 0.01 --max-hsps 1
  """
}

uniref_diamond.collectFile().set { uniref_diamond_results_sqlitedb }


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

  container "$params.annotation_container"

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
	import annotations
	annotations.get_uniparc_mapping("${params.databases_dir}", "$fasta_file")
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

uniprot_mapping_tab2.splitCsv(header: true, sep: '\t')
.map{row -> "\"$row.uniprot_accession\"" }.unique().set { uniprot_nr_accessions }

process get_uniprot_data {
  container "$params.annotation_container"

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
	#!/usr/bin/env python
	import annotations
	annotations.get_uniprot_data("${params.databases_dir}", "${table}")
	"""
}



//process get_uniprot_proteome_data {

//  container "$params.annotation_container"

//  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true
//  echo true

//  when:
////  params.uniprot_data == true

//  input:
//  file "uniprot_data.tab" from uniprot_data

//  output:
//  file 'uniprot_match_annotations.db' into uniprot_db

//  script:
//  """
//  uniprot_annotations.py
//  """
//}


process get_string_mapping {

  container "$params.annotation_container"

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

    container "$params.annotation_container"

    publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

    when:
    params.string_pmid

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
    container "$params.annotation_container"

    publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

    when:
    params.string

    input:
    file 'string_mapping_PMID.tab' from string_mapping_BMID

    output:
    file 'string_mapping_PMID.db' into PMID_db

    script:
    """
    #!/usr/bin/env python
    import annotations
    annotations.get_PMID_data()
    """
}


process get_tcdb_mapping {

    container "$params.annotation_container"

    publishDir 'annotation/tcdb_mapping/', mode: 'copy', overwrite: true

    when:
    params.tcdb

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

    container "$params.chlamdb_container"

    publishDir 'annotation/tcdb_mapping', mode: 'copy', overwrite: true

    cpus 1
    maxForks 20

    echo false

    when:
    params.tcdb_gblast

    input:
    file(seq) from faa_tcdb_chunks

    output:
    file 'TCDB_RESULTS_*' into tcdb_results

    script:
    n = seq.name
    """
    export ALLOW_HTTP=true
    /usr/local/bin/BioVx/scripts/gblast3.py -i ${seq} -o TCDB_RESULTS_${seq} --db_path ${params.tcdb_db_path}
    """
}

process get_pdb_mapping {

    container "$params.annotation_container"

    publishDir 'annotation/pdb_mapping/', mode: 'copy', overwrite: true

    when:
    params.pdb

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

    container "$params.annotation_container"

    publishDir 'annotation/oma_mapping/', mode: 'copy', overwrite: true

    when:
    params.oma

    input:
    file(seq) from merged_faa5


    output:
    file 'oma_mapping.tab' into oma_mapping
    file 'no_oma_mapping.faa' into no_oma_mapping

    script:
    fasta_file = seq.name
    """
    #!/usr/bin/env python
    import annotations
    annotations.get_oma_mapping("${params.databases_dir}", "$fasta_file")
    """
}

process execute_interproscan_no_uniparc_matches {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: true
  container "$params.annotation_container"

  cpus 20
  memory '16 GB'

  when:
  params.interproscan

  input:
  file(seq) from no_uniparc_mapping_faa.splitFasta( by: 2000, file: "no_uniparc_match_chunk_" )

  output:
  file '*gff3' into interpro_gff3_no_uniparc
  file '*tsv' into interpro_tsv_no_uniparc
  file '*xml' into interpro_xml_no_uniparc
  file '*log' into interpro_log_no_uniparc

  script:
  n = seq.name
  """
  bash "$params.interproscan_home"/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3 -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus} >> ${n}.log
  """
}


process execute_interproscan_uniparc_matches {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: true
  container "$params.annotation_container"

  cpus	20
  memory '8 GB'

  when:
  params.interproscan

  input:
  file(seq) from uniparc_mapping_faa.splitFasta( by: 1000, file: "uniparc_match_chunk_" )

  output:
  file '*gff3' into interpro_gff3_uniparc
  file '*tsv' into interpro_tsv_uniparc
  file '*xml' into interpro_xml_uniparc
  file '*log' into interpro_log_uniparc

  script:
  n = seq.name
  """
  bash "$params.interproscan_home"/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3 -i ${n} -d . -T . -iprlookup -cpu ${task.cpus} >> ${n}.log
  """
}


process execute_kofamscan {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  container "$params.kegg_container"

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

  """
  exec_annotation ${n} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o ${n}.tab
  """
}


process execute_PRIAM {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  container "$params.chlamdb_container"

  cpus 2
  memory '4 GB'

  when:
  params.PRIAM

  input:
  file(seq) from faa_chunks8

  output:
  file 'results/PRIAM_*/ANNOTATION/sequenceECs.txt' into PRIAM_results

  script:
  n = seq.name
  """
  java -jar  /usr/local/bin/PRIAM_search.jar -i ${n} -o results -p $params.databases_dir/PRIAM/PRIAM_JAN18 --num_proc ${task.cpus}
  """
}


PRIAM_results.collectFile(name: 'annotation/PRIAM/sequenceECs.txt')


process setup_orthology_db {
  publishDir 'orthology/', mode: 'link', overwrite: true
  container "$params.annotation_container"

  cpus 4
  memory '8 GB'

  when:
  params.uniref_diamond_BBH_phylogeny

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


process setup_diamond_uniref_db {

  container "$params.annotation_container"
  publishDir 'annotation/diamond_uniref', mode: 'copy', overwrite: true
  echo false
  cpus 4
  memory '8 GB'

  when:
  params.uniref_diamond_BBH_phylogeny == true

  input:
  file diamond_tsv_list from uniref_diamond_results_sqlitedb.collect()

  output:
  file 'diamond_uniref.db' into diamond_uniref_db
  file 'nr_uniref_hits.tab' into uniref_diamond_nr

  script:
  """
	#!/usr/bin/env python
	import annotations
	annotations.setup_diamond_uniref_db("$diamond_tsv_list")
  """
}


process get_diamond_uniref_top_hits {

  container "$params.annotation_container"
  publishDir 'annotation/diamond_uniref_BBH_phylogenies', mode: 'copy', overwrite: true
  echo false
  cpus 4
  memory '8 GB'

  when:
  params.uniref_diamond_BBH_phylogeny

  input:
  file 'orthology.db' from orthology_db
  file 'diamond_uniref.db' from diamond_uniref_db

  output:
  file '*_nr_hits.faa' optional true into diamond_uniref_hits_fasta

  script:
  """
	#!/usr/bin/env python
	import annotations
	annotations.get_diamond_uniref_top_hits("$params.databases_dir",
	$params.uniref_diamond_BBH_phylogeny_phylum_filter,
	$params.uniref_diamond_BBH_phylogeny_top_n_hits)
  """
}

process align_uniref_BBH_with_mafft {

  container "$params.mafft_container"

  publishDir 'orthology/orthogroups_uniref_diamond_BBH_alignments', mode: 'copy', overwrite: true

  when:
  params.uniref_diamond_BBH_phylogeny == true

  input:
  file og from diamond_uniref_hits_fasta.flatten().collate( 20 )

  output:
  file "*_mafft.faa" into mafft_alignments_uniref_BBH

  script:
  """
  unset MAFFT_BINARIES
  for faa in ${og}; do
  mafft --anysymbol \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}

mafft_alignments_uniref_BBH.flatten().set{ diamond_BBH_alignments }

process orthogroup_uniref_BBH_phylogeny_with_fasttree {

  // Note : not sure fasttree is actually included in
  // annotation_container
  container "$params.fasttree_container"
  cpus 4
  publishDir 'orthology/orthogroups_uniref_diamond_BBH_phylogenies', mode: 'copy', overwrite: true

  when:
  params.uniref_diamond_BBH_phylogeny == true

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

  container "$params.annotation_container"
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

  container "$params.chlamdb_container"

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

  container "$params.chlamdb_container"

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

  container "$params.chlamdb_container"

  when:
  params.effector_prediction == true

  input:
  file "nr_fasta.faa" from nr_faa_large_sequences_chunks3

  output:
  file 'DeepT3_results.tab' into DeepT3_results

  script:
  """
  export PYTHONPATH=$PYTHONPATH:/usr/local/bin/DeepT3/DeepT3/DeepT3-Keras/
  DeepT3_scores.py -f nr_fasta.faa -o DeepT3_results.tab -d /usr/local/bin/DeepT3/DeepT3/DeepT3-Keras/
  """
}

DeepT3_results.collectFile(name: 'annotation/T3SS_effectors/DeepT3_results.tab')

process execute_T3_MM {

  container "$params.chlamdb_container"

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

  container "$params.chlamdb_container"

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

// inclusion membrane T3SS effector prediction,
// based on membrane domains as detected by interproscan
/*
process T3SS_inc_protein_prediction_interproscan_domains {
  container "$params.annotation_container"
  publishDir "annotation/T3SS_inc_effectors/", mode: "copy", overwrite: true

  input:

  output:

  when:

  script:
  """
  #!/usr/bin/env python
    
  import annotations
  
  """
}
*/

// inclusion membrane T3SS effector prediction
process execute_T3SS_inc_protein_prediction {
  container "$params.annotation_container"

  publishDir "annotation/T3SS_inc_effectors/", mode: "copy", overwrite: true

  input:
  String suffix = "_PREDICTED_INC"

  // algorithm fast enough to avoid the overhead cost of 
  // running it in separate processes
  file genomes from inc_effectors_prediction.collect()

  output:
  file "*" + suffix

  when:
  params.effector_prediction

  script:
  """
  #!/usr/bin/env python

  import annotations
  genomes_list = "${genomes}".split()
  for genome in genomes_list:
    annotations.T3SS_inc_proteins_detection(genome, genome + "${suffix}")
  """
}

// TODO see if this is interesting to make smaller chunks
// to speed up execution
process execute_psortb {
  container "$params.psort_container"
  memory '2 GB'

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
  publishDir 'annotation/pdb_mapping', mode: 'copy', overwrite: true
  container "$params.blast_container"

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
  container "$params.annotation_container"

  when:
  params.uniparc

  input:
  file(table) from uniparc_mapping_tab1

  output:
  file 'uniparc_crossreferences.tab'

  script:

  """
	#!/usr/bin/env python
	import annotations
	annotations.get_uniparc_crossreferences("${params.databases_dir}", "${table}")
  """
}


process get_idmapping_crossreferences {

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true
  container "$params.annotation_container"

  when:
  params.uniprot_idmapping

  input:
  file(table) from uniparc_mapping_tab2

  output:
  file 'idmapping_crossreferences.tab'

  script:

	"""
	#!/usr/bin/env python
	import annotations
	annotations.get_idmapping_crossreferences("${params.databases_dir}", "${table}")
	"""
}


process get_uniprot_goa_mapping {

  publishDir 'annotation/goa/', mode: 'copy', overwrite: true
  container "$params.annotation_container"
  echo true

  when:
  params.uniprot_goa

  input:
  val (uniprot_acc_list) from uniprot_nr_accessions.collect()

  output:
  file 'goa_uniprot_exact_annotations.tab'

  script:

  """
	#!/usr/bin/env python
	import annotations
	annotations.get_uniprot_goa_mapping("$params.databases_dir", $uniprot_acc_list)
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
