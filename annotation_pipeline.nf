#!/usr/bin/env nextflow
/*
 * Author:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 *
 */

/*
* Helper function for gen_python_args
*/
def gen_arg_string_s(String name, String element) {
    return "\"" + name + "\" : " + element
}

def gen_arg_string_m(String name, Map m) {
    acc = []
    for (el in m) {
        acc += gen_arg_string(name + "." + el.key, el.value)
    }
    return acc
}

// really ugly, but since nextflow does not allow
// multiple function with the same name but different signature, 
// this is the only solution
def gen_arg_string(String s, Object o) {
    prefix = "\"" + s + "\" : "
    if(o instanceof String){
        return prefix + "\"" + (String)o + "\""
    } else if(o instanceof Boolean) {
        Boolean b = (Boolean)o
        return prefix + (b ? "True" : "False")
    } else if(o instanceof Integer){
        Integer i = (Integer)o
        return prefix + i
    } else if(o instanceof Map) {
        return gen_arg_string_m(s, o)
    } else if(o instanceof List) {
        // for now, only list of string are supported
        acc = []
        for (el in (List)o) {
            acc += "\"" + (String)el + "\""
        }
        return prefix + "[" + acc.join(", ") + "]"
    }else {
        print("wrong type for " + s + " " + o)
        assert(false)
    }
}

/*
* This function generates a string that can be interpreted as map by python, from params.
* The idea is to be able to pass the parameters to external python scripts
*/
def gen_python_args() {
    acc = []
    for (parameter in params) {
        acc += gen_arg_string(parameter.key, parameter.value)
    }
    return "{" + acc.join(", ") + "}"
}

str_pythonized_params = gen_python_args()

if (params.local_assemblies) {
    Channel.fromPath(params.local_assemblies_tsv)
        .splitCsv(header: true, sep: '\t')
        .map { row -> file(params.local_assemblies_link_dir + "/" + row.draft_genome) }
        .set { to_prokka_local_assemblies }
} else {
    // declare an empty channel to avoid any complaints from nextflow
    Channel.empty().set { to_prokka_local_assemblies }
}

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
	container "$params.prokka_container"

	// NOTE : according to Prokka documentation, a linear acceleration
	// is obtained up to 8 processes and after that, the overhead becomes
	// more important. 
    cpus 4

	input:
	    file genome_fasta from to_prokka_local_assemblies

	output:
	    file "$output_dir/${output_file_prefix}.gbk" into gbk_from_local_assembly

	script:
    output_dir = "prokka_results"
    output_file_prefix = "${genome_fasta.baseName}"
	"""
	prokka $genome_fasta --outdir $output_dir --prefix ${output_file_prefix} \\
		 --centre X --compliant --cpus ${task.cpus}
	"""
}

// leave all the contigs with no coding region (check_gbk doesn't like them)
process prokka_filter_CDS {
	publishDir 'data/prokka_output_filtered', overwrite: true
	container "$params.annotation_container"

	input:
	    file prokka_file from gbk_from_local_assembly.collect()

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

    when:
    params.local_sample_sheet

    input:
        file local_gbk from local_genomes

    output:
        file "${local_gbk}.gz" into raw_local_gbffs

    script:
    """
    gzip -f ${local_gbk}
    """
}

if(!params.local_assemblies && !params.local_sample_sheet) {
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

    publishDir 'data/refseq_corresp', mode: 'copy', overwrite: true

    maxForks 2
    maxRetries 3
    //errorStrategy 'ignore'

    echo false

    when:
    params.ncbi_sample_sheet && params.get_refseq_locus_corresp

    input:
    each accession from assembly_accession_list_refseq

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

    maxForks 2
    maxRetries 3
    //errorStrategy 'ignore'

    echo false

    when:
    params.get_refseq_locus_corresp == true

    input:
    file accession_list from raw_ncbi_gbffs_refseq.collect()

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

all_raw_gbff = raw_ncbi_gbffs.mix(raw_local_gbffs).mix(gbk_prokka_filtered)

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

edited_gbks.into { to_load_gbk_into_db; to_convert_gbk_to_faa }

process convert_gbk_to_faa {

  publishDir 'data/faa_locus', mode: 'copy', overwrite: true

  container "$params.annotation_container"

  echo false

  input:
  each file(edited_gbk) from to_convert_gbk_to_faa

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
                  to_checkm
                  to_macsysfinder 
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
  file 'nr_mapping.tab' into nr_mapping_to_db_setup

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
        to_uniparc_mapping
        to_string_mapping
        to_tcdb_mapping
        to_pdb_mapping
        to_oma_mapping
        to_db_setup
        to_filter_sequences }

merged_faa_chunks.splitFasta( by: 300, file: "chunk_" )
.into { to_rpsblast_COG
        to_blast_swissprot
        to_plat_refseq
        to_diamond_refseq
        to_kofamscan
        to_PRIAM
        to_psortb }

process prepare_orthofinder {

  container "$params.annotation_container"

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

  container "$params.annotation_container"

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
        orthogroups_2
        to_load_orthofinder_in_db }

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
  container "$params.annotation_container"

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
mafft_alignments.collect().into { all_alignments_1
                                 all_alignments_2
                                 all_alignments_3
                                 all_alignments_4
                                 to_load_alignment }

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


// ToDo: group alignments to avoid using a process per alignment
process orthogroups_phylogeny_with_fasttree3 {

  container "$params.annotation_container"
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


process orthogroups_phylogeny_with_iqtree {

  container "$params.iqtree_container"
  cpus 2
  publishDir 'orthology/orthogroups_phylogenies_iqtree', mode: 'copy', overwrite: true

  when:
  params.orthogroups_phylogeny_with_iqtree

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
  iqtree -nt ${task.cpus} -s ${og} -alrt 1000 -bb 1000 -pre ${og.getBaseName()}
  """
}

process orthogroups_phylogeny_with_iqtree_no_boostrap {

  container "$params.iqtree_container"
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
  iqtree -nt ${task.cpus} -s ${og} -pre ${og.getBaseName()}
  """
}

process get_core_orthogroups {

  container "$params.annotation_container"
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

  # bugfix
  import annotations
  fasta_files = "${core_groups}".split(" ")
  annotations.concatenate_core_orthogroups(fasta_files)

  """
}

process build_core_phylogeny_with_fasttree {

  container "$params.annotation_container"

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: true

  when:
  params.core_genome_phylogeny_with_fasttree

  input:
  file 'msa.faa' from core_msa

  output:
    file 'core_genome_phylogeny.nwk' into core_genome_phylogeny

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
  file genome_list from to_checkm.collect()

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

  when:
  params.cog == true

  input:
  file 'seq' from to_rpsblast_COG

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

  when:
  params.blast_swissprot

  input:
  file(seq) from to_blast_swissprot

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

  when:
  params.plast_refseq == true

  input:
  file(seq) from to_plat_refseq

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

  container "$params.annotation_container"

  when:
  params.diamond_refseq

  input:
  file(seq) from to_diamond_refseq

  output:
  file '*tab' into refseq_diamond

  script:

  n = seq.name
  """
  # new version of the database
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

  container "$params.annotation_container"

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true

  when:
  params.uniparc == true

  input:
  file(seq) from to_uniparc_mapping

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
	#!/usr/bin/env python3.6
	import annotations
	annotations.get_uniprot_data("${params.databases_dir}", "${table}")
	"""
}


process get_uniprot_proteome_data {

  container "$params.annotation_container"

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

  container "$params.annotation_container"

  publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

  when:
  params.string == true

  input:
  file(seq) from to_string_mapping


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

  container "$params.annotation_container"

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
	import annotations
	annotations.get_PMID_data()

  """
}


process get_tcdb_mapping {

  container "$params.annotation_container"

  publishDir 'annotation/tcdb_mapping/', mode: 'copy', overwrite: true

  when:
  params.tcdb == true

  input:
  file(seq) from to_tcdb_mapping


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
  /usr/local/bin/BioVx/scripts/gblast3.py -i ${seq} -o TCDB_RESULTS_${seq} --db_path /usr/local/bin/tcdb_db/
  """
}

process get_pdb_mapping {

  container "$params.annotation_container"

  publishDir 'annotation/pdb_mapping/', mode: 'copy', overwrite: true

  when:
  params.pdb == true

  input:
  file(seq) from to_pdb_mapping


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
  params.oma == true

  input:
  file(seq) from to_oma_mapping


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

  when:
  params.interproscan

  input:
  file(seq) from no_uniparc_mapping_faa.splitFasta( by: 200, file: "no_uniparc_match_chunk_" )

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
  bash "$params.interproscan_home"/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . --disable-precalc -cpu ${task.cpus} >> ${n}.log
  """
}


process execute_interproscan_uniparc_matches {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: true
  container "$params.annotation_container"

  when:
  params.interproscan

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
  bash "$params.interproscan_home"/interproscan.sh --pathways --enable-tsv-residue-annot -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . -iprlookup -cpu ${task.cpus} >> ${n}.log
  """
}


process execute_kofamscan {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  container "$params.kegg_container"

  when:
  params.ko == true

  input:
  file(seq) from to_kofamscan

  output:
  file '*tab'

  script:
  n = seq.name

  """
  exec_annotation ${n} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list.txt --cpu ${task.cpus} -o ${n}.tab
  """
}


process execute_PRIAM {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  container "$params.chlamdb_container"

  when:
  params.PRIAM

  input:
  file(seq) from to_PRIAM

  output:
  file 'results/PRIAM_*/ANNOTATION/sequenceECs.txt' into PRIAM_results

  script:
  n = seq.name
  """
  java -jar  /usr/local/bin/PRIAM_search.jar -i ${n} -o results -p $params.databases_dir/PRIAM/PRIAM_JAN18 --num_proc ${task.cpus}
  """
}


PRIAM_results.collectFile(name: 'annotation/PRIAM/sequenceECs.txt')


// create a database containing the following tables:
// - sequences hash -> sequence aa
// - locus tag -> sequence hash
// - locus tag -> orthogroup
/* process setup_orthology_db {
  publishDir 'orthology/', overwrite: true
  container "$params.annotation_container"

  when:
  params.refseq_diamond_BBH_phylogeny

  input:
  file nr_mapping_file from nr_mapping
  file orthogroup from orthogroups_2
  file nr_fasta from to_setup_orthology_db

  output:
  file 'orthology.db' into orthology_db

  script:
  """
  #!/usr/bin/env python

  import annotations
  annotations.setup_orthology_db("${nr_fasta}", "${nr_mapping_file}", "${orthogroup}")
  """
}
 */

// transfers the diamond matches to a sqlite database
// NOTE: may be worth to only store the needed fields and limit
// the number of hits to what is needed
/*
process setup_diamond_refseq_db {

  container "$params.annotation_container"
  publishDir 'annotation/diamond_refseq', mode: 'copy', overwrite: true

  when:
  params.refseq_diamond_BBH_phylogeny

  input:
  file diamond_tsv_list from refseq_diamond_results_sqlitedb.collect()

  output:
  file 'diamond_refseq.db' into diamond_refseq_db
  file 'nr_refseq_hits.tab' into refseq_diamond_nr

  script:
  """
	#!/usr/bin/env python
    
    # bugfix
	import annotations
	annotations.setup_diamond_refseq_db("$diamond_tsv_list")
  """
}
*/


// Get the taxids of all the accessions that were matched by diamond
// insert them in refseq_taxonomy.db
// Table is :
//  accession, taxid, description and length
// 
// Note: this may be slow as it relies on db queries + web queries for 
// the accession number that couldn't be resolved locally

/*
process get_refseq_hits_taxonomy {

  publishDir 'annotation/diamond_refseq/', mode: 'copy', overwrite: true
  container "$params.annotation_container"

  when:
  params.diamond_refseq_taxonomy

  input:
  file refseq_hit_table from refseq_diamond_nr

  output:
  file 'refseq_taxonomy.db' into refseq_hit_taxid_mapping_db

  script:
  """
	#!/usr/bin/env python
    # test

	import annotations
	annotations.get_refseq_hits_taxonomy("$refseq_hit_table", "$params.databases_dir")
  """
}
*/


/*
process get_diamond_refseq_top_hits {

  container "$params.annotation_container"
  publishDir 'annotation/diamond_refseq_BBH_phylogenies', mode: 'copy', overwrite: true

  when:
  params.refseq_diamond_BBH_phylogeny

  input:
  file 'orthology.db' from orthology_db
  file 'refseq_taxonomy.db' from refseq_hit_taxid_mapping_db
  file 'diamond_refseq.db' from diamond_refseq_db

  output:
  file '*_nr_hits.faa' optional true into diamond_refseq_hits_fasta

  script:
  """
    #!/usr/bin/env python
    import annotations
    kwargs = $str_pythonized_params
    annotations.get_diamond_refseq_top_hits(kwargs)
  """
}
*/



// Filter out small sequences and ambiguous AA
process filter_sequences {

  container "$params.annotation_container"
  publishDir 'data/', mode: 'copy', overwrite: true

  when:
  params.effector_prediction

  input: 
    file nr_fasta from to_filter_sequences

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
  file(genome) from to_macsysfinder

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

  publishDir 'annotation/psortb/', mode: 'copy', overwrite: true

  when:
  params.psortb == true

  input:
  file(chunk) from to_psortb

  output:
  file "psortb_${chunk}.txt" into psortb_results

  script:
  """
  /usr/local/psortb/bin/psort --negative ${chunk} > psortb_${chunk}.txt
  """
}

process blast_pdb {
  publishDir 'annotation/pdb_mapping', mode: 'copy', overwrite: true
  container "$params.annotation_container"

  when:
  params.pdb

  input:
  file(seq) from no_pdb_mapping.splitFasta( by: 500, file: "chunk_" )

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
	annotations.get_uniparc_crossreferences("${params.databases_dir}" , "${table}")
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

process create_db {
    publishDir "db"

    input:
		file gbks from to_load_gbk_into_db
        file orthofinder from to_load_orthofinder_in_db
        file alignments from to_load_alignment
        file diamond_tsv_list from refseq_diamond_results_sqlitedb.collect()
        file nr_mapping_file from nr_mapping_to_db_setup
        file nr_fasta from to_db_setup

    output:
        file db_name into db_gen

    when:
        params.chlamdb_setup

    script:
    db_name="$params.chlamdb.db_name"
    """
    #!/usr/bin/env python

    import setup_chlamdb
    
    kwargs = ${gen_python_args()}
    gbk_list = "${gbks}".split()
    alignments_lst = "$alignments".split()
    diamond_tab_files = "$diamond_tsv_list".split()

    setup_chlamdb.setup_chlamdb(**kwargs)
    print("Loading gbks")
    setup_chlamdb.load_gbk(gbk_list, kwargs)

    # kept for now, need to check whether this is really necessary to keep the hash
    print("Loading seq hashes")
    setup_chlamdb.load_seq_hashes(kwargs, "$nr_mapping_file", "$nr_fasta")

    print("Loading orthofinder results")
    setup_chlamdb.load_orthofinder_results("$orthofinder", kwargs)

    print("Loading alignments")
    setup_chlamdb.load_alignments_results(kwargs, alignments_lst)

    print("Loading refseq matches")
    hsh_sseqid = setup_chlamdb.load_refseq_matches(kwargs, diamond_tab_files)

    print("Loading refseq matches infos")
    setup_chlamdb.load_refseq_matches_infos(kwargs, hsh_sseqid)

    print("Loading refseq matches taxonomy")
    setup_chlamdb.load_refseq_matches_linear_taxonomy(kwargs)
    """
}

db_gen.into { to_BBH; db_taxo }

process extract_non_PVC_best_hits_sequences {
    input:
        file curr_db from to_BBH

    when:
        params.refseq_diamond_BBH_phylogeny

    output:
        file "*_nr_hits.faa" into diamond_best_hits

    script:
    """
    #!/usr/bin/env python
    # final version?

    import annotations

    kwargs = ${gen_python_args()}
    annotations.get_diamond_top_hits(kwargs)
    """
}

process align_refseq_BBH_with_mafft {
  container "$params.annotation_container"

  publishDir 'orthology/orthogroups_refseq_diamond_BBH_alignments', mode: 'copy', overwrite: true

  input:
    file og from diamond_best_hits.flatten().collate( 20 )

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
  container "$params.annotation_container"
  publishDir 'orthology/orthogroups_refseq_diamond_BBH_phylogenies', mode: 'copy', overwrite: true

  input:
  each file(og) from diamond_BBH_alignments

  output:
    file "${og.baseName}.nwk" into BBH_phylogenies

  script:
  """
  FastTree ${og} > ${og.baseName}.nwk
  """
}

BBH_phylogenies.collect().into {BBH_phylogenies_to_db}

process load_taxo_stats_into_db {

    input:
        file db from db_taxo
        file BBH_phylogeny_trees from BBH_phylogenies_to_db
        file core_phylogeny from core_genome_phylogeny

    output:
        file db

    script:
    """
    #!/usr/bin/env python

    import setup_chlamdb

    kwargs = ${gen_python_args()}

    BBH_list = "$BBH_phylogeny_trees".split(" ")
    gene_list = "$gene_phylogeny".split(" ")

    setup_chlamdb.load_reference_phylogeny(kwargs, "$core_phylogeny")
    setup_chlamdb.load_gene_phylogenies(kwargs, gene_list)
    setup_chlamdb.load_BBH_phylogenies(kwargs, BBH_list)
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
