#!/usr/bin/env nextflow

/*
 * Authors:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 * - Alessia Carrara
 * - Bastian Marquis
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



// Input processing
Channel.fromPath(params.local_assemblies)
    .splitCsv(header: true, strip: true)
    .map { row -> file(row.file) }
    .into { gbk_from_local_assembly_f; error_search }

gbk_from_local_assembly_f.filter { it.extension == "gbk" }
    .set { gbk_from_local_assembly }

error_search.filter { it.extension!="gbk" }
    .subscribe { error "Unsupported file extension" }


process check_gbk {
    container "$params.annotation_container"

	input:
	    file gbk from gbk_from_local_assembly.collect()

	output:
	    file "*_filtered.gbk" into checked_gbks

	script:
	"""
        #!/usr/bin/env python
        import annotations

        gbk_files = "${gbk}".split() 
        annotations.check_gbk(gbk_files)
	"""
}


checked_gbks.into {
    to_load_gbk_into_db
    to_convert_gbk }
 

process convert_gbk {
  container "$params.annotation_container"
  input:
      each file(edited_gbk) from to_convert_gbk

  output:
      file "*.faa" into faa_files
      file "*.ffn" into ffn_files_seq
      file "*.fna" into fna_files_SEQ

  script:
  """
    #!/usr/bin/env python
    import annotations

    annotations.convert_gbk_to_fasta("${edited_gbk}", "${edited_gbk.baseName}.faa")
    annotations.convert_gbk_to_fna("${edited_gbk}", "${edited_gbk.baseName}.ffn")
    annotations.convert_gbk_to_fasta("${edited_gbk}", "${edited_gbk.baseName}.fna",
        output_fmt="fna", keep_pseudo=True)
  """
}


faa_files.into{ faa_locus1; faa_locus2; faa_files_SEQ }
faa_locus2.collectFile(name: 'merged.faa', newLine: true)
    .set { to_get_nr_sequences }

faa_locus1.into { faa_genomes1
                  faa_genomes2
                  faa_genomes3
                  faa_genomes4
                  to_checkm }

process get_nr_sequences {
  container "$params.annotation_container"

  input:
    file(seq) from to_get_nr_sequences
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


nr_seqs.splitFasta( by: 300, file: "chunk_" )
.into { to_rpsblast_COG
        to_blast_swissprot
        to_diamond_refseq
        to_kofamscan
        to_pfam_scan }

if(params.pfam_scan) {
    process pfam_scan {
        container "$params.pfam_scan_container"

        input:
            file faa_chunk from to_pfam_scan

        output:
            file pfam_result_file into pfam_results

        script:
        pfam_result_file="${faa_chunk}_results"
        """
            pfam_scan.pl -f $faa_chunk -d $params.pfam_db > $pfam_result_file
        """
    }
} else {
    Channel.value("void").set { pfam_results }
}


fna_files_SEQ.into {fna_files_SEQ_1; fna_files_SEQ_2}
faa_files_SEQ.into {faa_files_SEQ_1; faa_files_SEQ_2}
ffn_files_seq.into {ffn_files_seq_1; ffn_files_seq_2}

faa_files_SEQ_1.collectFile(name: 'merged.faa', newLine: true).set { merged_faa_makeblastdb }
fna_files_SEQ_1.collectFile(name: "merged.fna", newLine: true).set { merged_fna_makeblastdb }
ffn_files_seq_1.collectFile(name: 'merged.ffn', newLine: true).set { merged_ffn_makeblastdb }

fna_files_SEQ_2.mix(faa_files_SEQ_2, ffn_files_seq_2, merged_ffn_makeblastdb,
                    merged_fna_makeblastdb, merged_faa_makeblastdb).set { to_makeblastdb }


process makeblastdb {
    container "$params.blast_container"
    publishDir "${params.results_dir}/blast_DB/$workflow.runName/${file_type}"

    input:
        file(input_file) from to_makeblastdb

    output:
        file "${input_file.baseName}*"
        file input_file

    script:
    file_type="${input_file.extension}"
    dbtype=(file_type=="ffn" || file_type=="fna")?"nucl":"prot"
    """
    makeblastdb -in ${input_file} -input_type fasta -parse_seqids \
        -dbtype $dbtype -out ${input_file.baseName}
    """
}


process prepare_orthofinder {
  container "$params.orthofinder_container"

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

  input:
  file complete_dir from result_dir
  each seq from species_fasta
  each blastdb from species_blastdb

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
  container "$params.orthofinder_container"

  cpus 2

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
  publishDir "${params.results_dir}/alignments/$workflow.runName"

  input:
      file og from orthogroups_fasta.toSortedList().flatten().collate(20)

  output:
      file "*_mafft.faa" into mafft_alignments

  script:
  """
  unset MAFFT_BINARIES
  for faa in ${og}; do
      mafft --anysymbol \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}


mafft_alignments.into { to_identity_calculation; to_phylogeny; to_cleanup }

to_cleanup.collect().set { to_alignment_gather }
to_phylogeny.collect().into { all_alignments_3; all_alignments_4 }


to_identity_calculation.toSortedList().flatten().collate(50).
    set { to_identity_calculation_split }


process identity_calculation {
    container "$params.annotation_container"
    
    input:
        file input_fasta from to_identity_calculation_split

    output:
        file "*${suffix}" into to_load_alignment

    script:
    suffix = "_ident.csv"
    """
    #!/usr/bin/env python
    import annotations
    
    input_fasta = "${input_fasta}".split()
    for file in input_fasta:
        output_filename = file+"${suffix}"
        annotations.calculate_og_identities(file, output_filename)
    """
}


all_alignments_4.flatten().map { it }.filter { (it.text =~ /(>)/).size() > 2 }.set { alignement_larger_than_2_seqs }


alignement_larger_than_2_seqs.toSortedList().flatten().collate(50).set { to_fasttree_orthogroups }


process orthogroups_phylogeny_with_fasttree3 {
  container "$params.fasttree_container"
  publishDir "${params.results_dir}/gene_phylogenies/$workflow.runName"

  input:
    file og from to_fasttree_orthogroups

  output:
    file "*.nwk" into gene_phylogeny

  script:
  """
    for i in ${og}; do
        FastTree \$i > \${i/.faa/.nwk}
    done
  """
}


gene_phylogeny.collect().set { all_og_phylogeny }


process get_core_orthogroups {
    container "$params.annotation_container"

  input:
  file 'Orthogroups.txt' from orthogroups
  file genomes_list from faa_genomes3.collect()
  file fasta_files from all_alignments_3.collect()

  output:
  file '*_taxon_ids.faa' into core_orthogroups
  file "orthogroups_summary_info.tsv" into og_summary

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
  publishDir "${params.results_dir}/gene_phylogenies/$workflow.runName"

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

  input:
  file genome_list from to_checkm.collect()

  output:
    file "checkm_results.tab" into checkm_table
  """
  checkm taxonomy_wf ${params.checkm_args} -x faa . checkm_results --tab_table  \
        --genes -f checkm_results.tab
  """
}


if(params.cog) {
    process rpsblast_COG {
      container "$params.blast_container"


      input:
      file seq from to_rpsblast_COG

      output:
      file result_file into COG_to_load_db

      script:
      n = seq.name
      result_file = "${n}.tab"
      """
      rpsblast -db $params.cog_db/cog_db -query $seq -outfmt 6 -evalue 0.001 \
                -num_threads ${task.cpus} > ${result_file}
      """
    }
} else {
    Channel.value("void").set { COG_to_load_db }
}


if (params.blast_swissprot) {
    process blast_swissprot {

      container "$params.blast_container"

      input:
      file(seq) from to_blast_swissprot

      output:
      file '*tab' into swissprot_blast

      script:

      n = seq.name
      """
      blastp -db $params.swissprot_db/swissprot.fasta -query ${n} \
                -outfmt 6 -evalue 0.001 > ${n}.tab
      """
    }
} else {
    Channel.value("void").set { swissprot_blast }
}

if(params.diamond_refseq) {
    process diamond_refseq {
      container "$params.diamond_container"

      input:
      file(seq) from to_diamond_refseq

      output:
      file '*tab' into refseq_diamond

      script:

      n = seq.name
      """
      # new version of the database
      diamond blastp -p ${task.cpus} -d $params.refseq_db/merged_refseq.dmnd \
            -q ${n} -o ${n}.tab --max-target-seqs 200 -e 0.01 --max-hsps 1
      """
    }

    refseq_diamond.collectFile().set { refseq_diamond_results_sqlitedb }
}


if(params.ko) {
    process execute_kofamscan {
      container "$params.kegg_container"

      input:
      file(seq) from to_kofamscan

      output:
      file '*tab' into to_load_KO

      script:
      n = seq.name

      """
      exec_annotation ${n} -p ${params.ko_db}/profiles/prokaryote.hal \
            -k ${params.ko_db}/ko_list --cpu ${task.cpus} -o ${n}.tab
      """
    }
} else {
    Channel.value("void").set { to_load_KO }
}

Channel.fromPath("${params.zdb.file}").set { db_skeleton }


process setup_db {
    container "$params.annotation_container"
    input:
        file db_skeleton

    output:
        file output_file into db_base

    script:
    output_file = "$workflow.runName"
    """
    cp $db_skeleton $output_file
    """
}

process load_base_db {
    container "$params.annotation_container"
    publishDir "${params.results_dir}/db"

    // Necessary to prevent segfaults due to the large size
    // of the stage script
    beforeScript "ulimit -Ss unlimited"

    input:
        file db_base
		file gbks from to_load_gbk_into_db
        file orthofinder from to_load_orthofinder_in_db
        file alignments from to_load_alignment.collect()
        file nr_mapping_file from nr_mapping_to_db_setup
        file checkm_results from checkm_table
        file core_phylogeny from core_genome_phylogeny
        file og_phylogenies_list from all_og_phylogeny
        file og_summary

    output:
        file db_name into db_gen

    script:
    db_name="$db_base"
    """
    #!/usr/bin/env python

    import setup_chlamdb
    
    kwargs = ${gen_python_args()}
    gbk_list = "${gbks}".split()
    alignments_lst = "$alignments".split()

    print("Loading gbks", flush=True)
    setup_chlamdb.load_gbk(gbk_list, kwargs, "$db_base")

    # kept for now, need to check whether this is really necessary to keep the hash
    print("Loading seq hashes", flush=True)
    setup_chlamdb.load_seq_hashes(kwargs, "$nr_mapping_file", "$db_base")

    print("Loading orthofinder results", flush=True)
    setup_chlamdb.load_orthofinder_results("$orthofinder", kwargs, "$db_base")

    print("Loading alignments", flush=True)
    setup_chlamdb.load_alignments_results(kwargs, alignments_lst, "$db_base")

    print("Loading checkm results", flush=True)
    setup_chlamdb.load_genomes_info(kwargs, gbk_list, "$checkm_results", "$db_base")

    print("Loading phylogenies")
    gene_list = "$og_phylogenies_list".split(" ")

    setup_chlamdb.load_reference_phylogeny(kwargs, "$core_phylogeny", "$db_base")
    setup_chlamdb.load_gene_phylogenies(kwargs, "$og_summary", gene_list, "$db_base")
    """
}


if(!params.diamond_refseq) {
    db_gen.set { to_load_BBH_phylo }
    Channel.empty().set { diamond_best_hits }
} else {
    process load_refseq_results {
        container "$params.annotation_container"
        input:
            file diamond_tsv_list from refseq_diamond_results_sqlitedb.collect()
            file curr_db from db_gen

        output:
            file "*_nr_hits.faa" into diamond_best_hits
            file curr_db into to_load_BBH_phylo

        script:
        """
        #!/usr/bin/env python
        # small modif'
        import setup_chlamdb
        
        kwargs = ${gen_python_args()}
        diamond_tab_files = "$diamond_tsv_list".split()
        setup_chlamdb.load_refseq_matches_infos(kwargs, diamond_tab_files, "$curr_db")
        """
    }
}


process align_refseq_BBH_with_mafft {
  container "$params.mafft_container"

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


process orthogroup_refseq_BBH_phylogeny_with_fasttree {
  container "$params.fasttree_container"

  input:
    file og from mafft_alignments_refseq_BBH

  output:
    file "*.nwk" into BBH_phylogenies

  script:
  """
    for orthogroup in $og;
    do
        FastTree \$orthogroup > \${orthogroup/.faa/.nwk}
    done
  """
}

BBH_phylogenies.collect().set { BBH_phylogenies_to_db }


if(!params.diamond_refseq) {
    Channel.value("dummy").set { BBH_phylogenies_to_db }
}


process load_BBH_phylogenies {
    container "$params.annotation_container"

    input:
        file db from to_load_BBH_phylo
        file BBH_phylogeny_trees from BBH_phylogenies_to_db

    output:
        file db into to_load_COG

    script:

    if(params.diamond_refseq)
    """
    #!/usr/bin/env python

    import setup_chlamdb

    kwargs = ${gen_python_args()}

    BBH_list = "$BBH_phylogeny_trees".split(" ")
    setup_chlamdb.load_BBH_phylogenies(kwargs, BBH_list, "$db")
    """
    else
    """
    echo "Nothing to load, continuing"
    """
}


process load_COG_into_db {
    container "$params.annotation_container"
    input:
        file db from to_load_COG
        file cog_file from COG_to_load_db.collect()

    output:
        file db into to_load_KO_db
    
    script:
    if(params.cog)
        """
        #!/usr/bin/env python

        import setup_chlamdb
        
        kwargs = ${gen_python_args()}
        cog_files = "${cog_file}".split()
        setup_chlamdb.load_cog(kwargs, cog_files, "$db", \
            "${params.cog_db}/cdd_to_cog")
        """
    else
        """
        echo \"Not supposed to load COG, passing\"
        """
}


process load_KO_into_db {
    container "$params.annotation_container"
    input:
        file KO_results from to_load_KO.collect()
        file db from to_load_KO_db

    output:
        file db into to_pfam_db

    script:
    if(params.ko)
        """
        #!/usr/bin/env python

        kwargs = ${gen_python_args()}
        ko_files = "${KO_results}".split()
        import setup_chlamdb

        # this last function should be exported in a separate script to generate
        # the scaffold of a database
        setup_chlamdb.load_KO(kwargs, ko_files, "$db")
        setup_chlamdb.load_module_completeness(kwargs, "$db")
        """
    else
        """
        echo \"Not supposed to load kegg orthologs, passing\"
        """
}

process load_PFAM_info_db {
    container "$params.annotation_container"
    input:
        file db from to_pfam_db
        file pfam_annot from pfam_results.collect()

    output:
        file db into to_load_swissprot_hits

    script:
    if(params.pfam_scan)
        """
#!/usr/bin/env python
import setup_chlamdb

kwargs = ${gen_python_args()}
pfam_files = "$pfam_annot".split()

setup_chlamdb.load_pfam(kwargs, pfam_files, "$db", "$params.pfam_db/Pfam-A.hmm.dat")
        """
    else
        """
        echo \"Not supposed to load pfam, passing\"
        """
}


process load_swissprot_hits_into_db {
    container "$params.annotation_container"
    input:
        file db from to_load_swissprot_hits
        file blast_results from swissprot_blast.collect()

    output:
        file db into to_create_index

    script:
    if(params.blast_swissprot)
    """
        #!/usr/bin/env python
        import setup_chlamdb

        kwargs = ${gen_python_args()}
        blast_results = "$blast_results".split()

        setup_chlamdb.load_swissprot(kwargs, blast_results, "$db", "$params.swissprot_db/swissprot.fasta")
    """
    else
    """
        echo \"Not supposed to load swissprot hits, passing\"
    """
}


process create_chlamdb_search_index {
    container "$params.annotation_container"
    publishDir "${params.results_dir}/search_index/"

    input:
        file db from to_create_index

    output:
        file index_name into to_index_cleanup
        file db into to_db_cleanup

    script:
    // add a prefix to differentiate it from the db
    index_name = "index_$workflow.runName"
    """
        #!/usr/bin/env python
        import setup_chlamdb

        params = ${gen_python_args()}
        setup_chlamdb.setup_chlamdb_search_index(params, "$db", "$index_name")
    """
}


process cleanup {
    input:
        file index from to_index_cleanup
        file db from to_db_cleanup
        file alignments from to_alignment_gather

    script:
    custom_run_name=(params.name)?params.name:""
    results_dir="$baseDir/${params.results_dir}"
    """
    ln -sf $index ${results_dir}/search_index/$workflow.runName

    if [ ! -d "${results_dir}/.completed_runs" ]; then
        mkdir ${results_dir}/.completed_runs
    fi

    if [ ! -z "${custom_run_name}" ]; then
        echo $workflow.runName > ${results_dir}/.completed_runs/$custom_run_name
    fi
    echo $workflow.runName > ${results_dir}/.completed_runs/$workflow.runName
    echo $workflow.runName > ${results_dir}/.completed_runs/latest
    """
}


workflow.onComplete {
    println "Annotation pipeline completed"
    println "You may now start the web server by running the start_webserver script"
}
