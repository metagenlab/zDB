#!/usr/bin/env nextflow

/*
 * Authors:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 * - Alessia Carrara
 * - Bastian Marquis
 * - Niklaus Johner
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
    if(o instanceof GString){
        o = o.toString()
    }
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

process check_reference_databases {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path ref_db

    output:
        val 'dummy'

    script:
    """
        #!/usr/bin/env python
        from annotations import check_reference_databases

        check_reference_databases(${gen_python_args()})
    """
}


process check_gbk {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        val dummy
        path gbk
        path input_file

    output:
        path "filtered/*"
        path "filtered"

    script:
    """
        #!/usr/bin/env python
        import annotations
        import os

        if not os.path.isdir("filtered"):
          os.mkdir("filtered")

        annotations.InputHandler("$input_file").check_and_revise_gbks()
    """
}

process convert_gbk {
  container "$params.annotation_container"
  conda "$baseDir/conda/annotation.yaml"

  input:
      path edited_gbk

  output:
      path "*.faa"
      path "*.ffn"
      path "*.fna"

  script:
  """
    #!/usr/bin/env python
    import annotations

    annotations.convert_gbk_to_fasta("${edited_gbk}", "${edited_gbk.baseName}.faa")
    annotations.convert_gbk_to_fna("${edited_gbk}", "${edited_gbk.baseName}.fna")
    annotations.convert_gbk_to_fasta("${edited_gbk}", "${edited_gbk.baseName}.ffn",
        output_fmt="fna", keep_pseudo=True)
  """
}


process get_nr_sequences {
  container "$params.annotation_container"
  conda "$baseDir/conda/annotation.yaml"

  input:
    path(seq)
    path genome_list

  output:

  path 'nr.faa'
  path 'nr_mapping.tab'

  script:
  fasta_file = seq.name
  """
    #!/usr/bin/env python

    import annotations
    annotations.get_nr_sequences("${fasta_file}", "${genome_list}".split())
  """
}


process pfam_scan {
    container "$params.pfam_scan_container"
    conda "$baseDir/conda/pfam_scan.yaml"

    input:
        tuple (path(pfam_db), path(faa_chunk) )

    output:
        path pfam_result_file

    script:
    pfam_result_file="${faa_chunk}_results"
    """
        pfam_scan.pl -f $faa_chunk -d pfam > $pfam_result_file
    """
}

process makeblastdb {
    container "$params.blast_container"
    conda "$baseDir/conda/blast.yaml"
    publishDir "${params.results_dir}/blast_DB/$workflow.runName/${file_type}", mode: 'copy'

    input:
        path(input_file)

    output:
        path "${input_file.baseName}*"
        path input_file

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
  conda "$baseDir/conda/orthofinder.yaml"

  input:
    path genome_list

  output:
    path "OrthoFinder/Results_$params.orthofinder_output_dir/WorkingDirectory/Species*.fa"
    path "OrthoFinder/Results_$params.orthofinder_output_dir/WorkingDirectory/BlastDBSpecies*.phr"
    path "OrthoFinder/Results_$params.orthofinder_output_dir/"

  script:
  """
  orthofinder -op -a 8 -n "$params.orthofinder_output_dir" -S blast -f . > of_prep.tab
  """
}

process blast_orthofinder {
  container "$params.orthofinder_container"
  conda "$baseDir/conda/orthofinder.yaml"

  input:
  path complete_dir
  each seq
  each blastdb

  output:
  path "${complete_dir.baseName}/WorkingDirectory/Blast${species_1}_${species_2}.txt"


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
  conda "$baseDir/conda/orthofinder.yaml"

  cpus 2

  input:
  path complete_dir
  val blast_results

  output:
  // for some reason, executing orthofinder on previous blast results makes it change directory
  // and output its results in the new directory... TODO : fixit!
  path "Results_$params.orthofinder_output_dir/WorkingDirectory/OrthoFinder/Results_$params.orthofinder_output_dir/Orthogroups/Orthogroups.txt"
  path "Results_$params.orthofinder_output_dir/WorkingDirectory/OrthoFinder/Results_$params.orthofinder_output_dir/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

  script:
  """
  orthofinder -og -t ${task.cpus} -a ${task.cpus} \
    -b ./Results_$params.orthofinder_output_dir/WorkingDirectory/ > of_grouping.txt \
    -n "$params.orthofinder_output_dir"
  """
}

process orthogroups2fasta {
  container "$params.annotation_container"
  conda "$baseDir/conda/annotation.yaml"

  input:
  path 'Orthogroups.txt'
  path genome_list

  output:
  path "*faa"

  """
    #!/usr/bin/env python
    import annotations
    annotations.orthogroups_to_fasta("$genome_list")
  """
}

process align_with_mafft {
  container "$params.mafft_container"
  conda "$baseDir/conda/mafft.yaml"
  publishDir "${params.results_dir}/alignments/$workflow.runName", mode: "copy"

  input:
      path og

  output:
      path "*_mafft.faa"

  script:
  """
  unset MAFFT_BINARIES
  for faa in ${og}; do
      mafft --anysymbol \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}

process identity_calculation {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path input_fasta

    output:
        path "*${suffix}"

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

process orthogroups_phylogeny_with_fasttree3 {
  container "$params.fasttree_container"
  conda "$baseDir/conda/fasttree.yaml"
  publishDir "${params.results_dir}/gene_phylogenies/$workflow.runName"

  input:
    path og

  output:
    path "*.nwk"

  script:
  """
    for i in ${og}; do
        FastTree \$i > \${i/.faa/.nwk}
    done
  """
}

process get_core_orthogroups {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

  input:
  path 'Orthogroups.txt'
  path genomes_list
  path fasta_files

  output:
  path '*_taxon_ids.faa'
  path "orthogroups_summary_info.tsv"

  script:

  """
    #!/usr/bin/env python
    import annotations
    genomes_list = "$genomes_list".split()
    annotations.get_core_orthogroups(genomes_list, int("${params.core_missing}"))
  """
}

// TODO: merge with get_core_orthogroups
process concatenate_core_orthogroups {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

  input:
  path core_groups

  output:
  path 'msa.faa'

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
  conda "$baseDir/conda/fasttree.yaml"
  publishDir "${params.results_dir}/gene_phylogenies/$workflow.runName"

  input:
  path 'msa.faa'

  output:
    path 'core_genome_phylogeny.nwk'

  script:
  '''
  FastTree -gamma -spr 4 -mlacc 2 -slownni msa.faa > core_genome_phylogeny.nwk
  '''
}

process checkm_analyse {
  container "$params.checkm_container"
  conda "$baseDir/conda/checkm.yaml"

  input:
  path genome_list

  output:
    path "checkm_results.tab"
  """
  checkm taxonomy_wf ${params.checkm_args} -x faa . checkm_results --tab_table  \
        --genes -f checkm_results.tab
  """
}

process rpsblast_COG {
  container "$params.blast_container"
  conda "$baseDir/conda/blast.yaml"

  input:
    tuple (file(cog_db), file(seq))

  output:
  path result_file

  script:
  n = seq.name
  result_file = "${n}.tab"
  """
  rpsblast -db cog/cog_db -query $seq -outfmt 6 -evalue 0.001 \
            -num_threads ${task.cpus} > ${result_file}
  """
}

process blast_swissprot {
  container "$params.blast_container"
  conda "$baseDir/conda/blast.yaml"

  input:
      tuple (file(swissprot_db), file(seq))

  output:
      path '*tab'

  script:

  n = seq.name
  """
  blastp -db $swissprot_db/swissprot.fasta -query ${n} \
  -outfmt 6 -evalue 0.001 -num_threads ${task.cpus} > ${n}.tab
  """
}

process diamond_refseq {
  container "$params.diamond_container"

  input:
  file(seq)

  output:
  path '*tab'

  script:

  n = seq.name
  """
  # new version of the database
  diamond blastp -p ${task.cpus} -d $params.refseq_db/merged_refseq.dmnd \
        -q ${n} -o ${n}.tab --max-target-seqs 200 -e 0.01 --max-hsps 1
  """
}

process execute_kofamscan {
  container "$params.kegg_container"
  conda "$baseDir/conda/kofamscan.yaml"

  input:
    tuple (file(ko_db), file(seq))

  output:
      path '*tab'

  script:
  n = seq.name

  """
  exec_annotation ${n} -p ${ko_db}/profiles/prokaryote.hal \
        -k ${ko_db}/ko_list --cpu ${task.cpus} -o ${n}.tab
  """
}

process prepare_amrscan {
  container "$params.ncbi_amr_container"
  conda "$baseDir/conda/amrfinderplus.yaml"

  output:
      path 'versions.txt'

  script:
  conda = params.conda
  """
  if $conda; then
      # amrfinder from conda comes without database so we need to download it
      # update (-u) fails on Linux because folder latest is missing so we need
      # to force update in that case.
      amrfinder -u || amrfinder -U
  fi
  amrfinder -V > versions.txt
  """
}

process execute_amrscan {
  container "$params.ncbi_amr_container"
  conda "$baseDir/conda/amrfinderplus.yaml"

  input:
    file(seq)
    file(version)

  output:
      file "amrfinder_results*.tab"

  script:
  n = seq.name
  """
  amrfinder --plus -p ${n} > amrfinder_results_${n}.tab
  """
}

process blast_vfdb {
  container "$params.blast_container"
  conda "$baseDir/conda/blast.yaml"

  input:
      tuple (path(vf_db), path(seq))

  output:
      path '*tab'

  script:

  n = seq.name
  """
  blastp -db $vf_db/vfdb.fasta -query ${n} \
  -outfmt "6 qaccver saccver pident length evalue bitscore qcovs" \
  -evalue ${params.vf_evalue} -qcov_hsp_perc ${params.vf_coverage} \
  -max_hsps 1 -num_threads ${task.cpus} > ${n}.tab
  """
}

process execute_islandpath {
    conda "$baseDir/conda/islandpath.yaml"
    container "$params.islandpath_container"

    input:
    path(genome)

    output:
    tuple (path(genome), path("*.gff"))

    script:
    
    n = genome.baseName
    """
    islandpath $genome ${n}.gff
    """
}

process gff_to_fasta {
    conda "$baseDir/conda/annotation.yaml"
    container "$params.annotation_container"

    input:
    tuple (path(genome), path(gff))

    output:
    path("*.fasta")

    script:
    """
    #!/usr/bin/env python
    import setup_chlamdb
    setup_chlamdb.gis_to_fasta("${genome}", "${gff}", "${genome.baseName}.fasta")
    """
}

process blast_gis {
  container "$params.blast_container"
  conda "$baseDir/conda/blast.yaml"

  input:
      tuple (path(gi_fasta), path(blast_db))

  output:
      path '*tab'

  script:

  n = gi_fasta.name
  """
  blastn -db $blast_db/merged -query $n \
  -outfmt "6 qaccver saccver pident length evalue bitscore qcovs sstart send" \
  -evalue 1.6e-7 -perc_identity 90 -qcov_hsp_perc 95 \
  -num_threads ${task.cpus} > ${n}.tab
  """
}

process extract_gis_hits {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path gis_predictions
        path gis_hits

    output:
        path result_file

    script:
    result_file = "all_gis.csv"
    """
        #!/usr/bin/env python
        import setup_chlamdb

        gis_predictions = "$gis_predictions".split()
        gis_hits = "$gis_hits".split()
        setup_chlamdb.extract_gis_hits(gis_predictions, gis_hits, "${result_file}")
    """
}

process gi_hits_to_fasta {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path gbk_files
        path gi_hits

    output:
        path result_file

    script:
    result_file = "all_gis.fasta"
    """
        #!/usr/bin/env python
        import setup_chlamdb
        gbk_files = "$gbk_files".split()
        setup_chlamdb.gi_hits_to_fasta(gbk_files, "${gi_hits}", "${result_file}")
    """
}

process compare_gis {
    label 'mount_basedir'
    container "$params.sourmash_container"
    conda "$baseDir/conda/sourmash.yaml"

    input:
        path gi_fastas

    output:
        path result_file

    script:
    result_file = "gi_comp.csv"
    """
    mkdir sigs
    sourmash sketch dna $gi_fastas --singleton
    sourmash compare ${gi_fastas}.sig --csv $result_file -p ${task.cpus}
    """
}

process setup_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    publishDir "${params.results_dir}/db"

    input:
        path base_dir

    output:
        path output_file

    script:
    output_file = "$workflow.runName"
    """
    cp $base_dir $output_file
    """
}

process load_base_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    // Necessary to prevent segfaults due to the large size
    // of the stage script
    // TODO: remove this by putting the numerous into a single directory
    beforeScript "ulimit -Ss unlimited"

    input:
        path db_base
        path input_file
        path gbks
        path orthofinder
        path alignments
        path nr_mapping_file
        path checkm_results
        path core_phylogeny
        path og_phylogenies_list
        path og_summary

    output:
        path db_name

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

    print("Loading groups", flush=True)
    setup_chlamdb.load_groups("$input_file", kwargs, "$db_base")

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

process load_refseq_results {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"
    input:
        path diamond_tsv_list
        path curr_db

    output:
        path "*_nr_hits.faa"
        path curr_db

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

process align_refseq_BBH_with_mafft {
  container "$params.mafft_container"
  conda "$baseDir/conda/mafft.yaml"

  input:
    path og

  output:
  path "*_mafft.faa"

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
  conda "$baseDir/conda/fasttree.yaml"

  input:
    path og

  output:
    path "*.nwk"

  script:
  """
    for orthogroup in $og;
    do
        FastTree \$orthogroup > \${orthogroup/.faa/.nwk}
    done
  """
}

process load_BBH_phylogenies {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db
        path BBH_phylogeny_trees

    output:
        path db

    script:
    """
    #!/usr/bin/env python

    import setup_chlamdb

    kwargs = ${gen_python_args()}

    BBH_list = "$BBH_phylogeny_trees".split(" ")
    setup_chlamdb.load_BBH_phylogenies(kwargs, BBH_list, "$db")
    """
}

process load_COG_into_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db
        path cog_file
        path cdd_to_cog
        path cog_db_dir

    output:
        path db

    script:
        """
        #!/usr/bin/env python

        import setup_chlamdb

        kwargs = ${gen_python_args()}
        cog_files = "${cog_file}".split()
        setup_chlamdb.load_cog(kwargs, cog_files, "$db", "$cdd_to_cog", "$cog_db_dir")
        """
}

process load_KO_into_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path KO_results
        path db
        path ko_db_dir

    output:
        path db

    script:
        """
        #!/usr/bin/env python

        kwargs = ${gen_python_args()}
        ko_files = "${KO_results}".split()
        import setup_chlamdb

        # this last function should be exported in a separate script to generate
        # the scaffold of a database
        setup_chlamdb.load_KO(kwargs, ko_files, "$db", "$ko_db_dir")
        setup_chlamdb.load_module_completeness(kwargs, "$db")
        """
}

process load_PFAM_info_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db
        path pfam_annot
        path pfam_dat
        path pfam_db

    output:
        path db

    script:
        """
        #!/usr/bin/env python
        import setup_chlamdb

        kwargs = ${gen_python_args()}
        pfam_files = "$pfam_annot".split()

        setup_chlamdb.load_pfam(kwargs, pfam_files, "$db", "$pfam_dat", "$pfam_db")
        """
}

process load_swissprot_hits_into_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db
        path blast_results
        path swissprot_db
        path swissprot_db_dir

    output:
        path db

    script:
    """
        #!/usr/bin/env python
        import setup_chlamdb

        kwargs = ${gen_python_args()}
        blast_results = "$blast_results".split()

        setup_chlamdb.load_swissprot(kwargs, blast_results, "$db", "swissprot.fasta", "$swissprot_db_dir")
    """
}

process load_amr_into_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        file collected_amr_files
        file db
        file version

    output:
        file db

    script:
        """
        #!/usr/bin/env python

        import setup_chlamdb

        kwargs = ${gen_python_args()}
        amr_files = "${collected_amr_files}".split()
        setup_chlamdb.load_amr(kwargs, amr_files, "$db", "$version")
        """
}

process load_vfdb_hits_into_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db
        path blast_results
        path vf_db_fasta
        path vf_db_defs

    output:
        path db

    script:
    min_seqid = params.vf_seqid
    """
        #!/usr/bin/env python
        import setup_chlamdb

        kwargs = ${gen_python_args()}
        blast_results = "$blast_results".split()

        setup_chlamdb.load_vfdb_hits(kwargs, blast_results, "$db", "$vf_db_fasta", "$vf_db_defs", float("$min_seqid"))
    """
}

process load_gis_into_db {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db
        path gis_predictions
        path gis_hits

    output:
        path db

    script:
    """
        #!/usr/bin/env python
        import setup_chlamdb

        kwargs = ${gen_python_args()}
        gis_predictions = "$gis_predictions".split()
        gis_hits = "$gis_hits".split()
        setup_chlamdb.load_gis(kwargs, gis_predictions, gis_hits, "$db")
    """
}

process create_chlamdb_search_index {
    label 'mount_basedir'
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    input:
        path db

    output:
        path index_name
        path db

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
        path index
        path db
        path alignments
        path results_dir
        path gbks

    script:
    custom_run_name=(params.name)?params.name:""
    db_path="${results_dir}/db/$workflow.runName"
    gbk_dir="${results_dir}/blast_DB/$workflow.runName/gbk"
    """
    mv $db_path ${db_path}_backup
    mv \$(readlink ${db_path}_backup) $db_path
    rm ${db_path}_backup

    if [ ! -d "${results_dir}/search_index" ]; then
        mkdir ${results_dir}/search_index/
    fi

    mkdir -p "${gbk_dir}"
    for gbk in filtered/*.gb*; do
    cp \$gbk ${gbk_dir}
    done

    mv \$(readlink $index) ${results_dir}/search_index/index_${workflow.runName}

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

workflow {

    dummy = check_reference_databases(Channel.fromPath(params.base_db))

    Channel.fromPath(params.input)
        .set { input_file }

    input_file
        .splitCsv(header: true, strip: true)
        .map { row -> file(row.file) }
        .set { gbk_files }

    (checked_gbks, to_cleanup_gbks) = check_gbk(dummy, gbk_files.collect(), input_file)

    (faa_files, ffn_files_seq, fna_files_SEQ) = convert_gbk(checked_gbks.flatten())

    merged_faa_files = faa_files.collectFile(name: 'merged.faa', newLine: true)

    (nr_seqs, nr_mapping_to_db_setup) = get_nr_sequences(merged_faa_files, faa_files.collect())

    split_nr_seqs = nr_seqs.splitFasta( by: 300, file: "chunk_" )

    fna_files_SEQ.collectFile(name: "merged.fna", newLine: true).set { merged_fna_makeblastdb }
    ffn_files_seq.collectFile(name: 'merged.ffn', newLine: true).set { merged_ffn_makeblastdb }

    fna_files_SEQ.mix(faa_files, ffn_files_seq, merged_ffn_makeblastdb,
                      merged_fna_makeblastdb, merged_faa_files).set { to_makeblastdb }

    makeblastdb(to_makeblastdb)

    (species_fasta, species_blastdb, result_dir) = prepare_orthofinder(faa_files.collect())
    blast_results = blast_orthofinder(result_dir, species_fasta, species_blastdb)
    (orthogroups, singletons) = orthofinder_main(result_dir, blast_results.collect())

    orthogroups_fasta = orthogroups2fasta(orthogroups, faa_files.collect())

    mafft_alignments = align_with_mafft(orthogroups_fasta.toSortedList().flatten().collate(20))

    mafft_alignments.toSortedList().flatten().collate(50).set { to_identity_calculation_split }

    to_load_alignment = identity_calculation(to_identity_calculation_split)


    mafft_alignments.collect().flatten().map { it }.filter { (it.text =~ /(>)/).size() > 2 }.set { alignement_larger_than_2_seqs }
    alignement_larger_than_2_seqs.toSortedList().flatten().collate(50).set { to_fasttree_orthogroups }
    gene_phylogeny = orthogroups_phylogeny_with_fasttree3(to_fasttree_orthogroups)


    (core_orthogroups, og_summary) = get_core_orthogroups(orthogroups, faa_files.collect(), mafft_alignments.collect())

    core_msa = concatenate_core_orthogroups(core_orthogroups.collect())

    core_genome_phylogeny = build_core_phylogeny_with_fasttree(core_msa)

    checkm_table = checkm_analyse(faa_files.collect())

    db = setup_db(Channel.fromPath("${params.zdb.file}"))
    db = load_base_db(db, input_file, checked_gbks, orthogroups, to_load_alignment.collect(), nr_mapping_to_db_setup, checkm_table, core_genome_phylogeny, gene_phylogeny.collect(), og_summary)

    if(params.pfam) {
        Channel.fromPath("${params.pfam_db}", type: "dir").set { pfam_db }
        pfam_db.combine(split_nr_seqs).set { to_pfam_scan_combined }
        pfam_results = pfam_scan(to_pfam_scan_combined)
        db = load_PFAM_info_db(db, pfam_results.collect(), Channel.fromPath("$params.pfam_db/Pfam-A.hmm.dat"), pfam_db)
    }

    if(params.cog) {
        Channel.fromPath("$params.cog_db", type: "dir").set { to_cog_multi }
        to_cog_multi.combine(split_nr_seqs).set { to_rpsblast_COG_multi }
        COG_to_load_db = rpsblast_COG(to_rpsblast_COG_multi)
        db = load_COG_into_db(db, COG_to_load_db.collect(), Channel.fromPath("$params.cog_db/cdd_to_cog"),  Channel.fromPath("$params.cog_db"))
    }

    if (params.blast_swissprot) {
        Channel.fromPath("$params.swissprot_db", type: "dir").set { to_swissprot_multi }
        to_swissprot_multi.combine(split_nr_seqs).set { to_blast_swissprot_multi }
        swissprot_blast = blast_swissprot(to_blast_swissprot_multi)
        db = load_swissprot_hits_into_db(db, swissprot_blast.collect(), Channel.fromPath("$params.swissprot_db/swissprot.fasta"), Channel.fromPath("$params.swissprot_db"))
    }

    if(params.diamond_refseq) {
        refseq_diamond = diamond_refseq(split_nr_seqs)
        refseq_diamond.collectFile().set { refseq_diamond_results_sqlitedb }
        (diamond_best_hits, db) = load_refseq_results(refseq_diamond_results_sqlitedb.collect(), db_gen)
        mafft_alignments_refseq_BBH = align_refseq_BBH_with_mafft(diamond_best_hits.flatten().collate( 20 ))
        BBH_phylogenies = orthogroup_refseq_BBH_phylogeny_with_fasttree(mafft_alignments_refseq_BBH)
        db = load_BBH_phylogenies(db, BBH_phylogenies.collect())
    }

    if(params.ko) {
        Channel.fromPath("$params.ko_db", type: "dir").set { to_ko_multi }
        to_ko_multi.combine(split_nr_seqs).set { to_kofamscan_multi }
        to_load_KO = execute_kofamscan(to_kofamscan_multi)
        db = load_KO_into_db(to_load_KO.collect(), db, Channel.fromPath("$params.ko_db"))
    }

    if(params.amr) {
        amr_version = prepare_amrscan()
        amr_table = execute_amrscan(split_nr_seqs, amr_version)
        db = load_amr_into_db(amr_table.collect(), db, amr_version)
    }

    if(params.vfdb) {
        vf_ref_db = Channel.fromPath("$params.vf_db", type: "dir")
        db_seq_combined = vf_ref_db.combine(split_nr_seqs)
        vfdb_blast = blast_vfdb(db_seq_combined)
        db = load_vfdb_hits_into_db(db, vfdb_blast.collect(), Channel.fromPath("$params.vf_db/vfdb.fasta"), Channel.fromPath("$params.vf_db/VFs.xls"))
    }

    if(params.gi) {
        genome_and_gis = execute_islandpath(checked_gbks.flatten())
        gi_fastas = gff_to_fasta(genome_and_gis)
        fna_db = channel.fromPath("${params.results_dir}/blast_DB/$workflow.runName/fna/")
        gi_fasta_and_fna_db = gi_fastas.combine(fna_db).combine(makeblastdb.out[1].collect()).map { tuple(it[0], it[1]) }
        gis_hits = blast_gis(gi_fasta_and_fna_db)
        extract_gis_hits(genome_and_gis.map { it[1] }.collect(), gis_hits.collect())
        gi_hits_to_fasta(checked_gbks.flatten().collect(), extract_gis_hits.out)
        compare_gis(gi_hits_to_fasta.out)
        // db = load_gis_into_db(db, genome_and_gis.collect().map { it[1] }, gis_hits.collect())
    }

    (to_index_cleanup, to_db_cleanup) = create_chlamdb_search_index(db)
    cleanup(to_index_cleanup, to_db_cleanup, mafft_alignments.collect(), Channel.fromPath("${params.results_dir}", type: "dir"), to_cleanup_gbks)
}

workflow.onComplete {
    println "Annotation pipeline completed"
    println "You may now start the web server with the zdb webapp command"
}
