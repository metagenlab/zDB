#!/usr/bin/env nextflow
/*
 * Author:
 * - Trestan Pillonel <trestan.pillonel@gmail.com>
 *
 */


/*
 * Pipeline input params
 */

params.input = "faa/*.faa" 	// input sequences
log.info params.input
params.databases_dir = "$PWD/databases"
params.blast_cog = true
params.orthofinder = true
params.genome_faa_folder = "$PWD/faa"
params.executor = 'local'

log.info "====================================="
log.info "input                  : ${params.input}"
log.info "Blast COG              : ${params.blast_cog}"
log.info "Orthofinder            : ${params.orthofinder}"
log.info "Orthofinder path       : ${params.genome_faa_folder}"
log.info "Executor               : ${params.executor}"

  Channel
    .fromPath(params.input)
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .into { faa_genomes1
            faa_genomes2
            faa_genomes3 }

process prepare_orthofinder {

  echo true
  conda 'bioconda::orthofinder=2.2.7'

  input:
    file genome_list from faa_genomes1.collect()

  output:
    file 'Results_*/WorkingDirectory/Species*.fa' into species_fasta
    file 'Results_*/WorkingDirectory/BlastDBSpecies*.phr' into species_blastdb
    file 'Result*' into result_dir

  script:
  """
  orthofinder -op -a 8 -f . > of_prep.tab
  """
}

process blast_orthofinder {

  echo true

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
  blastp -outfmt 6 -evalue 0.001 -query $seq -db $blastdb_path/$blastdb_name > ${complete_dir}/WorkingDirectory/Blast${species_1}_${species_2}.txt
  """

}

process orthofinder_main {
  echo true
  conda 'bioconda::orthofinder=2.2.7'

  publishDir 'orthology', mode: 'copy', overwrite: false

  input:
  file complete_dir from result_dir
  file blast_results from blast_results.collect()

  output:
  file 'Results_*/WorkingDirectory/Orthogroups.txt' into orthogroups
  file 'Results_*/WorkingDirectory/SingleCopyOrthogroups.txt' into singletons

  script:
  """
  echo orthofinder -og -a 8 -b ./Results*/WorkingDirectory
  orthofinder -og -a 8 -b ./Results*/WorkingDirectory/ > of_grouping.txt
  """
}

process orthogroups2fasta {
  '''
  Get fasta file of each orthogroup
  '''

  echo true
  conda 'bioconda::biopython=1.70'

  publishDir 'orthology/orthogroups_fasta', mode: 'copy', overwrite: false

  input:
  file 'Orthogroups.txt' from orthogroups
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

/*rthogroups_fasta.flatten().collate( 2 ).set {chunks}
*/

process align_with_mafft {

  echo true
  conda 'bioconda::mafft=7.407'

  publishDir 'orthology/orthogroups_alignments', mode: 'copy', overwrite: false

  input:
  file og from orthogroups_fasta.flatten().collate( 20 )

  output:
  file "*_mafft.faa" into mafft_alignments

  script:
  """
  unset MAFFT_BINARIES
  for faa in ${og}; do
  echo mafft \$faa  \${faa/.faa/_mafft.faa};
  mafft \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}

/* Get only alignments with more than than two sequences */
mafft_alignments.collect().flatten().map { it }.filter { (it.text =~ /(>)/).size() > 3 }.set { large_alignments }

/*
process orthogroups_phylogeny_with_raxml {

  echo true
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

process orthogroups_phylogeny_with_iqtree {

  conda 'bioconda::iqtree=1.6.8'
  cpus 2
  publishDir 'orthology/orthogroups_phylogenies_iqtree', mode: 'copy', overwrite: false

  input:
  each file(og) from large_alignments

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

process get_core_orthogroups {

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'

  publishDir 'orthology/core_groups', mode: 'copy', overwrite: false

  input:
  file 'Orthogroups.txt' from orthogroups
  file genome_list from faa_genomes3.collect()

  output:
  file '*_taxon_ids.faa' into core_orthogroups

  script:

  """
  #!/usr/bin/env python

  from Bio import SeqIO
  import os
  import pandas as pd

  def orthofinder2core_groups(fasta_list, mcl_file, n_missing=0,orthomcl=False):
    n_missing = 0

    orthogroup2locus_list = {}

    with open(mcl_file, 'r') as f:
        all_grp = [i for i in f]
        for n, line in enumerate(all_grp):
            if orthomcl:
                groups = line.rstrip().split('\t')
                groups = [i.split('|')[1] for i in groups]
            else:
                groups = line.rstrip().split(' ')
                groups = groups[1:len(groups)]
            orthogroup2locus_list["group_%s" % n] = groups

    locus2genome = {}

    for fasta in fasta_list:
        genome = os.path.basename(fasta).split('.')[0]
        for seq in SeqIO.parse(fasta, "fasta"):
            locus2genome[seq.name] = genome

    df = pd.DataFrame(index=orthogroup2locus_list.keys(), columns=set(locus2genome.values()))
    df = df.fillna(0)

    for group in orthogroup2locus_list:
        genome2count = {}
        for locus in orthogroup2locus_list[group]:
            if locus2genome[locus] not in genome2count:
                genome2count[locus2genome[locus]] = 1
            else:
                genome2count[locus2genome[locus]] += 1

        for genome in genome2count:
            df.loc[group, genome] = genome2count[genome]
    df =df.apply(pd.to_numeric, args=('coerce',))

    n_genomes = len(set(locus2genome.values()))
    n_minimum_genomes = n_genomes-n_missing
    freq_missing = (n_genomes-float(n_missing))/n_genomes
    limit = freq_missing*n_genomes
    print ('limit', limit)

    groups_with_paralogs = df[(df > 1).sum(axis=1) > 0].index
    df = df.drop(groups_with_paralogs)

    core_groups = df[(df == 1).sum(axis=1) >= limit].index.tolist()

    return df, core_groups, orthogroup2locus_list, locus2genome


  orthology_table, core_groups, orthogroup2locus_list, locus2genome = orthofinder2core_groups("${genome_list}".split(" "),
                                                                                              'Orthogroups.txt',
                                                                                              0,
                                                                                              False)
  sequence_data = {}
  for fasta_file in "${genome_list}".split(" "):
      sequence_data.update(SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta")))

  for one_group in core_groups:
    dest = '%s_taxon_ids.faa' % one_group
    new_fasta = []
    for locus in orthogroup2locus_list[one_group]:
        tmp_seq = sequence_data[locus]
        tmp_seq.name = locus2genome[locus]
        tmp_seq.id = locus2genome[locus]
        tmp_seq.description = locus2genome[locus]
        new_fasta.append(tmp_seq)

    out_handle = open(dest, 'w')
    SeqIO.write(new_fasta, out_handle, "fasta")
    out_handle.close()

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
