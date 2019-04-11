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
params.blast_cog = false
params.orthofinder = true
params.interproscan = false
params.uniparc = true
params.blast_swissprot = false
params.string = true
params.core_missing = 0
params.genome_faa_folder = "$PWD/faa"
params.executor = 'local'

log.info "====================================="
log.info "input                  : ${params.input}"
log.info "Blast COG              : ${params.blast_cog}"
log.info "Orthofinder            : ${params.orthofinder}"
log.info "Orthofinder path       : ${params.genome_faa_folder}"
log.info "Core missing           : ${params.core_missing}"
log.info "Executor               : ${params.executor}"

Channel
  .fromPath(params.input)
  .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
  .into { faa_genomes1
          faa_genomes2
          faa_genomes3 }

Channel
    .fromPath(params.input)
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .collectFile(name: 'merged.faa', newLine: true)
    .splitFasta( by: 1000, file: "chunk_" )
    .into { faa_chunks1
            faa_chunks2
            faa_chunks3
            faa_chunks4}

Channel
    .fromPath(params.input)
    .ifEmpty { error "Cannot find any input sequence files matching: ${params.input}" }
    .collectFile(name: 'merged.faa', newLine: true)
    .into { merged_faa1
            merged_faa2 }



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

  cpus 2

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
  mafft \$faa > \${faa/.faa/_mafft.faa}
  done
  """
}

/* Get only alignments with more than than two sequences */
mafft_alignments.collect().into {all_alignments_1
                                 all_alignments_2
                                 all_alignments_3}

all_alignments_1.flatten().map { it }.filter { (it.text =~ /(>)/).size() > 3 }.set { large_alignments }
all_alignments_1.flatten().map { it }.filter { (it.text =~ /(>)/).size() == 3 }.set { small_alignments }
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

process orthogroups_phylogeny_with_iqtree_no_boostrap {

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

  script:
  """
  iqtree -nt 2 -s ${og} -pre ${og.getBaseName()}
  """
}

process get_core_orthogroups {

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'

  publishDir 'orthology/core_groups', mode: 'copy', overwrite: false

  input:
  file 'Orthogroups.txt' from orthogroups
  file genome_list from faa_genomes3.collect()
  file fasta_files from all_alignments_2.collect()

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
                                                                                              int(${params.core_missing}),
                                                                                              False)
  sequence_data = {}
  for fasta_file in "${fasta_files}".split(" "):
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

process concatenate_core_orthogroups {

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: false

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

  publishDir 'orthology/core_alignment_and_phylogeny', mode: 'copy', overwrite: false

  input:
  file 'msa.faa' from core_msa

  output:
  file 'core_genome_phylogeny.nwk'

  script:
  '''
  FastTree -gamma -spr 4 -mlacc 2 -slownni msa.faa > core_genome_phylogeny.nwk
  '''
}


process blast_COG {

  conda 'bioconda::blast=2.7.1'

  cpus 4

  when:
  params.blast_cog == true

  input:
  file 'seq' from faa_chunks1

  output:
  file 'blast_result' into blast_result

  script:
  n = seq.name
  """
  blastp -db $params.databases_dir/COG/prot2003-2014_test.fa -query seq -outfmt 6 -num_threads ${task.cpus} > blast_result
  """
}

blast_result.collectFile(name: 'annotation/COG/blast_COG.tab')

process blast_swissprot {

  conda 'bioconda::blast=2.7.1'

  publishDir 'annotation/blast_swissprot', mode: 'copy', overwrite: false

  cpus 4

  when:
  params.blast_swissprot == true

  input:
  file(seq) from faa_chunks3

  output:
  file '*tab' into swissprot_blast

  script:

  n = seq.name
  """
  # chunk_.6
  blastp -db $params.databases_dir/uniprot/swissprot/uniprot_sprot.fasta -query ${n} -outfmt 6 -num_threads ${task.cpus} > ${n}.tab
  """
}

process get_uniparc_mapping {

  conda 'bioconda::biopython=1.68'

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: false

  when:
  params.uniparc == true

  input:
  file(seq) from merged_faa1


  output:
  file 'uniparc_mapping.tab' into uniparc_mapping
  file 'uniprot_mapping.tab' into uniprot_mapping
  file 'no_uniprot_mapping.faa' into no_uniprot_mapping
  file 'no_uniparc_mapping.faa' into no_uniparc_mapping

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

uniparc_map.write("locus_tag\\tuniparc_id\\tuniparc_accession\\n")
uniprot_map.write("locus_tag\\tuniprot_accession\\ttaxon_id\\tdescription\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_mapping_uniprot_records = []
no_mapping_uniparc_records = []
for record in records:
    match = False
    sql = 'select t1.uniparc_id,uniparc_accession,accession,taxon_id,description, db_name from uniparc_accession t1 inner join uniparc_cross_references t2 on t1.uniparc_id=t2.uniparc_id inner join crossref_databases t3 on t2.db_id=t3.db_id where sequence_hash=? and status=1'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_mapping_uniparc_records.append(record)
        no_mapping_uniprot_records.append(record)
    else:
        uniparc_map.write("%s\\t%s\\t%s\\n" % (record.id,
                                            hits[0][0],
                                            hits[0][1]))
        for uniprot_hit in hits:
            if uniprot_hit[5] in ["UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL"]:
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
  """
}

process get_string_mapping {

  conda 'bioconda::biopython=1.68'

  publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: false

  when:
  params.string == true

  input:
  file(seq) from merged_faa2


  output:
  file 'string_mapping.tab' into string_mapping
  file 'no_string_mapping.faa' into no_string_mapping

  script:
  fasta_file = seq.name
  """
#!/usr/bin/env python

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("${params.databases_dir}/string/string_proteins.db")
cursor = conn.cursor()

fasta_file = "${fasta_file}"

string_map = open('string_mapping.tab', 'w')
no_string_mapping = open('no_string_mapping.faa', 'w')

string_map.write("locus_tag\\tstring_id\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_mapping_string_records = []
for record in records:
    sql = 'select accession from hash_table where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_mapping_string_records.append(record)
    else:
        for hit in hits:
          string_map.write("%s\\t%s\\n" % (record.id,
                                              hit[0]))


SeqIO.write(no_mapping_string_records, no_string_mapping, "fasta")

  """
}


process execute_interproscan {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: false

  when:
  params.interproscan == true

  input:
  file(seq) from faa_chunks2

  output:
  file '*gff3' into interpro_gff3
  file '*html.tar.gz' into interpro_html
  file '*svg.tar.gz' into interpro_svg
  file '*tsv' into interpro_tsv
  file '*xml' into interpro_xml

  script:
  n = seq.name
  """
  interproscan.sh -appl ProDom,HAMAP,TIGRFAM,SUPERFAMILY,PRINTS,PIRSF,COILS,ProSiteProfiles,PfamA,SMART,Phobius,SMART,SignalP_GRAM_NEGATIVE -goterms -pa -f TSV,XML,GFF3,HTML,SVG -i ${n} -d . -T . -iprlookup
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
