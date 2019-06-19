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
params.cog = true
params.orthofinder = true
params.interproscan = true
params.uniparc = true
params.uniprot_data = false
params.tcdb = true
params.blast_swissprot = true
params.plast_refseq = false
params.diamond_refseq = true
params.diamond_refseq_taxonomy = true
params.refseq_diamond_BBH_phylogeny = true
params.refseq_diamond_BBH_phylogeny_top_n_hits = 4
params.refseq_diamond_BBH_phylogeny_phylum_filter = '["Chlamydiae", "Verrucomicrobia", "Planctomycetes", "Kiritimatiellaeota", "Lentisphaerae"]'
params.string = true
params.pdb = true
params.oma = true
params.ko = true
params.tcdb_gblast = false
params.PRIAM = false
params.orthogroups_phylogeny_with_iqtree = false
params.orthogroups_phylogeny_with_fasttree = true
params.core_missing = 6
params.genome_faa_folder = "$PWD/faa"
params.executor = 'local'

params.local_sample_sheet = false //"local_assemblies.tab"
params.ncbi_sample_sheet = "ncbi_assemblies.tab"

log.info "====================================="
log.info "input                  : ${params.input}"
log.info "COG                    : ${params.cog}"
log.info "Orthofinder            : ${params.orthofinder}"
log.info "Orthofinder path       : ${params.genome_faa_folder}"
log.info "Core missing           : ${params.core_missing}"
log.info "Executor               : ${params.executor}"



// Each Sample
if (params.ncbi_sample_sheet != false){
  Channel.fromPath( file(params.ncbi_sample_sheet) )
                      .splitCsv(header: true, sep: '\t')
                      .map{row ->
                          // get the list of accessions
                          def assembly_accession = row."Genbank"
                          return "${assembly_accession}"
                      }
                      .into{
                          assembly_accession_list
                      }
}
if (params.local_sample_sheet != false){
  Channel.fromPath( file(params.local_sample_sheet) )
                      .splitCsv(header: true, sep: '\t')
                      .map{row ->
                          // get the list of accessions
                          def gbk_path = row."gbk_path"
                          return gbk_path
                      }
                      .map { file(it) }
                      .set { local_gbk_list }
}

// only define process if nedded
if (params.local_sample_sheet != false){
  process copy_local_assemblies {

    publishDir 'data/gbk_local', mode: 'copy', overwrite: true

    cpus 1

    when:
    params.local_sample_sheet != false

    input:
    file(local_gbk) from local_gbk_list

    output:
    file "${local_gbk.name}.gz" into raw_local_gbffs

    script:

    """
    gzip -f ${local_gbk.name}
    """
  }
}

// only define process if nedded
if (params.ncbi_sample_sheet != false){
  process download_assembly {

    conda 'bioconda::biopython=1.68'

    publishDir 'data/gbk_ncbi', mode: 'copy', overwrite: true

    when:
    params.ncbi_sample_sheet != false

    input:
    val assembly_accession_list from assembly_accession_list.collect()

    cpus 1

    output:
    file '*.gbff.gz' into raw_ncbi_gbffs

    script:
    //accession = assembly_accession[0]
    """
  #!/usr/bin/env python

  import re
  from ftplib import FTP
  from Bio import Entrez, SeqIO
  Entrez.email = "trestan.pillonel@chuv.ch"
  Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"

  accession_list = "${assembly_accession_list}".split(' ')
  for accession in accession_list:
    handle1 = Entrez.esearch(db="assembly", term="%s" % accession)
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
      raise("%s assembly not annotated! --- exit ---" % accession)
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
}

// merge local and ncbi gbk into a single channel
if (params.ncbi_sample_sheet != false && params.local_sample_sheet == false) {
println "ncbi"
raw_ncbi_gbffs.collect().into{all_raw_gbff}
}
else if(params.ncbi_sample_sheet == false && params.local_sample_sheet != false) {
println "local"
raw_local_gbffs.collect().into{all_raw_gbff}
}
else {
println "both"
raw_ncbi_gbffs.mix(raw_local_gbffs).into{all_raw_gbff}
}

process gbk_check {

  publishDir 'data/gbk_edited', mode: 'copy', overwrite: true

  conda 'bioconda::biopython=1.68'

  cpus 2

  input:
  file(all_gbff) from all_raw_gbff.collect()

  output:
  file "*merged.gbk" into edited_gbks

  script:
  println all_gbff
  """
  gbff_check.py -i ${all_gbff} -l 1000
  """
}

process convert_gbk_to_faa {

  publishDir 'data/faa_locus', mode: 'copy', overwrite: true

  conda 'bioconda::biopython=1.68'

  echo false

  cpus 1

  input:
  each file(edited_gbk) from edited_gbks

  output:
  file "*.faa" into faa_files

  script:
  """
#!/usr/bin/env python
print("${edited_gbk}")
from Bio import Entrez, SeqIO

records = SeqIO.parse("${edited_gbk}", 'genbank')
edited_records = open("${edited_gbk.baseName}.faa", 'w')
for record in records:
  protein2count = {}
  for feature in record.features:
      if feature.type == 'CDS' and 'pseudo' not in feature.qualifiers and 'pseudogene' not in feature.qualifiers:
          try:
            locus_tag = feature.qualifiers["locus_tag"][0]
          except KeyError:
            protein_id = feature.qualifiers["protein_id"][0].split(".")[0]
            if protein_id not in protein2count:
                protein2count[protein_id] = 1
                locus_tag = protein_id
            else:
                protein2count[protein_id] += 1
                locus_tag = "%s_%s" % (protein_id, protein2count[protein_id])
          try:
            edited_records.write(">%s %s\\n%s\\n" % (locus_tag,
                                                     record.description,
                                                     feature.qualifiers['translation'][0]))
          except KeyError:
              print("problem with feature:", feature)
  """
}

faa_files.into{ faa_locus1
                faa_locus2
              }


faa_locus1.into { faa_genomes1
                 faa_genomes2
                 faa_genomes3 }

faa_locus2.collectFile(name: 'merged.faa', newLine: true)
    .into { merged_faa0 }


process get_nr_sequences {

  conda 'bioconda::biopython=1.68'

  publishDir 'data/', mode: 'copy', overwrite: true

  input:
  file(seq) from merged_faa0

  output:

  file 'nr.faa' into nr_seqs
  file 'nr_mapping.tab' into nr_mapping

  script:
  fasta_file = seq.name
  """
#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqUtils import CheckSum

fasta_file = "${fasta_file}"

nr_fasta = open('nr.faa', 'w')
nr_mapping = open('nr_mapping.tab', 'w')

checksum_nr_list = []

records = SeqIO.parse(fasta_file, "fasta")
updated_records = []

for record in records:

    checksum = CheckSum.crc64(record.seq)
    nr_mapping.write("%s\\t%s\\n" % (record.id,
                                   checksum))
    if checksum not in checksum_nr_list:
      checksum_nr_list.append(checksum)
      record.id = checksum
      record.name = ""
      updated_records.append(record)

SeqIO.write(updated_records, nr_fasta, "fasta")

  """
}

nr_seqs.collectFile(name: 'merged_nr.faa', newLine: true)
.into { merged_faa_chunks
        merged_faa1
        merged_faa2
        merged_faa3
        merged_faa4
        merged_faa5
        merged_faa6 }

merged_faa_chunks.splitFasta( by: 1000, file: "chunk_" )
.into { faa_chunks1
        faa_chunks2
        faa_chunks3
        faa_chunks4
        faa_chunks5
        faa_chunks6
        faa_chunks7
        faa_chunks8 }

process prepare_orthofinder {

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

  conda 'bioconda::blast=2.7.1'

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

  conda 'bioconda::orthofinder=2.2.7'

  publishDir 'orthology', mode: 'copy', overwrite: true
  echo true

  input:
  file complete_dir from result_dir
  file blast_results from blast_results.collect()

  output:
  file 'Results_*/WorkingDirectory/Orthogroups.txt' into orthogroups
  file 'Results_*/WorkingDirectory/SingleCopyOrthogroups.txt' into singletons

  script:
  """
  echo "${complete_dir.baseName}"
  orthofinder -og -a 8 -b ./Results*/WorkingDirectory/ > of_grouping.txt
  """
}

orthogroups
.into { orthogroups_1
        orthogroups_2}

process orthogroups2fasta {
  '''
  Get fasta file of each orthogroup
  '''

  conda 'bioconda::biopython=1.70'

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

  conda 'bioconda::mafft=7.407'

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

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'

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

  input:
  file 'msa.faa' from core_msa

  output:
  file 'core_genome_phylogeny.nwk'

  script:
  '''
  FastTree -gamma -spr 4 -mlacc 2 -slownni msa.faa > core_genome_phylogeny.nwk
  '''
}


process rpsblast_COG {

  conda 'bioconda::blast=2.7.1'

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
  rpsblast -db $params.databases_dir/cdd/Cog -query seq -outfmt 6 -evalue 0.001 -num_threads ${task.cpus} > blast_result
  """
}

blast_result.collectFile(name: 'annotation/COG/blast_COG.tab')

process blast_swissprot {

  conda 'bioconda::blast=2.7.1'

  publishDir 'annotation/blast_swissprot', mode: 'copy', overwrite: true

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
  blastp -db $params.databases_dir/uniprot/swissprot/uniprot_sprot.fasta -query ${n} -outfmt 6 -evalue 0.001  -num_threads ${task.cpus} > ${n}.tab
  """
}


process plast_refseq {

  publishDir 'annotation/plast_refseq', mode: 'copy', overwrite: true

  cpus 4

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
  /home/tpillone/work/dev/annotation_pipeline_nextflow/bin/plast -p plastp -a ${task.cpus} -d $params.databases_dir/refseq/merged.faa.pal -i ${n} -M BLOSUM62 -s 75 -seeds-use-ratio 20 -max-database-size 50000000 -e 1e-5 -G 11 -E 1 -o ${n}.tab -F F -bargraph -verbose -force-query-order 1000 -max-hit-per-query 100 -max-hsp-per-hit 1 > ${n}.log;
  """
}

process diamond_refseq {

  publishDir 'annotation/diamond_refseq', mode: 'copy', overwrite: true

  cpus 4
  conda 'bioconda::diamond=0.9.24'

  when:
  params.diamond_refseq == true

  input:
  file(seq) from faa_chunks6

  output:
  file '*tab' into refseq_diamond

  script:

  n = seq.name
  """
  diamond blastp -p ${task.cpus} -d $params.databases_dir/refseq/merged_refseq.dmnd -q ${n} -o ${n}.tab --max-target-seqs 100 -e 0.01 --max-hsps 1
  """
}

refseq_diamond.collectFile()
.into { refseq_diamond_results_taxid_mapping
        refseq_diamond_results_sqlitedb }

refseq_diamond_results_taxid_mapping
.splitCsv(header: false, sep: '\t')
.map{row ->
    def protein_accession = row[1]
    return "${protein_accession}"
}
.unique()
.collectFile(name: 'nr_refseq_hits.tab', newLine: true)
.set {refseq_diamond_nr}

//.collate( 300 )
//.set {
//    nr_refseq_hits_chunks
//}

process get_refseq_hits_taxonomy {

  conda 'biopython=1.73=py36h7b6447c_0'

  publishDir 'annotation/diamond_refseq/', mode: 'copy', overwrite: true

  echo false

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
import sqlite3
conn = sqlite3.connect("refseq_taxonomy.db")
cursor = conn.cursor()
from Bio.SeqUtils import CheckSum
from Bio import Entrez
Entrez.email = "trestan.pillonel@chuv.ch"
Entrez.api_key = "719f6e482d4cdfa315f8d525843c02659408"

sql = 'create table refseq_hits (accession varchar(200), taxid INTEGER, description TEXT, length INTEGER)'
cursor.execute(sql,)
conn.commit()

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i+n]

f = open("${refseq_hit_table}", 'r')
hit_list = [i.rstrip() for i in f]

accession_lists = chunks(hit_list, 300)

def accession2taxon(gi_list, database):
    from socket import error as SocketError
    import errno

    hit_annotation = []

    try:
        handle = Entrez.esummary(db=database, id=','.join(gi_list), retmax=len(gi_list))
    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise('error connexion with %s' % ','.join(gi_list))
        else:
            import time
            print ('connexion error, trying again...')
            time.sleep(60)
            accession2taxon(gi_list, database)

    if isinstance(gi_list, list) and len(gi_list) == 1:
        record = Entrez.read(handle, validate=False)
        hit_annotation.append([
                               record['AccessionVersion'],
                               record['TaxId'],
                               record['Title'],
                               record['Length']
                               ]
                               )
        return hit_annotation
    else:
        record = Entrez.parse(handle, validate=False)
        try:
            for i in record:
                hit_annotation.append([
                                       i['AccessionVersion'],
                                       i['TaxId'],
                                       i['Title'],
                                       i['Length']
                                       ])
        except RuntimeError:
            print ('RuntimeError error, trying again... %s' % gi_list)
            accession2taxon(gi_list, database=database)
        except http.client.IncompleteRead:
            print ('IncompleteRead error, trying again...%s' % gi_list)
            accession2taxon(gi_list, database=database)
        return hit_annotation

for n, one_list in enumerate(accession_lists):
  data = accession2taxon(one_list, "protein")
  sql_template = 'insert into refseq_hits values (?,?,?,?)'
  cursor.executemany(sql_template, data)
  conn.commit()

# index accession and taxid columns
sql_index_1 = 'create index acc on refseq_hits (accession);'
sql_index_2 = 'create index taxid on refseq_hits (taxid);'

cursor.execute(sql_index_1)
cursor.execute(sql_index_2)
conn.commit()

  """
}


process get_uniparc_mapping {

  conda 'bioconda::biopython=1.68'

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


process get_uniprot_data {

  conda 'biopython=1.73=py36h7b6447c_0'

  publishDir 'annotation/uniparc_mapping/', mode: 'copy', overwrite: true
  echo false

  when:
  params.uniprot_data == true

  input:
  file(table) from uniprot_mapping_tab

  output:
  file 'uniprot_data.tab' into uniprot_data

  script:

  """
#!/usr/bin/env python3.6

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

def uniprot_record2annotations(record):

    data2info = {}
    if 'recommendedName_ecNumber' in record.annotations:
        data2info['ec_number'] = record.annotations['recommendedName_ecNumber'][0]
    else:
        data2info['ec_number'] = '-'
    if 'gene_name_primary' in record.annotations:
        data2info['gene'] = record.annotations['gene_name_primary']
    else:
        data2info['gene'] = '-'
    if 'comment_catalyticactivity' in record.annotations:
        data2info['comment_catalyticactivity'] =  record.annotations['comment_catalyticactivity'][0]
    else:
        data2info['comment_catalyticactivity'] = '-'
    if 'comment_subunit' in record.annotations:
        data2info['comment_subunit'] =  record.annotations['comment_subunit'][0]
    else:
        data2info['comment_subunit'] = '-'
    if 'comment_function' in record.annotations:
        data2info['comment_function'] =  record.annotations['comment_function'][0]
    else:
        data2info['comment_function'] = '-'
    if 'recommendedName_fullName' in record.annotations:
        data2info['recommendedName_fullName'] =  record.annotations['recommendedName_fullName'][0]
    else:
        data2info['recommendedName_fullName'] = '-'
    if 'comment_similarity' in record.annotations:
        data2info['comment_similarity'] =  record.annotations['comment_similarity'][0]
    else:
        data2info['comment_similarity'] = '-'
    if 'proteinExistence' in record.annotations:
        data2info['proteinExistence'] =  record.annotations['proteinExistence'][0]
    else:
        data2info['proteinExistence'] = '-'
    if 'keywords' in record.annotations:
        data2info['keywords'] =  ';'.join(record.annotations['keywords'])
    else:
        data2info['keywords'] = '-'
    if 'comment_pathway' in record.annotations:
        data2info['comment_pathway'] = record.annotations['comment_pathway'][0]
    else:
        data2info['comment_pathway'] = '-'
    if 'comment_developmentalstage' in record.annotations:
        data2info['developmentalstage'] = record.annotations['comment_developmentalstage'][0]
    else:
        data2info['developmentalstage'] = '-'
    data2info['proteome'] = '-'
    for ref in record.dbxrefs:
        if "Proteomes" in ref:
            data2info['proteome'] = ref.split(":")[1]

    return data2info

def uniprot_accession2go_and_status(uniprot_accession):

    import urllib.request
    from urllib.error import URLError

    uniprot_accession = uniprot_accession.split(".")[0]

    link = "http://www.uniprot.org/uniprot/?query=%s&columns=annotation%%20score,reviewed&format=tab" % (uniprot_accession)

    try:
        page = urllib.request.urlopen(link)
        data = page.read().decode('utf-8').split('\\n')
        rows = [i.rstrip().split('\\t') for i in data]
    except URLError:
        success = False
        while not success:
            import time
            print ('connection problem, trying again...')
            time.sleep(10)
            try:
                page = urllib.request.urlopen(req)
                data = page.read().decode('utf-8').split('\\n')
                rows = [i.rstrip().split('\\t') for i in data]
                success=True
            except:
                success = False
    print(uniprot_accession, rows)
    return (rows[1][0][0], rows[1][1])

def uniprot_id2record(uniprot_accession, n_trial=0):
    link = "http://www.uniprot.org/uniprot/%s.xml" % (uniprot_accession)
    from io import StringIO
    from Bio import SeqIO
    from urllib.error import URLError
    import urllib.request
    try:
        data = urllib.request.urlopen(link).read().decode('utf-8')
    except URLError:
        print ('echec', link)
        return False
    rec = StringIO(data)
    try:
        record = SeqIO.read(rec, 'uniprot-xml')
    except:
        import time
        print ('problem with %s, trying again, %s th trial' % (uniprot_accession, n_trial))
        time.sleep(50)
        n_trial+=1
        if n_trial < 5:
            record = uniprot_id2record(uniprot_accession, n_trial=n_trial)
        else:
            return False
    return record



conn = sqlite3.connect("${params.databases_dir}/uniprot/uniparc/uniparc.db")
cursor = conn.cursor()

uniprot_table = open("${table}", 'r')
uniprot_data = open('uniprot_data.tab', 'w')

uniprot_data.write("uniprot_accession\\tuniprot_score\\tuniprot_status\\tproteome\\tcomment_function\\tec_number\\tcomment_subunit\\tgene\\trecommendedName_fullName\\tproteinExistence\\tdevelopmentalstage\\tcomment_similarity\\tcomment_catalyticactivity\\tcomment_pathway\\tkeywords\\n")

for row in uniprot_table:
    data = row.rstrip().split("\\t")
    uniprot_accession = data[1]
    if uniprot_accession == 'uniprot_accession':
        continue
    uniprot_accession = uniprot_accession.split(".")[0]
    print("uniprot_accession",uniprot_accession)
    try:
        uniprot_score, uniprot_status = uniprot_accession2go_and_status(uniprot_accession)
        uniprot_record = uniprot_id2record(uniprot_accession)
        annotation = uniprot_record2annotations(uniprot_record)
    # deal with eventual removed entries
    except IndexError:
        uniprot_score = 0
        uniprot_status = 'Removed'
        annotation["proteome"] = '-'
        annotation["comment_function"] = '-'
        annotation["ec_number"] = '-'
        annotation["comment_subunit"] = '-'
        annotation["gene"] = '-'
        annotation["recommendedName_fullName"] = '-'
        annotation["proteinExistence"] = '-'
        annotation["developmentalstage"] = '-'
        annotation["comment_similarity"] = '-'
        annotation["comment_catalyticactivity"] = '-'
        annotation["comment_pathway"] = '-'
        annotation["keywords"] = '-'

    uniprot_data.write("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" % ( uniprot_accession,
                                                                                                         uniprot_score,
                                                                                                         uniprot_status,
                                                                                                         annotation["proteome"],
                                                                                                         annotation["comment_function"],
                                                                                                         annotation["ec_number"],
                                                                                                         annotation["comment_subunit"],
                                                                                                         annotation["gene"],
                                                                                                         annotation["recommendedName_fullName"],
                                                                                                         annotation["proteinExistence"],
                                                                                                         annotation["developmentalstage"],
                                                                                                         annotation["comment_similarity"],
                                                                                                         annotation["comment_catalyticactivity"],
                                                                                                         annotation["comment_pathway"],
                                                                                                         annotation["keywords"])
                                                                                                        )
  """
}



process get_string_mapping {

  conda 'bioconda::biopython=1.68'

  publishDir 'annotation/string_mapping/', mode: 'copy', overwrite: true

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

process get_string_PMID_mapping {

  conda 'bioconda::biopython=1.68'

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

import urllib2

def string_id2pubmed_id_list(accession):

    link = 'http://string-db.org/api/tsv/abstractsList?identifiers=%s' % accession
    print link
    try:
        data = urllib2.urlopen(link).read().rstrip().decode('utf-8').split('\\n')[1:]
    except urllib2.URLError:
        print ('echec', link)
        return False
    pid_list = [row.split(':')[1] for row in data]
    print ('list', pid_list)
    return pid_list

o = open("string_mapping_PMID.tab", "w")

string_mapping = "${string_map}"

with open(string_mapping, 'r') as f:
    for n, row in enumerate(f):
        if n == 0:
            continue
        else:
            data = row.rstrip().split("\t")
            pmid_list = string_id2pubmed_id_list(data[1])
            if pmid_list:
                for id in pmid_list:
                    o.write("%s\\t%s\\n" % (data[0], id))
            else:
                o.write("%s\\tNone\\n" % (data[0]))

  """
}


process get_tcdb_mapping {

  conda 'bioconda::biopython=1.68'

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

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("${params.databases_dir}/TCDB/tcdb.db")
cursor = conn.cursor()

fasta_file = "${fasta_file}"

tcdb_map = open('tcdb_mapping.tab', 'w')
no_tcdb_mapping = open('no_tcdb_mapping.faa', 'w')

tcdb_map.write("locus_tag\\ttcdb_id\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_tcdb_mapping_records = []
for record in records:
    sql = 'select accession from hash_table where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_tcdb_mapping_records.append(record)
    else:
        for hit in hits:
          tcdb_map.write("%s\\t%s\\n" % (record.id,
                                              hit[0]))


SeqIO.write(no_tcdb_mapping_records, no_tcdb_mapping, "fasta")

  """
}

no_tcdb_mapping.splitFasta( by: 1000, file: "chunk" )
.set { faa_tcdb_chunks }

process tcdb_gblast3 {

  publishDir 'annotation/tcdb_mapping', mode: 'copy', overwrite: true

  cpus 1
  maxForks 1
  conda 'anaconda::biopython=1.67=np111py27_0 conda-forge::matplotlib=2.2.3 biobuilds::fasta'
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
  export "PATH=\$HMMTOP_PATH:\$GBLAST3_PATH:\$PATH"
  gblast3.py -i ${seq} -o TCDB_RESULTS_${seq}
  """
}

process get_pdb_mapping {

  conda 'bioconda::biopython=1.68'

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

from Bio import SeqIO
import sqlite3
from Bio.SeqUtils import CheckSum

conn = sqlite3.connect("${params.databases_dir}/pdb/pdb.db")
cursor = conn.cursor()

fasta_file = "${fasta_file}"

pdb_map = open('pdb_mapping.tab', 'w')
no_pdb_mapping = open('no_pdb_mapping.faa', 'w')

pdb_map.write("locus_tag\\tpdb_id\\n")

records = SeqIO.parse(fasta_file, "fasta")
no_pdb_mapping_records = []
for record in records:
    sql = 'select accession from hash_table where sequence_hash=?'
    cursor.execute(sql, (CheckSum.seguid(record.seq),))
    hits = cursor.fetchall()
    if len(hits) == 0:
        no_pdb_mapping_records.append(record)
    else:
        for hit in hits:
          pdb_map.write("%s\\t%s\\n" % (record.id,
                                              hit[0]))


SeqIO.write(no_pdb_mapping_records, no_pdb_mapping, "fasta")

  """
}

process get_oma_mapping {

  conda 'bioconda::biopython=1.68'

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

no_uniparc_mapping_faa.splitFasta( by: 300, file: "no_uniparc_match_chunk_" )
.set { no_uniparc_match_chunks }
uniparc_mapping_faa.splitFasta( by: 1000, file: "uniparc_match_chunk_" )
.set { uniparc_match_chunks }

process execute_interproscan_no_uniparc_matches {

  publishDir 'annotation/interproscan', mode: 'copy', overwrite: true

  cpus 8
  memory '16 GB'
  conda 'anaconda::openjdk=8.0.152'

  when:
  params.interproscan == true

  input:
  file(seq) from no_uniparc_match_chunks

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
  file(seq) from uniparc_match_chunks

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

  conda 'hmmer=3.2.1 parallel ruby=2.4.5'

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
  export "PATH=\$KOFAMSCAN_HOME:\$PATH"
  exec_annotation ${n} -p ${params.databases_dir}/kegg/profiles/prokaryote.hal -k ${params.databases_dir}/kegg/ko_list --cpu ${task.cpus} -o ${n}.tab
  """
}


process execute_PRIAM {

  publishDir 'annotation/KO', mode: 'copy', overwrite: true

  cpus 4
  memory '8 GB'

  when:
  params.PRIAM == true

  input:
  file(seq) from faa_chunks8

  output:
  file '*tab'

  script:
  n = seq.name
  """
  export "PATH=\$KOFAMSCAN_HOME:\$PATH"
  java -jar  databases/PRIAM/PRIAM_search.jar -i group_333.faa -o test -p databases/PRIAM/PRIAM_JAN18
  """
}


process setup_orthology_db {
  /*
  - [ ] load locus_tag2hash table
  - [ ] load locus2orthogroup data into sqlite db (locus2orthogroup)
  */

  conda 'biopython=1.73=py36h7b6447c_0'
  publishDir 'orthology/', mode: 'link', overwrite: true

  cpus 4
  memory '8 GB'

  when:
  params.refseq_diamond_BBH_phylogeny == true

  input:
  file nr_mapping_file from nr_mapping
  file orthogroup from orthogroups_2
  file nr_fasta from merged_faa6

  output:
  file 'orthology.db' into orthology_db

  script:
  """
  #!/usr/bin/env python


  from Bio import SeqIO
  import sqlite3
  from Bio.SeqUtils import CheckSum

  fasta_dict = SeqIO.to_dict(SeqIO.parse("${nr_fasta}", "fasta"))

  conn = sqlite3.connect("orthology.db")
  cursor = conn.cursor()

  # sequence table
  sql0 = 'create table sequence_hash2aa_sequence (sequence_hash binary, sequence TEXT )'
  cursor.execute(sql0,)
  sql = 'insert into  sequence_hash2aa_sequence values (?, ?)'
  for hash in fasta_dict:
    cursor.execute(sql, (hash, str(fasta_dict[hash].seq)))

  # hash mapping table
  sql1 = 'create table locus_tag2sequence_hash (locus_tag varchar(200), sequence_hash binary)'
  cursor.execute(sql1,)

  sql = 'insert into locus_tag2sequence_hash values (?,?)'
  with open("${nr_mapping_file}", 'r') as f:
      for row in f:
          data = row.rstrip().split("\\t")
          cursor.execute(sql, data)
  conn.commit()

  # orthogroup table
  sql2 = 'create table locus_tag2orthogroup (locus_tag varchar(200), orthogroup varchar(200))'
  cursor.execute(sql2,)
  sql = 'insert into locus_tag2orthogroup values (?, ?)'
  with open("${orthogroup}", 'r') as f:
      for row in f:
          data = row.rstrip().split(" ")
          for locus in data[1:]:
            cursor.execute(sql,(locus, data[0][0:-1]))
  conn.commit()

  # index hash, locus and orthogroup columns
  sql_index_1 = 'create index hash1 on sequence_hash2aa_sequence (sequence_hash);'
  sql_index_2 = 'create index hash2 on locus_tag2sequence_hash (sequence_hash);'
  sql_index_3 = 'create index locus1 on locus_tag2sequence_hash (locus_tag);'
  sql_index_4 = 'create index locus2 on locus_tag2orthogroup (locus_tag);'
  sql_index_5 = 'create index og on locus_tag2orthogroup (orthogroup);'

  cursor.execute(sql_index_1)
  cursor.execute(sql_index_2)
  cursor.execute(sql_index_3)
  cursor.execute(sql_index_4)
  cursor.execute(sql_index_5)
  conn.commit()
  """
}


process setup_diamond_refseq_db {
  /*
  - [ ] load blast results into sqlite db
  */

  conda 'bioconda::biopython=1.68 anaconda::pandas=0.23.4'
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
sql_index_2 = 'create index acc on diamond_refseq (qseqid);'

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
  mafft \$faa > \${faa/.faa/_mafft.faa}
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
