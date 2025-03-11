
process setup_base_db {
    container "$params.annotation_container"
    conda "$baseDir/conda/annotation.yaml"

    publishDir "$params.zdb_base_db", mode: "move"

    output:
        path("zdb_base")

    script:
    """
    python "$params.setup_base_db_script" --load_cog --load_kegg
    """
}


process download_cog_cdd {

    output:
        tuple path("COG*.smp"), path("Cog.pn"), path("cdd.info")

    script:
    """
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.info
    tar xvf cdd.tar.gz && rm cdd.tar.gz
    """
}


process setup_cog_cdd {
    container "$params.blast_container"
    conda "$baseDir/conda/blast.yaml"

    publishDir "$params.cog_db", mode: "move"

    input:
        tuple path(smps), path(pn), path(info)

    output:
        file "cog_db*"
        file "cdd_to_cog"
        path("${info}", includeInputs: true)

    script:
    """
    makeprofiledb -title COG -in Cog.pn -out cog_db -threshold 9.82 \
        -scale 100.0 -dbtype rps -index true
    grep "tag id" COG* | sed 's/.smp:.*tag id//' | \
                        sed 's/COG//' > cdd_to_cog
    """
}

process diamond_refseq {
    publishDir "$params.refseq_db", mode: "move"
    container "$params.diamond_container"

    input:
        path nr_refseq

    output:
        path "refseq_nr.dmnd"

    script:
    """
        diamond makedb --in $nr_refseq -d refseq_nr
    """
}


process download_pfam_db {

    output:
        tuple path("Pfam-A.hmm"), path("Pfam-A.hmm.dat"), path("Pfam.version")

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz
    gunzip < Pfam-A.hmm.gz > Pfam-A.hmm && rm -f Pfam-A.hmm.gz
    gunzip < Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat && rm -f Pfam-A.hmm.dat.gz
    gunzip < Pfam.version.gz > Pfam.version && rm -f Pfam.version.gz
    """
}


process prepare_hmm {
    container "$params.pfam_scan_container"
    conda "$baseDir/conda/pfam_scan.yaml"

    publishDir "$params.pfam_db", mode: "move"

    input:
        tuple path(pfam_hmm), path(pfam_defs), path(pfam_version)

    output:
        path "${pfam_hmm}.h3*"
        path("${pfam_hmm}", includeInputs: true)
        path("${pfam_defs}", includeInputs: true)
        path("${pfam_version}", includeInputs: true)

    script:
    """
        hmmpress $pfam_hmm
    """
}


process download_ko_profiles {
    publishDir "$params.ko_db", mode: "move"

    output:
        path "ko_list"
        path "profiles/*.hmm"
        path "profiles/prokaryote.hal"
        path "version.txt"

    script:
    version="2022-03-01"
    """
    wget https://www.genome.jp/ftp/db/kofam/archives/$version/ko_list.gz
    wget https://www.genome.jp/ftp/db/kofam/archives/$version/profiles.tar.gz

    gunzip < ko_list.gz > ko_list && rm -f ko_list.gz
    tar xvf profiles.tar.gz
    rm profiles.tar.gz
    echo $version > version.txt
    """
}


process download_swissprot {

    output:
        tuple path("swissprot.fasta"), path("relnotes.txt")

    script:
    """
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt
    gunzip < uniprot_sprot.fasta.gz > swissprot.fasta
    rm uniprot_sprot.fasta.gz
    """
}


process prepare_swissprot {
    container "$params.blast_container"
    conda "$baseDir/conda/blast.yaml"

    publishDir "$params.swissprot_db", mode: "move"

    input:
        tuple path(swissprot_fasta), path(relnotes)

    output:
        path("*", includeInputs: true)

    script:
    """
    makeblastdb -dbtype prot -in swissprot.fasta
    """
}

process download_vfdb {

    output:
        tuple path("vfdb.fasta"), path("VFs.xls")

    script:
    """
    wget http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
    wget http://www.mgc.ac.cn/VFs/Down/VFs.xls.gz
    gunzip < VFDB_setB_pro.fas.gz > vfdb.fasta
    gunzip < VFs.xls.gz > VFs.xls
    rm VFDB_setB_pro.fas.gz
    rm VFs.xls.gz
    """
}

process prepare_vfdb {
    container "$params.blast_container"
    conda "$baseDir/conda/blast.yaml"

    publishDir "$params.vf_db", mode: "move"

    input:
        tuple path(vfdb_fasta), path(vf_descr)

    output:
        path("*", includeInputs: true)

    script:
    """
    makeblastdb -dbtype prot -in $vfdb_fasta
    """
}

workflow setup_cogg_db {
    download_cog_cdd() | setup_cog_cdd
}

/**
 * NOTE: this assumes that the user already downloaded the refseq
 * database (nr_refseq) and concatenated all the files into a single
 * fasta file called "refseq_nr.fasta".
 */
workflow setup_refseq_db {
    Channel.fromPath("${params.refseq_db}/refseq_nr.fasta") | diamond_refseq
}

workflow setup_pfam_db {
    download_pfam_db() | prepare_hmm
}

workflow setup_ko_db {
    download_ko_profiles()
}

workflow setup_swissprot_db {
    download_swissprot() | prepare_swissprot
}

workflow setup_vfdb {
    download_vfdb() | prepare_vfdb
}

workflow {
    if( params.setup_base_db )
        setup_base_db()

    if( params.cog )
        setup_cogg_db()

    if ( params.diamond_refseq )
        setup_refseq_db()

    if ( params.pfam )
        setup_pfam_db()

    if ( params.ko )
        setup_ko_db()

    if ( params.blast_swissprot )
        setup_swissprot_db()

    if ( params.vfdb )
        setup_vfdb()
}
