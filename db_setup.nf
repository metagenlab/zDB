

process download_cog_cdd {

    output:
        tuple path("COG*.smp"), path("Cog.pn")

    script:
    """
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
    tar xvf cdd.tar.gz && rm cdd.tar.gz
    """
}


process setup_cog_cdd {
    container "$params.blast_container"
    conda "$baseDir/conda/blast.yaml"

    publishDir "$params.cog_db", mode: "move"

    input:
        tuple path(smps), path(pn)

    output:
        file "cog_db*"
        file "cdd_to_cog"

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
    publishDir "$params.pfam_db"

    output:
        tuple path("Pfam-A.hmm"), path("Pfam-A.hmm.dat")

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
    gunzip < Pfam-A.hmm.gz > Pfam-A.hmm && rm -f Pfam-A.hmm.gz
    gunzip < Pfam-A.hmm.dat.gz > Pfam-A.hmm.dat && rm -f Pfam-A.hmm.dat.gz
    """
}


process prepare_hmm {
    container "$params.pfam_scan_container"
    conda "$baseDir/conda/pfam_scan.yaml"

    publishDir "$params.pfam_db", mode: "copy"

    input:
        tuple path(pfam_hmm), path(pfam_defs)

    output:
        path "${pfam_hmm}.h3*"
        path "Pfam-A.hmm.dat"
        path "Pfam-A.hmm"

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

    script:
    version="2022-03-01"
    """
    wget https://www.genome.jp/ftp/db/kofam/archives/$version/ko_list.gz
    wget https://www.genome.jp/ftp/db/kofam/archives/$version/profiles.tar.gz

    gunzip < ko_list.gz > ko_list && rm -f ko_list.gz
    tar xvf profiles.tar.gz
    rm profiles.tar.gz
    """
}


process download_swissprot {
    // not optimal... might be a bit slow to move
    publishDir "$params.swissprot_db", mode: "copy"

    output:
        path "swissprot.fasta"

    script:
    """
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip < uniprot_sprot.fasta.gz > swissprot.fasta
    rm uniprot_sprot.fasta.gz
    """
}


process prepare_swissprot {
    container "$params.blast_container"
    conda "$baseDir/conda/blast.yaml"

    publishDir "$params.swissprot_db", mode: "move"

    input:
        path swissprot_fasta

    output:
        file "*"

    script:
    """
    makeblastdb -dbtype prot -in swissprot.fasta
    rm $swissprot_fasta
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

workflow {
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
}
