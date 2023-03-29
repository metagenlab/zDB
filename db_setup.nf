

process download_cog_cdd {
    when:
        params.cog
        
    output:
        file "COG*.smp" into to_make_cdd_profiles
        file "Cog.pn" into to_make_cdd_profiles_pn

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
        file smps from to_make_cdd_profiles
        file pn from to_make_cdd_profiles_pn

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


Channel.fromPath("${params.refseq_db}/refseq_nr.fasta"). set { nr_refseq  }


/**
 * NOTE: this assumes that the user already downloaded the refseq 
 * database (nr_refseq) and concatenated all the files into a single
 * fasta file called "refseq_nr.fasta".
 */
process diamond_refseq {
    publishDir "$params.refseq_db", mode: "move"
    container "$params.diamond_container"

    when:
        params.diamond_refseq

    input:
        file nr_refseq

    output:
        file "refseq_nr.dmnd"
        
    script:
    """
        diamond makedb --in $nr_refseq -d refseq_nr
    """
}


process download_pfam_db {
    publishDir "$params.pfam_db"
    when:
        params.pfam
    
    output:
        file "Pfam-A.hmm" into pfam_hmm
        file "Pfam-A.hmm.dat" into pfam_defs

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
        file pfam_hmm
        file pfam_defs

    output:
        file "${pfam_hmm}.h3*"
        file "Pfam-A.hmm.dat"
        file "Pfam-A.hmm"
    
    script:
    """
        hmmpress $pfam_hmm
    """
}


process download_ko_profiles {
    publishDir "$params.ko_db", mode: "move"

    when:
        params.ko

    output:
        file "ko_list"
        file "profiles/*.hmm"
        file "profiles/prokaryote.hal"

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

    when:
        params.blast_swissprot

    output:
        file "swissprot.fasta" into swissprot_fasta

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
        file swissprot_fasta

    output:
        file "*"

    script:
    """
    makeblastdb -dbtype prot -in swissprot.fasta
    rm $swissprot_fasta
    """
}
