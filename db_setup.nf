// Note: we will probably need to use a fixed version
// of the different databases to avoid bugs due to 
// changes in formatting

/*
process setup_cog_db {
    when:
        params.cog
}


process setup_refseq_db {
    when:
        params.diamond_refseq
}

process setup_ko_db {
    when:
        params.ko
} */


process download_cog_defs {
    publishDir "$params.cog_db"
    when:
        params.cog

    output:
        file "cog-20.def.tab" into cog_defs
        file "fun-20.tab" into cog_funcs
    
    script:
    """
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/cog-20.def.tab
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/fun-20.tab

    iconv -f cp1252 -t utf-8 cog-20.def.tab > cog_translated.def.tab
    mv cog_translated.def.tab cog-20.def.tab
    """
}

cog_defs.mix(cog_funcs).set { to_setup_db_cog }


process download_cog_cdd {
    when:
        params.cog && !params.only_base_db

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
    publishDir "$params.cog_db"

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


process download_pfam_db {
    publishDir "$params.pfam_db"
    when:
        params.pfam_scan && !params.only_base_db
    
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
    publishDir "$params.pfam_db"

    input:
        file pfam_hmm

    output:
        file "$pfam_hmm.h3*"
    
    script:
    """
        hmmpress $pfam_hmm_profile
    """
}


process download_swissprot_db {
    when:
        params.blast_swissprot && !params.only_base_db

    output:
        file "swissprot.fasta" into swissprot_fasta

    script:
    """
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
    gunzip < uniprot_sprot.fasta.gz > swissprot.fasta && rm -f uniprot_sprot.fasta
    """
}


process setup_swissprot {
    container "$params.blast_container"
    publishDir "$params.swissprot_db"

    input:
        file swissprot_fasta

    output:
        file "*"
        file "$swissprot_fasta"

    script:
    """
    makeblastdb -dbtype prot -in $swissprot_fasta
    """
}

to_setup_db_cog.set { to_setup_db }


process setup_base_db {
    container "$params.annotation_container"
    publishDir "$params.base_db"

    input:
        file setup from to_setup_db.collect()

    output:
        file "George"

    script:
    command_line_params = ""
    command_line_params += (params.cog)? "--load_cog " : ""
    command_line_params += (params.cog)? "--cog_dir ${baseDir}/${params.cog_db} " : ""
    """
        /home/metagenlab/metagenlab_libs/chlamdb/chlamdb.py \
            --biosql_schema /home/metagenlab/metagenlab_libs/biosql_schema/biosqldb-sqlite.sql \
            ${command_line_params}
    """
}

