// Note: we will probably need to use a fixed version
// of the different databases to avoid bugs due to 
// changes in formatting

/*
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


if(!params.cog_location && params.cog) {

    process download_cog_cdd {
        output:
            file "COG*.smp" into to_make_cdd_profiles
            file "Cog.pn" into to_make_cdd_profiles_pn
            
        script:
        """
        wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
        tar xvf cdd.tar.gz && rm cdd.tar.gz
        """
    }
} else if(params.cog) {
    Channel.fromPath(params.cog_location+"/COG*.smp")
        .collect().set { to_make_cdd_profiles }
    Channel.fromPath(params.cog_location+"/Cog.pn").set { to_make_cdd_profiles_pn }
}


if(params.cog) {
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
}


process download_pfam_db {
    publishDir "$params.pfam_db"
    when:
        params.pfam_scan
    
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
        file "${pfam_hmm}.h3*"
    
    script:
    """
        hmmpress $pfam_hmm
    """
}


if(!params.ko_location && params.ko) {
    process download_KO_data {
        publishDir "$params.ko_data"
        container "$params.annotation_container"

        output:
            file "ko*" into to_load_ko

        script:
        """
        /home/metagenlab/metagenlab_libs/chlamdb/chlamdb.py \
            --download_ko_files --ko_dir=${params.ko_db}
        """
    }

    process download_ko_profiles {
        publishDir "$params.ko_data", mode: "move"

        output:
            file "ko_list"
            file "profiles/*.hmm"
    
        script:
        """
        wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
        wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz

        gunzip < ko_list.gz > ko_list && rm -f ko_list.gz
        tar xvf profiles.tar.gz
        """
    }
} else if (params.ko) {
    Channel.fromPath(params.ko_location + "ko:*").collect()
        .set { to_load_ko }
} else {
    Channel.empty().set { to_load_ko }
}


if(!params.swissprot_location && params.blast_swissprot) {
    process download_swissprot_db {
        when:
            params.blast_swissprot

        output:
            file "swissprot.fasta" into swissprot_fasta

        script:
        """
        wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip < uniprot_sprot.fasta.gz > swissprot.fasta && rm -f uniprot_sprot.fasta
        """
    }
} else if(params.blast_swissprot) {
    Channel.fromPath(params.swissprot_location+"/swissprot.fasta")
        .set { swissprot_fasta }
} else {
    Channel.empty().set { swissprot_fasta }
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

to_setup_db_cog.mix(to_load_ko).set { to_setup_db }


process setup_base_db {
    container "$params.annotation_container"
    publishDir "$params.base_db", mode: "move"

    input:
        file setup from to_setup_db.collect()

    output:
        file "George"

    script:
    command_line_params = ""
    command_line_params += (params.cog)? "--load_cog " : ""
    command_line_params += (params.cog)? "--cog_dir ${params.cog_db} " : ""
    command_line_params += (params.ko)? "--load_kegg " : ""
    command_line_params += (params.ko)? "--ko_dir=${params.ko_location}" : ""
    """
        /home/metagenlab/metagenlab_libs/chlamdb/chlamdb.py \
            --biosql_schema /home/metagenlab/metagenlab_libs/biosql_schema/biosqldb-sqlite.sql \
            ${command_line_params}
    """
}
