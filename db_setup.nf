
process setup_base_db {
    label 'mount_basedir'
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

process download_refseq {
    publishDir "$params.refseq_db", mode: "copy"

    output:
        tuple path("refseq_nr.fasta"), path("RELEASE_NUMBER")

    script:
    refseq_version_link="https://ftp.ncbi.nlm.nih.gov/refseq/release/RELEASE_NUMBER"
    refseq_complete_base_link="https://ftp.ncbi.nlm.nih.gov/refseq/release/complete"
    refseq_html_file="refseq.html"
    refseq_filenames_file="refseq-genome-files.txt"
    """
    # Adapted from script from Mike Lee https://hackmd.io/@AstrobioMike/RefSeq-db-download

    # Download the version file
    curl -L -s -O "${refseq_version_link}"

    # downloading html page (using this to get all the files we want to download)
    curl -L -s -o ${refseq_html_file} ${refseq_complete_base_link}

    # parsing out genomic.fna.gz filenames (which are also their link suffixes)
    grep "complete.wp_protein." ${refseq_html_file} | grep "protein.faa.gz" | cut -f 2 -d '"' > ${refseq_filenames_file}

    num_files=\$(wc -l ${refseq_filenames_file})

    printf "\n  We are beginning the download of \${num_files} files now...\n"
    printf "  See you in a bit :)\n\n"

    # downloading in parallel with xargs (num run in parallel is set with -P option)
    xargs -I % -P 10 curl -L -s -O "${refseq_complete_base_link}/%" < ${refseq_filenames_file}

    # Merge the files and delete them as we go, to avoid using too much disk space
    touch refseq_nr.fasta.gz
    while read -r line; do
        echo "Adding and removing: \$line"
        cat \$line >> refseq_nr.fasta.gz
        rm \$line
    done < ${refseq_filenames_file}

    gunzip refseq_nr.fasta.gz
    """
}


process diamond_refseq {
    publishDir "$params.refseq_db", mode: "move"
    container "$params.diamond_container"

    input:
        tuple path(nr_refseq), path(nr_version)

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
    """
    wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
    wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz

    gunzip < ko_list.gz > ko_list && rm -f ko_list.gz
    tar xvf profiles.tar.gz
    rm profiles.tar.gz
    echo "\$(date +'%Y-%m-%d')" > version.txt
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

workflow setup_refseq_db {
    download_refseq() | diamond_refseq
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
