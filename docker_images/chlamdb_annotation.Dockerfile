# gblast3 OK
# plast OK
# PRIAM OK
# kofamscan OK

# BPBAac ok
# T3_MM ok
# DeepT3
# effectiveT3


FROM continuumio/miniconda3:4.6.14

ENV TZ Europe/Zurich

RUN echo $TZ > /etc/timezone && \
    apt-get update && apt-get install -y tzdata && \
    rm /etc/localtime && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    apt-get clean

RUN conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda

RUN apt-get update && apt-get -yq install curl \
	wget \
	build-essential \
    git \ 
    procps \ 
    ncbi-blast+-legacy && \ 
    apt-get clean

WORKDIR /usr/local/bin 

# export "PATH=\$HMMTOP_PATH:\$GBLAST3_PATH:\$PATH"
RUN git clone -b dev --single-branch https://github.com/metagenlab/annotation_pipeline_nextflow.git && echo ok

ENV PATH=/usr/local/bin/annotation_pipeline_nextflow/bin/hmmtop_2.1/:$PATH
ENV HMMTOP_ARCH=/usr/local/bin/annotation_pipeline_nextflow/bin/hmmtop_2.1/hmmtop.arch
ENV HMMTOP_PSV=/usr/local/bin/annotation_pipeline_nextflow/bin/hmmtop_2.1/hmmtop.psv

RUN mkdir -p /usr/local/bin/tcdb_db/

RUN curl -S http://www.tcdb.org/public/tcdb > /usr/local/bin/tcdb_db/tcdb
RUN curl -S http://tcdb.org/public/betabarrel > /usr/local/bin/tcdb_db/betabarrel

RUN formatdb -i /usr/local/bin/tcdb_db/tcdb -p T
RUN formatdb -i /usr/local/bin/tcdb_db/betabarrel -p T

RUN git clone https://github.com/metagenlab/BioVx.git

ENV PATH=/usr/local/bin/BioVx/scripts:$PATH

RUN wget http://plast.gforge.inria.fr/files/plastbinary_linux_v2.3.1.tar.gz && tar zxvf plastbinary_linux_v2.3.1.tar.gz \
&& mv plastbinary_linux_20160121/build/bin/plast . && rm -rf plastbinary_linux*

RUN conda env create -f /usr/local/bin/annotation_pipeline_nextflow/docker_images/gblast.yml

RUN mkdir -p /usr/local/bin/PRIAM/

RUN wget http://priam.prabi.fr/utilities/PRIAM_search.jar

RUN conda env create -f /usr/local/bin/annotation_pipeline_nextflow/docker_images/priam.yml

RUN wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan.tar.gz && tar zxvf plastbinary_linux_v2.3.1.tar.gz && rm plastbinary_linux_v2.3.1.tar.gz

RUN wget https://biocomputer.bio.cuhk.edu.hk/softwares/T3_MM/T3_MM.tar.gz && tar zxvf T3_MM.tar.gz && rm T3_MM.tar.gz

RUN wget https://biocomputer.bio.cuhk.edu.hk/softwares/BPBAac/BPBAac.tar.gz && tar zxvf BPBAac.tar.gz && rm BPBAac.tar.gz

RUN conda env create -f /usr/local/bin/annotation_pipeline_nextflow/docker_images/BPBAac.yml

RUN conda env create -f /usr/local/bin/annotation_pipeline_nextflow/docker_images/deep_t3.yml

RUN conda env create -f /usr/local/bin/annotation_pipeline_nextflow/docker_images/effective.yml

RUN mkdir effective 

WORKDIR /usr/local/bin/effective

RUN mkdir module

RUN wget http://effectors.csb.univie.ac.at/sites/eff/files/others/TTSS_GUI-1.0.1.jar

RUN curl -L http://effectors.csb.univie.ac.at/sites/eff/files/others/TTSS_STD-2.0.2.jar > /usr/local/bin/effective/module/TTSS_STD-2.0.2.jar

WORKDIR /usr/local/bin

RUN git clone https://github.com/lje00006/DeepT3.git

ENV PYTHONPATH=/usr/local/bin/DeepT3/DeepT3-Keras:$PYTHONPATH

RUN conda clean --all --yes

ENTRYPOINT ["/bin/bash"]
CMD ["/bin/bash"]