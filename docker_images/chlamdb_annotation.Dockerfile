# gblast3 OK
# plast OK
# PRIAM

# BPBAac
# DeepT3
# effectiveT3
# T3_MM

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

RUN conda env create -f /usr/local/bin/annotation_pipeline_nextflow/docker_images/gblast.yml /gblast3.yml

RUN conda clean --all --yes

ENTRYPOINT ["/bin/bash"]
CMD ["/bin/bash"]