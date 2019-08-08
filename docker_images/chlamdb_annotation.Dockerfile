# gblast3 OK
# plast OK
# PRIAM OK
# kofamscan OK

# macsyfinder
# genome-properties

# BPBAac ok
# T3_MM ok
# DeepT3 ok
# effectiveT3 ok


FROM ubuntu:18.04

ENV TZ Europe/Zurich

RUN echo $TZ > /etc/timezone && \
    apt-get update && apt-get install -y tzdata && \
    rm /etc/localtime && \
    ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata && \
    apt-get clean

## Enable Ubuntu Universe, Multiverse, and deb-src for main.
RUN sed -i 's/^#\s*\(deb.*main restricted\)$/\1/g' /etc/apt/sources.list
RUN sed -i 's/^#\s*\(deb.*universe\)$/\1/g' /etc/apt/sources.list
RUN sed -i 's/^#\s*\(deb.*multiverse\)$/\1/g' /etc/apt/sources.list

RUN apt-get update && apt-get -yq install curl \
	wget \
    git \
    build-essential \
    procps \
    ncbi-blast+-legacy \
    python3.6 \
    python3-distutils \
    ncbi-blast+ \
    ncbi-blast+-legacy \
    r-base \
    r-cran-e1071 \
    python-dev \
    openjdk-8-jdk \
    perl && apt-get clean

WORKDIR /usr/local/bin

RUN wget https://github.com/wrpearson/fasta36/releases/download/fasta-v36.3.8g/fasta-36.3.8g-linux64.tar.gz && tar -xvzf fasta-36.3.8g-linux64.tar.gz && mv fasta-36.3.8g/bin/* . && rm -rf fasta-36.3.8g*

RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python2.7 get-pip.py
RUN python3.6 get-pip.py

RUN pip2 install numpy
RUN pip2 install biopython==1.67 matplotlib==2.2.3

RUN pip3 install tensorflow keras==2.2.4 biopython==1.73

 #- fasta
 #- blast=2.7.1
 #- openjdk=8.0.152
 #- blast-legacy=2.2.26
 #- r-base
 #- r-e1071
 #- keras=2.2.4
 #- biopython=1.73

RUN git clone -b dev --single-branch https://github.com/metagenlab/annotation_pipeline_nextflow.git

ENV PATH "/usr/local/bin/annotation_pipeline_nextflow/bin/hmmtop_2.1/:${PATH}"
ENV HMMTOP_ARCH "/usr/local/bin/annotation_pipeline_nextflow/bin/hmmtop_2.1/hmmtop.arch"
ENV HMMTOP_PSV "/usr/local/bin/annotation_pipeline_nextflow/bin/hmmtop_2.1/hmmtop.psv"

RUN mkdir -p /usr/local/bin/tcdb_db/

RUN curl -S http://www.tcdb.org/public/tcdb > /usr/local/bin/tcdb_db/tcdb
RUN curl -S http://tcdb.org/public/betabarrel > /usr/local/bin/tcdb_db/betabarrel

RUN formatdb -i /usr/local/bin/tcdb_db/tcdb -p T
RUN formatdb -i /usr/local/bin/tcdb_db/betabarrel -p T

RUN git clone https://github.com/metagenlab/BioVx.git

ENV PATH "/usr/local/bin/BioVx/scripts:${PATH}"

RUN wget http://plast.gforge.inria.fr/files/plastbinary_linux_v2.3.1.tar.gz && tar zxvf plastbinary_linux_v2.3.1.tar.gz \
&& mv plastbinary_linux_20160121/build/bin/plast . && rm -rf plastbinary_linux*

RUN mkdir -p /usr/local/bin/PRIAM/

RUN wget http://priam.prabi.fr/utilities/PRIAM_search.jar

RUN git -C /usr/local/bin/annotation_pipeline_nextflow pull

RUN wget ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan.tar.gz && tar zxvf kofamscan.tar.gz && rm kofamscan.tar.gz

RUN wget https://biocomputer.bio.cuhk.edu.hk/softwares/T3_MM/T3_MM.tar.gz && tar zxvf T3_MM.tar.gz && rm T3_MM.tar.gz

RUN wget https://biocomputer.bio.cuhk.edu.hk/softwares/BPBAac/BPBAac.tar.gz && tar zxvf BPBAac.tar.gz && rm BPBAac.tar.gz

RUN mkdir effective

WORKDIR /usr/local/bin/effective

RUN mkdir module

RUN wget http://effectors.csb.univie.ac.at/sites/eff/files/others/TTSS_GUI-1.0.1.jar

RUN curl -L http://effectors.csb.univie.ac.at/sites/eff/files/others/TTSS_STD-2.0.2.jar > /usr/local/bin/effective/module/TTSS_STD-2.0.2.jar

WORKDIR /usr/local/bin

RUN git clone https://github.com/lje00006/DeepT3.git

ENV PYTHONPATH "/usr/local/bin/DeepT3/DeepT3/DeepT3-Keras:${PYTHONPATH}"

RUN wget https://github.com/gem-pasteur/macsyfinder/archive/macsyfinder-1.0.5.tar.gz && tar zxvf macsyfinder-1.0.5.tar.gz && mv macsyfinder-macsyfinder-1.0.5 macsyfinder

WORKDIR /usr/local/bin/macsyfinder

RUN python setup.py build

RUN python setup.py test -vv

ENV MACSY_HOME "/usr/local/bin/macsyfinder/"

RUN apt install -y hmmer && apt-get clean

WORKDIR /usr/local/bin

ENTRYPOINT ["/bin/bash"]
#CMD ["/bin/bash"]
