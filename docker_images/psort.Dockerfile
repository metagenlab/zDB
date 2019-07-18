
FROM ubuntu
LABEL source="https://github.com/brinkmanlab/psortb-docker/"

# Install packages then remove cache package list information
RUN DEBIAN_FRONTEND=noninteractive apt-get update && apt-get -yq install curl \
	wget \
	build-essential \
	net-tools \
	librpc-xml-perl \
	bioperl \
	ncbi-blast+-legacy \
	nano \
	libf2c2 \
	libjson-rpc-perl
#	fort77 \

WORKDIR /usr/local/src

RUN echo '/usr/local/lib64' >>/etc/ld.so.conf

RUN wget http://www.psort.org/download/docker/pft2.3.4.docker64bit.tar.gz && tar zxvf pft2.3.4.docker64bit.tar.gz && cp pftools/pfscan /usr/local/bin/

RUN wget http://www.psort.org/download/libpsortb-1.0.tar.gz && tar zxvf libpsortb-1.0.tar.gz && cd libpsortb-1.0 && ./configure && make && make install && ldconfig

RUN wget http://www.psort.org/download/bio-tools-psort-all.3.0.6.tar.gz && tar zxvf bio-tools-psort-all.3.0.6.tar.gz

WORKDIR /usr/local/src/bio-tools-psort-all

RUN wget http://www.psort.org/download/docker/psortb.defaults

RUN perl Makefile.PL && make && make install && cp -r psort /usr/local/psortb

RUN formatdb -i /usr/local/psortb/conf/analysis/sclblast/gramneg/sclblast -p T

RUN apt-get clean && apt-get update && apt-get install -y locales
RUN locale-gen en_US.UTF-8
ENV LANG en_US.UTF-8 
ENV LANGUAGE en_us:en  
RUN export LC_ALL=en_US.UTF-8
RUN export LANG=en_US.UTF-8

CMD ["/bin/bash"]