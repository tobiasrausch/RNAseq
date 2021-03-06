SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .envs .eagle
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/bin && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/bin/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.envs: .conda .channels
	export PATH=${PBASE}/bin/bin:${PATH} && conda create -y --prefix=${PBASE}/bin/envs/rna samtools=1.7 igvtools=2.3.93 alfred=0.1.8 fastqc=0.11.7 bcftools=1.7 star=2.6.0c bedtools=2.27.1 freebayes=1.2.0 && conda install -y --prefix=${PBASE}/bin/envs/rna -c conda-forge ncurses=5.9 && touch .envs

.eagle: .conda
	cd ${PBASE}/bin && wget 'https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.tar.gz' && tar -xzf Eagle_v2.4.tar.gz && rm Eagle_v2.4.tar.gz && cd ../ && touch .eagle

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/
