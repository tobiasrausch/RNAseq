SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .check .strelka
PBASE=$(shell pwd)

all:   	$(TARGETS)

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(shell uname)-$(shell uname -m).sh" && bash Mambaforge-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Mambaforge-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y -c conda-forge -c bioconda samtools bcftools bedtools htslib bwa alfred star fastqc delly shapeit4 && touch .tools

.strelka: .mamba
	 wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2 && tar xvjf strelka-2.9.10.centos6_x86_64.tar.bz2 && rm strelka-2.9.10.centos6_x86_64.tar.bz2 && touch .strelka

.check: .mamba .tools .strelka
	export PATH=${PBASE}/mamba/bin:${PATH} && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && delly --version && echo "Installation complete!" && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) mamba/ strelka-2.9.10.centos6_x86_64/
