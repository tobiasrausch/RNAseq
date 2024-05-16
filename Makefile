SHELL := /bin/bash

# Targets
TARGETS = .mamba .tools .check
PBASE=$(shell pwd)

all:   	$(TARGETS)

.mamba:
	curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(shell uname)-$(shell uname -m).sh" && bash Mambaforge-$(shell uname)-$(shell uname -m).sh -b -p mamba && rm "Mambaforge-$(shell uname)-$(shell uname -m).sh" && touch .mamba

.tools: .mamba
	export PATH=${PBASE}/mamba/bin:${PATH} && mamba install -y -c conda-forge -c bioconda samtools bcftools bedtools htslib bwa igvtools alfred star fastqc freebayes delly shapeit4 && touch .tools

.check: .mamba .tools
	export PATH=${PBASE}/mamba/bin:${PATH} && samtools --version && bcftools --version && bedtools --version && bgzip --version && tabix --version && delly --version && echo "Installation complete!" && touch .check

clean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/
