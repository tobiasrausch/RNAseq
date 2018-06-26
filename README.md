RNA-Seq analysis workflow
-------------------------

Installation

`git clone https://github.com/tobiasrausch/RNAseq.git`

`cd RNAseq`

`make all`

If one of the above commands fail your operating system probably lacks some build essentials. These are usually pre-installed but if you lack them you need to install these. For instance, for Ubuntu this would require:

`apt-get install build-essential g++ git wget unzip`


Prerequisites
-------------

Download STAR GRCh37.75 genome

`cd genome/ && ./starGenome.sh`

Download the GTF annotation file

`cd gtf/ && ./downloadGTF.sh`

Download the 1000 Genomes Reference Panel for allele-specific expression

`cd refpanel/ && ./download_1kGP_hg19.sh`


RNA-Seq alignment, gene counting and variant calling
----------------------------------------------------

To perform alignment, generate browser tracks, and call variants you can use the RNA wrapper script.

`./src/rna.sh ERR188040_1.fastq.gz ERR188040_2.fastq.gz HG00096`

This script also phases variants against the 1000 Genomes reference panel and generates an allele-specific count table.


Example Workflow
----------------

To run an entire sample from the gEUVADIS data set

`cd example/ && ./gEUVADIS.sh`

