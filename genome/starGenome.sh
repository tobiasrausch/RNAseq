#!/bin/bash

# Download hg19 reference genome
wget -r -nH -nd -np -R index.html* 'http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/ENSEMBL/homo_sapiens/ENSEMBL.homo_sapiens.release-75/'
