#!/bin/bash

# Files
cDNAFile="Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
ncRNAFile="Arabidopsis_thaliana.TAIR10.ncrna.fa.gz"
GenomeIndexFile="TAIR10_Index_cDNAs_ncRNAs"

# Indexing genome
if [[ -s ${cDNAFile} ]] &&  [[ -s ${ncRNAFile} ]] ; then
	kallisto/kallisto index -i ${GenomeIndexFile} ${cDNAFile} ${ncRNAFile}
else
	echo "${cDNAFile} and/or ${ncRNAFile} not present or empty"
	exit 0
fi

