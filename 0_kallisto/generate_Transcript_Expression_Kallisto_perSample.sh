#!/bin/bash

sample=$1
GenomeIndexFile=$2
ResultsFolder=$3
FASTQFolder=$4
SingleOrPaired="Single"
Strandness="R"
fragmlen=200
fragmlensd=20

if [[ ${SingleOrPaired} == "Single" ]]; then
	fastqfile=$(ls ${FASTQFolder}/${sample}*fastq.gz 2> ~/null)
	if [[ ${fastqfile} != "" ]]; then
		if [[ ${Strandness} == "R" ]]; then
			kallisto/kallisto quant -i ${GenomeIndexFile} -o ${ResultsFolder}/${sample} --single --rf-stranded -l ${fragmlen} -s ${fragmlensd} ${fastqfile}
		elif [[ ${Strandness} == "F" ]]; then
			kallisto/kallisto quant -i ${GenomeIndexFile} -o ${ResultsFolder}/${sample} --single --fr-stranded -l ${fragmlen} -s ${fragmlensd} ${fastqfile}
		else
			echo "Strandness must be either R or F for Single-end reads"
			exit 0
		fi
	fi
else
	echo "Script not adapted to paired-end data"
fi

