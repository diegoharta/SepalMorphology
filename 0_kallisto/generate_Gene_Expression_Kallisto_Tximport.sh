#!/bin/bash
export PATH=/usr/local/R-3.6.0/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/snap/bin

# Scripts
PerSampleKallistoScript="generate_Transcript_Expression_Kallisto_perSample.sh"
TximportRScript="transcript_to_gene_tximport.R"

# Files & parameters
samples=$1
GenomeIndexFile="TAIR10_Index_cDNAs_ncRNAs"
ResultsFolder="Results"
mkdir -p ${ResultsFolder}
tmpMetafile=${ResultsFolder}/Metafile.tmp
echo "ID" > ${tmpMetafile}
FASTQFolder="FASTQ"
AnnotationsTrasncToGene="transcript_to_gene_EnsemblPlants47.txt.gz"

kallistodone=1
# Kallisto
while IFS= read -r sample
do
	echo ${sample} >> ${tmpMetafile}
	if [[ ! -s ${ResultsFolder}/${sample}/abundance.tsv ]]
	then
		kallistodone=0
	fi
done < "$samples"

if [[ ${kallistodone} -eq 0 ]] 
then
	echo "Running kallisto"
	for sample in $(tail -n +2 ${tmpMetafile})
	do
		echo "${sample}"
		mkdir -p ${ResultsFolder}/${sample}
		rm ${ResultsFolder}/${sample}/abundance.* 2> ~/null
		fastqpresent=$(ls ${FASTQFolder}/${sample}*fastq.gz 2> ~/null | wc -l)
		if [[ ${fastqpresent} -gt 0 ]]
		then
			echo "Running kallisto script"
			./${PerSampleKallistoScript} ${sample} ${GenomeIndexFile} ${ResultsFolder} ${FASTQFolder}
		else
			echo "${FASTQFolder}/${sample}*.fastq.gz missing"
			exit 0
		fi
	done
else
	echo "kallisto already done"
fi

# Tximport
if [[ ! -s ${ResultsFolder}/Counts_kallisto_tximport.tab ]] && [[ ! -s ${ResultsFolder}/TPM_kallisto_tximport.tab ]]
then
	echo "Running tximport"
	./${TximportRScript} ${tmpMetafile} ${AnnotationsTrasncToGene} ${ResultsFolder}
else
	echo "tximport already done"
fi






















