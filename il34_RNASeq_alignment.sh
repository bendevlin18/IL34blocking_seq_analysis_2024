#!/bin/bash
#SBATCH -n 4                    # Number of cores requested
#SBATCH --mem-per-cpu=64G       # 32 GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file

module load GCC/7.4.0
module load Python/2.7.11
module load samtools/1.3.1
module load Bedtools/2.27.1
module load STAR/2.7.5c
module load Subread/1.6.3

### setting the working directories on the cluster

######################################
######################################
######################################
### CHANGE THESE FOR EVERY PROJECT ###
### CHANGE THESE FOR EVERY PROJECT ###
### CHANGE THESE FOR EVERY PROJECT ###
######################################
######################################
######################################

## where all the fastq files are stored. NOTHING ELSE
datadir=/hpc/group/bilbolab/bad36/il34_RNASeq/FASTQ/

## make sure THIS script is saved in the home directory specified here, and run from that directory
homedir=/hpc/group/bilbolab/bad36/il34_RNASeq

## NEED TO ASK ABOUT THIS
star_reference=/hpc/group/bilbolab/bad36/DEPMS_RNASeq/ensembl/genome_dir

## specify path all the way to .gtf including the gunzipped file
gene_gtf=/hpc/group/bilbolab/bad36/DEPMS_RNASeq/ensembl/Mus_musculus.GRCm38.102.chr.gtf

ls $datadir > $homedir/allFiles.txt
awk '{gsub(".fastq", "")}1' $homedir/allFiles.txt > $homedir/fileList.txt
#awk '{gsub(/R./, "")}1' $homedir/fileList.txt > $homedir/uniqSamples.txt

##########################
# Align reads using STAR #
##########################

#mkdir $homedir/aligned/
dir=$datadir/

echo 'Starting Alignment: '$(date +%x_%r)

for i in `cat uniqSamples.txt`
do
	echo 'Aligning '$i' to mus musculus'

    STAR --runMode alignReads --genomeDir $star_reference  --readFilesIn $dir$i'R1.fastq' $dir$i'R2.fastq' --outFileNamePrefix $i --runThreadN 3 --outSAMtype BAM SortedByCoordinate --clip3pNbases 0 --clip5pNbases 3
    samtools index $homedir/$i'mapped.bam'
done

echo 'Finished Alignment: '$(date +%x_%r)

##################################
# Count reads with featureCounts #
##################################
mkdir $homedir/counts/

ls $homedir/bam_files/*.bam > bamFiles.txt

echo 'Starting Counting: '$(date +%x_%r)

for i in `cat bamFiles.txt`
do
featureCounts 	-p -F GTF -a $gene_gtf -s 2 -t gene -g gene_id -o $i'_counts_notrim.txt' $i
done

echo 'Finished Counting. Pipeline Complete: '$(date +%x_%r)
echo 'Enjoy your data!!'

