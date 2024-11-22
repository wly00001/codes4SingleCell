#!/bin/bash


#SBATCH --requeue                    # Return job to the queue if preempted
#SBATCH --job-name=sc_CR     # Assign an 8-character name to your job, no spaces, no special characters
#SBATCH --nodes=1                # Number of compute nodes
#SBATCH --ntasks=1              # Number of tasks to run (often = cores) on each node
#SBATCH --cpus-per-task=24
#SBATCH --time=100:00:00
#SBATCH --output=slurm.%x.%N.%j.o
#SBATCH --export=ALL
#SBATCH --mail-user=lwang@wistar.org
#SBATCH --mail-type=ALL

##### code starts here #####
echo Time is `date`, running is started
#biosoft path#
#/wistar/tian/lwang/module/samtools-1.11/samtools 
#/wistar/tian/lwang/module/STAR-2.7.7a/source/STAR
#/wistar/tian/lwang/module/FastQC/fastqc
#/wistar/tian/lwang/module/R-4.0.2/bin/Rscript

# export PATH=$PATH:/wistar/tian/lwang/module/samtools-1.11/ 
# export PATH=$PATH:/wistar/tian/lwang/module/STAR-2.7.7a/source/
# export PATH=$PATH:/wistar/tian/lwang/module/FastQC/
# export PATH=$PATH:/wistar/tian/lwang/module/R-4.0.2/bin/
# export PATH=$PATH:/wistar/tian/lwang/module/TrimGalore-0.6.6/




# /wistar/tian/lwang/module/cellranger-7.0.1/cellranger count --id=PBS \
#    --fastqs=/wistar/tian/lwang/scRNAseq/GSE160450/rawfastq \
#    --sample=SRR12937940 \
#    --transcriptome=/wistar/tian/lwang/reference/CellRanger/refdata-gex-mm10-2020-A



/wistar/tian/lwang/module/cellranger-7.0.1/cellranger count --id=LPS \
   --fastqs=/wistar/tian/lwang/scRNAseq/GSE160450/rawfastq \
   --sample=SRR12937941 \
   --transcriptome=/wistar/tian/lwang/reference/CellRanger/refdata-gex-mm10-2020-A


