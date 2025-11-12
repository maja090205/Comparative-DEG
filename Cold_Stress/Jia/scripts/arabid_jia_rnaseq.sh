#!/bin/bash

#SBATCH -J arabid_jia_rnaseq
#SBATCH --time=6-00:00
#SBATCH --chdir=/work/zarski/arabid_jia/skripts/
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=maja-katharina.zarski.idiv@ufz.de
#SBATCH --mail-type=BEGIN,END,FAIL

#output files
#SBATCH -o /work/%u/job_logs/%x-%j.out
#SBATCH -e /work/%u/job_logs/%x-%j.err

module load Nextflow/24.10.1

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
nextflow run nf-core/rnaseq \
    --input samplesheet_rnaseq.csv \
    --fasta /work/zarski/for_hiwi/arabid_genome/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
    --gtf /work/zarski/for_hiwi/arabid_genome/Arabidopsis_thaliana.TAIR10.60.gtf.gz \
    -profile singularity \
    --outdir /work/zarski/arabid_jia