#! /usr/bin/env bash

#BSUB -J cellranger
#BSUB -o logs/cellranger_%J.out
#BSUB -e logs/cellranger_%J.err
#BSUB -R "select[mem>4] rusage[mem=4]"
#BSUB -q rna

set -o nounset -o pipefail -o errexit -x

module load cellranger/6.0.1

mkdir -p logs

run_snakemake() {
    local config_file=$1
    
    drmaa_args='
        -o {log.out}
        -e {log.err}
        -J {params.job_name} 
        -R "{params.memory} span[hosts=1]"
        -n {threads} '

    snakemake \
        --snakefile Snakefile \
        --drmaa "$drmaa_args" \
        --jobs 200 \
        --latency-wait 60 \
        --configfile $config_file
}

# Dual tags run #1
run_snakemake src/configs/2022-03-11.yaml

# Dual tags run #2
run_snakemake src/configs/2022-10-28.yaml

# Original naive data
run_snakemake src/configs/2019-12-10.yaml


