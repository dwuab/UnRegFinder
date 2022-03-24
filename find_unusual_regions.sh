# An example script to run the pipeline
# It can be run on a single computing node or computer clusters using PBS/Torque, Slurm, etc.

set -euox

pipeline_dir=/data1/unusual_region_code/

# Run on a single computing node
snakemake -s ${pipeline_dir}/Snakefile -j 5 all -p --configfile config.json

# Run on PBS/Torque
# snakemake -s ${pipeline_dir}/Snakefile -j 5 all -p --configfile config.json \
# 	--cluster "qsub -q cu -N {params.job_name} -l mem={params.mem} -j oe -o logs/" \
# 	--latency-wait 120
