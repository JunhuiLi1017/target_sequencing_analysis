#bsub -W 8:00 -q long 'source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake; bash submit.target.sh &> submit.target.log'
source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake
snakemake -s Snakefile.hg38 --configfile ../config/config.hg38.yaml -p -j 99 --latency-wait 500 --cluster 'bsub -q long -o target_hg38.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 8:00'
conda deactivate
