#bsub -W 4:00 -q short 'source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake; bash submit.target.sh &> submit.target.log'
source ~/anaconda3/etc/profile.d/conda.sh; conda activate snakemake
snakemake -s /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_seq/workflow/Snakefile --configfile /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_seq/config/config.yaml -p -j 99 --latency-wait 50 --cluster 'bsub -q short -o /pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_seq/workflow/target.log -R "rusage[mem=4000]" -n 8 -R span[hosts=1] -W 4:00'
conda deactivate
