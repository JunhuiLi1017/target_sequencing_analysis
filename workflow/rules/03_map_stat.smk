rule bqsr_insert:
	input:
		bam="{outpath}/02_map/05_apply_bqsr/{sample}/{sample}.bam"
	output:
		png="{outpath}/02_map/06_stat/01_insert/01_insert/{sample}.bqsr.insert.png",
		table="{outpath}/02_map/06_stat/01_insert/01_insert/{sample}.bqsr.fragment.txt"
	log:
		"{outpath}/02_map/logs/{sample}.deeptools.log"
	params:
		maxFragmentLength=2000
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image['deeptools']
	shell:
		"""
		source ~/anaconda3/etc/profile.d/conda.sh
		conda activate deeptools
		bamPEFragmentSize -b {input.bam} -o {output.png} --maxFragmentLength 2000 --table {output.table}
		conda deactivate
		"""

rule targt_intersect_bam:
	input:
		bam="{outpath}/02_map/{bam_type}/{sample}/{sample}.{bam_type}.bam"
	output:
		bed=temp("{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.intersect.bed")
	log:
		"{outpath}/02_map/logs/{sample}.{bam_type}.target.log"
	params:
		TargeRegion=config['TargeRegion']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		container_image['bedtools']
	shell:
		"""
		bedtools intersect -abam {input.bam} -b {params.TargeRegion} -wa -bed > {output.bed}
		"""

rule targt_intersectV_bam:
	input:
		bam="{outpath}/02_map/{bam_type}/{sample}/{sample}.{bam_type}.bam"
	output:
		bed=temp("{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.v.intersect.bed")
	log:
		"{outpath}/02_map/logs/{sample}.{bam_type}.target.log"
	params:
		TargeRegion=config['TargeRegion']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/bedtools.sif"
	shell:
		"""
		bedtools intersect -abam {input.bam} -b {params.TargeRegion} -v -bed > {output.bed}
		"""


rule targt_cov_bam:
	input:
		bam="{outpath}/02_map/{bam_type}/{sample}/{sample}.{bam_type}.bam"
	output:
		coverage="{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.coverage"
	params:
		TargeRegion=config['TargeRegion'],
		sample="{sample}",
		versionsorted=config['versionsorted']
	log:
		"{outpath}/02_map/logs/{sample}.cov.target.{bam_type}.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/bedtools.sif"
	shell:
		"""
		bedtools coverage -a {params.TargeRegion} -b {input.bam} -sorted -g {params.versionsorted} | awk '{{print "{params.sample}\t" $4 "\t" $7}}' > {output.coverage}
		"""


rule targt_cov_d_bam:
	input:
		i1="{outpath}/02_map/{bam_type}/{sample}/{sample}.{bam_type}.bam"
	output:
		o1="{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.coverage.d"
	log:
		"{outpath}/02_map/logs/{sample}.targt_cov_d_{bam_type}bam.log"
	params:
		TargeRegion=config['TargeRegion'],
		versionsorted=config['versionsorted']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/bedtools.sif"
	shell:
		"""
		bedtools coverage -a {params.TargeRegion} -b {input.i1} -d -sorted -g {params.versionsorted} > {output.o1}
		"""

rule mosaic_cov_d_bam:
	input:
		i1="{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.coverage.d"
	output:
		o1="{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.coverage.mosaic.txt"
	log:
		"{outpath}/02_map/logs/{sample}.mosaic_cov_d_{bam_type}bam.log"		
	params:
		MosaicRegion=config['MosaicRegion'],
		sample="{sample}"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/bedtools.sif"
	shell:
		"""
		awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' {params.MosaicRegion} <(awk 'BEGIN{{OFS="\t"}}{{$2 = $2 + $7; $3=$2}}1' {input.i1}) | awk '{{print $8}}' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/g' | sort -k2,2n | awk '{{print "{params.sample}\t"$0}}'> {output.o1}
		"""

rule mosaic_cov_d_all_bam:
	input:
		i1=expand(["{outpath}/02_map/06_stat/02_target/{u.sample}.{bam_type}.target.coverage.mosaic.txt"],
		outpath=Outpath,
		bam_type=["02_sort","03_rmdup"],
		u=units.itertuples())
	output:
		o1="{outpath}/02_map/06_stat/02_target/all.{bam_type}.target.coverage.mosaic.txt",
		o2="{outpath}/02_map/06_stat/02_target/plot/all.{bam_type}.target.coverage.mosaic.pdf"
	log:
		"{outpath}/02_map/logs/all.mosaic_cov_d_all_{bam_type}bam.log"	 
	params:
		meta_sample=config['meta_sample'],
		depth_cumcov_mosaic="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/depth_cumcov_mosaic_v1.0.R"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/bedtools.sif"
	shell:
		"""
		awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1]"\t"$0}}}}' {params.meta_sample} <(cat {input.i1}) > {output.o1} && Rscript {params.depth_cumcov_mosaic} -i {output.o1} -d 1000 -o {output.o2}
		"""


rule targt_cov_hist_bam:
	input:
		i1="{outpath}/02_map/{bam_type}/{sample}/{sample}.{bam_type}.bam"
	output:
		o2="{outpath}/02_map/06_stat/02_target/{sample}.{bam_type}.target.coverage.hist"
	log:
		"{outpath}/02_map/logs/{sample}.targt_cov_hist.{bam_type}.bam.log"
	params:
		TargeRegion=config['TargeRegion'],
		versionsorted=config['versionsorted']
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/bedtools.sif"
	shell:
		"""
		bedtools coverage -a {params.TargeRegion} -b {input.i1} -hist -sorted -g {params.versionsorted} | grep "^all" > {output.o2}
		"""

rule exon_stat_sort:
	input:
		i0="{outpath}/02_map/06_stat/02_target/{sample}.02_sort.target.intersect.bed",
		i4="{outpath}/02_map/06_stat/02_target/{sample}.02_sort.target.v.intersect.bed"
	output:
		o1="{outpath}/02_map/06_stat/02_target/stat/sort/{sample}.stat.02_sort.txt"
	log:
		"{outpath}/02_map/logs/{sample}.exon_stat_sort.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		read_on_target_sort=$(cut -f 4 {input.i0} | sort | uniq | wc -l)
		read_off_target_sort=$(cut -f 4 {input.i4} | sort | uniq | wc -l)
		echo -e "No._read_on_target_before_dedup\t$read_on_target_sort\nNo._read_off_target_before_dedup\t$read_off_target_sort\n" > {output.o1}
		"""

rule exon_stat_dup:
	input:
		i01="{outpath}/02_map/06_stat/02_target/{sample}.03_rmdup.target.intersect.bed",
		i41="{outpath}/02_map/06_stat/02_target/{sample}.03_rmdup.target.v.intersect.bed"
	output:
		o1="{outpath}/02_map/06_stat/02_target/stat/dup/{sample}.stat.03_rmdup.txt"
	log:
		"{outpath}/02_map/logs/{sample}.exon_stat_dup.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		read_on_target_dup=$(cut -f 4 {input.i01} | sort | uniq | wc -l)
		read_off_target_dup=$(cut -f 4 {input.i41} | sort | uniq | wc -l)
		echo -e "No._read_on_target_after_dedup\t$read_on_target_dup\nNo._read_off_target_after_dedup\t$read_off_target_dup\n" > {output.o1}
		"""

rule exon_hist_stat:
	input:
		i1="{outpath}/02_map/06_stat/02_target/{sample}.02_sort.target.coverage.hist",
		i2="{outpath}/02_map/06_stat/02_target/{sample}.03_rmdup.target.coverage.hist",
		i3="{outpath}/02_map/06_stat/02_target/stat/sort/{sample}.stat.02_sort.txt",
		i4="{outpath}/02_map/06_stat/02_target/stat/dup/{sample}.stat.03_rmdup.txt"
	output:
		o3="{outpath}/02_map/06_stat/02_target/stat/{sample}.target.stat.txt",
		o1="{outpath}/02_map/06_stat/02_target/stat/sort/{sample}.hist.02_sort.txt",
		o2="{outpath}/02_map/06_stat/02_target/stat/dup/{sample}.hist.03_rmdup.txt"
	log:
		"{outpath}/02_map/logs/{sample}.exon_hist_stat.log"
	params:
		sample="{sample}",
		target_depth_stat="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/target_depth_stat_v2.R"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/r.sif"
	shell:
		"""
		Rscript {params.target_depth_stat} --shist {input.i1} --dhist {input.i2} --dup {input.i4} --sort {input.i3} -o {output.o3}
		awk '{{print "{params.sample}\t"$0}}' {input.i2} | cut -f 1,3- > {output.o2}
		awk '{{print "{params.sample}\t"$0}}' {input.i1} | cut -f 1,3- > {output.o1}
		"""

rule exon_hist_plot:
	input:
		i1=expand(["{outpath}/02_map/06_stat/02_target/stat/sort/{u.sample}.hist.02_sort.txt"], outpath=Outpath, u=units.itertuples()),
		i2=expand(["{outpath}/02_map/06_stat/02_target/stat/dup/{u.sample}.hist.03_rmdup.txt"], outpath=Outpath, u=units.itertuples())
	output:
		o1="{outpath}/02_map/06_stat/02_target/plot/allsample.02_sort.target.coverage.hist.pdf",
		o2="{outpath}/02_map/06_stat/02_target/plot/allsample.03_rmdup.target.coverage.hist.pdf",
		o3="{outpath}/02_map/06_stat/02_target/stat/allsample.02_sort.histgram.txt",
		o4="{outpath}/02_map/06_stat/02_target/stat/allsample.03_rmdup.histgram.txt"
	log:
		"{outpath}/02_map/logs/all.exon_hist_plot.log"
	params:
		meta_sample=config['meta_sample'],
		depth_cumcov_target="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/depth_cumcov_target_v1.0.R"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/r.sif"
	shell:
		"""
		awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1] "\t" $0}}}}' {params.meta_sample} <(cat {input.i1}) > {output.o3} && Rscript {params.depth_cumcov_target} -i {output.o3} -d 1000 -o {output.o1}
		awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1] "\t" $0}}}}' {params.meta_sample} <(cat {input.i2}) > {output.o4} && Rscript {params.depth_cumcov_target} -i {output.o4} -d 1000 -o {output.o2}	 
		"""

rule probe_plot:
	input:
		expand(["{outpath}/02_map/06_stat/02_target/{u.sample}.03_rmdup.target.coverage"], outpath=Outpath, u=units.itertuples())
	output:
		o1="{outpath}/02_map/06_stat/02_target/stat/all.03_rmdup.probe.coverage",
		o2="{outpath}/02_map/06_stat/02_target/plot/all.03_rmdup.probe.coverage.pdf"
	log:
		"{outpath}/02_map/logs/all.probe_plot.log"
	params:
		meta_sample=config['meta_sample'],
		probe_cov="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/probe_cov_v1.0.R"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	container:
		"../envs/r.sif"
	shell:
		"""
		awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1]"\t"$0}}}}' {params.meta_sample} <(cat {input}) > {output.o1} && Rscript {params.probe_cov} -i {output.o1} -o {output.o2}
		"""

rule map_stat_summary:
	input:
		i1=expand(
			["{outpath}/02_map/06_stat/02_target/stat/{u.sample}.target.stat.txt"],
			outpath=Outpath,
			u=units.itertuples()
		),
		i21=expand(
			["{outpath}/02_map/06_stat/02_target/plot/allsample.{bam_type}.target.coverage.hist.pdf"],
			outpath=Outpath,
			bam_type=["02_sort","03_rmdup"]
		),
		i3="{outpath}/02_map/06_stat/02_target/plot/all.03_rmdup.probe.coverage.pdf",
		i4="{outpath}/02_map/06_stat/02_target/plot/all.03_rmdup.target.coverage.mosaic.pdf",
		i5=expand(["{outpath}/02_map/06_stat/01_insert/01_insert/{u.sample}.bqsr.insert.png"],
			outpath=Outpath,
			u=units.itertuples())
	output:
		"{outpath}/02_map/07_summary/summary_target_dedup.txt"
	log:
		"{outpath}/02_map/logs/all.summary.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	shell:
		"""
		cat <(cat {input.i1} | grep "sample_id" | head -n 1) <(cat {input.i1} | grep -v "sample_id") > {output}
		"""
	