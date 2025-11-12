rule report:
	input:
		multiqc_fq="{outpath}/01_multiqc/multiqc_report.html",
		map_stat="{outpath}/02_map/07_summary/summary_target_dedup.txt",
		stats_report="{outpath}/02_map/07_summary/report_stats/multiqc_report.html",
		idxstats_report="{outpath}/02_map/07_summary/report_idxstats/multiqc_report.html",
		flagstat_report="{outpath}/02_map/07_summary/report_flagstat/multiqc_report.html",
		venn_diagram=expand(
			["{outpath}/04_overlap/01_plot/{u.sample}.venn.diagram.png"],
			outpath=Outpath,
			u=units.itertuples()
		),
		scatter_plot=expand(
			["{outpath}/04_overlap/01_plot/{u.sample}.scatter.png"],
			outpath=Outpath,
			u=units.itertuples()
		)
	output:
		o1="{outpath}/05_report/report.html"
	log:
		"{outpath}/05_report/logs/report.log"
	threads:
		resource['resource']['medium']['threads']
	resources:
		mem_mb=resource['resource']['medium']['mem_mb']
	params:
		report_ipynb = "/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/lib/report.ipynb"
	container:
		"../envs/r.sif"
	shell:
		'''
		jupyter nbconvert --to html --output {output.o1} {params.report_ipynb} > {log} 2>&1
		''' 