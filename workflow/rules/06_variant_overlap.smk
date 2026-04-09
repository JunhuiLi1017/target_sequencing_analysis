rule venn_input:
	input:
		get_vcf_inputs
	output:
		txt="{outpath}/04_overlap/01_plot/{sample}.venn.diagram.txt"
	log:
		"{outpath}/04_overlap/logs/{sample}.venn.txt.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	params:
		venn_input_cmd=get_venn_input
	shell:
		r"""
		{params.venn_input_cmd}
		"""

rule venn_input_filter_indel:
	input:
		txt="{outpath}/04_overlap/01_plot/{sample}.venn.diagram.txt"
	output:
		txt="{outpath}/04_overlap/01_plot/{sample}.venn.diagram.SNV.txt"
	log:
		"{outpath}/04_overlap/logs/{sample}.venn.snv.txt.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	shell:
		"""
		awk '{{split($2, a, "_");if (length(a[3]) == length(a[4]) ) {{ print}}}}' {input.txt} > {output.txt}
		"""

rule venn_plot:
	input:
		"{outpath}/04_overlap/01_plot/{sample}.venn.diagram.SNV.txt"
	output:
		venn="{outpath}/04_overlap/01_plot/{sample}.venn.diagram.png"
	log:
		"{outpath}/04_overlap/logs/{sample}.venn.diagram.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	params:
		venn_rscript="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/variants_venn_v1.0.R"
	shell:
		r"""
		Rscript {params.venn_rscript} --input {input} --outfile {output.venn} > {log} 2>&1
		"""

rule scatter_input:
	input:
		get_vcf_inputs
	output:
		txt="{outpath}/04_overlap/01_plot/{sample}.scatter.txt"
	log:
		"{outpath}/04_overlap/logs/{sample}.scatter.input.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	params:
		scatter_cmd=get_scatter_cmd,
		output_dir="{outpath}/04_overlap/01_plot"
	shell:
		r"""
		{params.scatter_cmd}
		"""

rule scatter_plot:
	input:
		txt="{outpath}/04_overlap/01_plot/{sample}.scatter.txt"
	output:
		plot="{outpath}/04_overlap/01_plot/{sample}.scatter.png"
	log:
		"{outpath}/04_overlap/logs/{sample}.scatter.plot.log"
	threads:
		resource['resource']['low']['threads']
	resources:
		mem_mb=resource['resource']['low']['mem_mb']
	params:
		variants_scatter_script="/pi/michael.lodato-umw/junhui.li11-umw/BautistaSotelo_Cesar/20201130_MosaicVariant_DNA/00script/00_pipeline/target_sequence_analysis/workflow/bin/variants_scatter_v1.0.R"
	shell:
		r"""
		Rscript {params.variants_scatter_script} --input {input.txt} --outfile {output.plot} > {log} 2>&1
		"""