rule targt_intersect_sortbam:
    input:
        i1="result/02_Map/bwa/{sample}.sort.bam"
    output:
        o1="result/02_Map/target/{sample}.sort.target.intersect.bed"
    params:
        TargeRegion=config['TargeRegion']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools intersect -abam {input.i1} -b {params.TargeRegion} -wa -bed > {output.o1}
        """

rule targt_intersectV_sortbam:
    input:
        i1="result/02_Map/bwa/{sample}.sort.bam"
    output:
        o1="result/02_Map/target/{sample}.sort.target.v.intersect.bed"
    params:
        TargeRegion=config['TargeRegion']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools intersect -abam {input.i1} -b {params.TargeRegion} -v -bed > {output.o1}
        """

rule targt_intersect_dupbam:
    input:
        i1="result/02_Map/dup/{sample}.sort.rmdup.bam"
    output:
        o1="result/02_Map/target/{sample}.dup.target.intersect.bed"
    params:
        TargeRegion=config['TargeRegion']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools intersect -abam {input.i1} -b {params.TargeRegion} -wa -bed > {output.o1}
        """

rule targt_intersectV_dupbam:
    input:
        i1="result/02_Map/dup/{sample}.sort.rmdup.bam"
    output:
        o1="result/02_Map/target/{sample}.dup.target.v.intersect.bed"
    params:
        TargeRegion=config['TargeRegion']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools intersect -abam {input.i1} -b {params.TargeRegion} -v -bed > {output.o1}
        """

rule targt_cov_sortbam:
    input:
        i1="result/02_Map/bwa/{sample}.sort.bam"
    output:
        o1="result/02_Map/target/{sample}.sort.target.coverage"
    params:
        TargeRegion=config['TargeRegion'],
        sample="{sample}",
        versionsorted=config['versionsorted']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools coverage -a {params.TargeRegion} -b {input.i1} -sorted -g {params.versionsorted} | awk '{{print "{params.sample}\t" $4 "\t" $7}}' > {output.o1}
        """

rule targt_cov_dupbam:
    input:
        i1="result/02_Map/dup/{sample}.sort.rmdup.bam"
    output:
        o1="result/02_Map/target/{sample}.dup.target.coverage"
    params:
        TargeRegion=config['TargeRegion'],
        sample="{sample}",
        versionsorted=config['versionsorted']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools coverage -a {params.TargeRegion} -b {input.i1} -sorted -g {params.versionsorted} | awk '{{print "{params.sample}\t" $4 "\t" $7}}' > {output.o1}
        """

rule targt_cov_d_sortbam:
    input:
        i1="result/02_Map/bwa/{sample}.sort.bam"
    output:
        o1="result/02_Map/target/{sample}.sort.target.coverage.d"
    params:
        TargeRegion=config['TargeRegion'],
        versionsorted=config['versionsorted']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools coverage -a {params.TargeRegion} -b {input.i1} -d -sorted -g {params.versionsorted} > {output.o1}
        """

rule targt_cov_d_dupbam:
    input:
        i1="result/02_Map/dup/{sample}.sort.rmdup.bam"
    output:
        o1="result/02_Map/target/{sample}.dup.target.coverage.d"
    params:
        TargeRegion=config['TargeRegion'],
        versionsorted=config['versionsorted']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools coverage -a {params.TargeRegion} -b {input.i1} -d -sorted -g {params.versionsorted} > {output.o1}
        """

rule mosaic_cov_d_dupbam:
    input:
        i1="result/02_Map/target/{sample}.dup.target.coverage.d"
    output:
        o1="result/02_Map/target/{sample}.dup.target.coverage.mosaic.txt"
    params:
        MosaicRegion=config['MosaicRegion'],
        sample="{sample}"
    shell:
        """
        awk 'NR==FNR{{c[$1,$2]=$0}}NR!=FNR{{if(c[$1,$2]){{print $0}}}}' {params.MosaicRegion} <(awk 'BEGIN{{OFS="\t"}}{{$2 = $2 + $7; $3=$2}}1' {input.i1}) | awk '{{print $8}}' | sort | uniq -c | sed 's/^[ \t]*//' | sed 's/ /\t/g' | sort -k2,2n | awk '{{print "{params.sample}\t"$0}}'> {output.o1}
        """

rule mosaic_cov_d_all_dupbam:
    input:
        i1=expand(["result/02_Map/target/{u.sample}.dup.target.coverage.mosaic.txt"], u=units.itertuples())
    output:
        o1="result/02_Map/target/all.dup.target.coverage.mosaic.txt",
        o2="result/02_Map/target/plot/all.dup.target.coverage.mosaic.pdf"
    params:
        meta_sample=config['meta_sample']
    shell:
        """
        awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1]"\t"$0}}}}' {params.meta_sample} <(cat {input.i1}) > {output.o1} && Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/depth_cumcov_mosaic_v1.0.R -i {output.o1} -d 1000 -o {output.o2}
        """


rule targt_cov_hist_sortbam:
    input:
        i1="result/02_Map/bwa/{sample}.sort.bam"
    output:
        o2="result/02_Map/target/{sample}.sort.target.coverage.hist"
    params:
        TargeRegion=config['TargeRegion'],
        versionsorted=config['versionsorted']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools coverage -a {params.TargeRegion} -b {input.i1} -hist -sorted -g {params.versionsorted} | grep "^all" > {output.o2}
        """

rule targt_cov_hist_dupbam:
    input:
        i1="result/02_Map/dup/{sample}.sort.rmdup.bam"
    output:
        o2="result/02_Map/target/{sample}.dup.target.coverage.hist"
    params:
        TargeRegion=config['TargeRegion'],
        versionsorted=config['versionsorted']
    log:
        "logs/bwa/{sample}.sort.target.log"
    shell:
        """
        bedtools coverage -a {params.TargeRegion} -b {input.i1} -hist -sorted -g {params.versionsorted} | grep "^all" > {output.o2}
        """

rule exon_stat_sort:
    input:
        i0="result/02_Map/target/{sample}.sort.target.intersect.bed",
        i4="result/02_Map/target/{sample}.sort.target.v.intersect.bed"
    output:
        o1="result/02_Map/target/stat/sort/{sample}.stat.sort.txt"
    log:
        "logs/bwa/{sample}.sort.exon.stat.log"
    shell:
        """
        read_on_target_sort=$(cut -f 4 {input.i0} | sort | uniq | wc -l)
        read_off_target_sort=$(cut -f 4 {input.i4} | sort | uniq | wc -l)
        echo -e "No._read_on_target_before_dedup\t$read_on_target_sort\nNo._read_off_target_before_dedup\t$read_off_target_sort\n" > {output.o1}
        """

rule exon_stat_dup:
    input:
        i01="result/02_Map/target/{sample}.dup.target.intersect.bed",
        i41="result/02_Map/target/{sample}.dup.target.v.intersect.bed"
    output:
        o1="result/02_Map/target/stat/dup/{sample}.stat.dup.txt"
    log:
        "logs/bwa/{sample}.sort.exon.stat.log"
    shell:
        """
        read_on_target_dup=$(cut -f 4 {input.i01} | sort | uniq | wc -l)
        read_off_target_dup=$(cut -f 4 {input.i41} | sort | uniq | wc -l)
        echo -e "No._read_on_target_after_dedup\t$read_on_target_dup\nNo._read_off_target_after_dedup\t$read_off_target_dup\n" > {output.o1}
        """

rule exon_hist_stat:
    input:
        i1="result/02_Map/target/{sample}.sort.target.coverage.hist",
        i2="result/02_Map/target/{sample}.dup.target.coverage.hist",
        i3="result/02_Map/target/stat/sort/{sample}.stat.sort.txt",
        i4="result/02_Map/target/stat/dup/{sample}.stat.dup.txt"
    output:
        o3="result/02_Map/target/stat/{sample}.target.stat.txt",
        o1="result/02_Map/target/stat/sort/{sample}.hist.sort.txt",
        o2="result/02_Map/target/stat/dup/{sample}.hist.dup.txt"
    params:
        sample="{sample}"
    shell:
        """
        Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/target_depth_stat_v2.R --shist {input.i1} --dhist {input.i2} --dup {input.i4} --sort {input.i3} -o {output.o3}
        awk '{{print "{params.sample}\t"$0}}' {input.i2} | cut -f 1,3- > {output.o2}
        awk '{{print "{params.sample}\t"$0}}' {input.i1} | cut -f 1,3- > {output.o1}
        """

rule exon_hist_plot:
    input:
        i1=expand(["result/02_Map/target/stat/sort/{u.sample}.hist.sort.txt"], u=units.itertuples()),
        i2=expand(["result/02_Map/target/stat/dup/{u.sample}.hist.dup.txt"], u=units.itertuples())
    output:
        o1="result/02_Map/target/plot/allsample.sort.target.coverage.hist.pdf",
        o2="result/02_Map/target/plot/allsample.dup.target.coverage.hist.pdf",
        o3="result/02_Map/target/stat/allsample.sort.histgram.txt",
        o4="result/02_Map/target/stat/allsample.dup.histgram.txt"
    params:
        meta_sample=config['meta_sample']
    shell:
        """
        awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1] "\t" $0}}}}' {params.meta_sample} <(cat {input.i1}) > {output.o3} && Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/depth_cumcov_target_v1.0.R -i {output.o3} -d 1000 -o {output.o1}
        awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1] "\t" $0}}}}' {params.meta_sample} <(cat {input.i2}) > {output.o4} && Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/depth_cumcov_target_v1.0.R -i {output.o4} -d 1000 -o {output.o2}     
        """

rule probe_plot:
    input:
        expand(["result/02_Map/target/{u.sample}.dup.target.coverage"], u=units.itertuples())
    output:
        o1="result/02_Map/target/stat/all.dup.probe.coverage",
        o2="result/02_Map/target/plot/all.dup.probe.coverage.pdf",
    params:
        meta_sample=config['meta_sample']
    shell:
        """
        awk 'NR==FNR{{c[$1]=$0}}NR!=FNR{{if(c[$1]){{print c[$1]"\t"$0}}}}' {params.meta_sample} <(cat {input}) > {output.o1} && Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/probe_cov_v1.0.R -i {output.o1} -o {output.o2}
        """

rule summary:
    input:
        i1=expand(
            ["result/02_Map/target/stat/{u.sample}.target.stat.txt"],
             u=units.itertuples()
        ),
        i21="result/02_Map/target/plot/allsample.sort.target.coverage.hist.pdf",
        i22="result/02_Map/target/plot/allsample.dup.target.coverage.hist.pdf",
        i3="result/02_Map/target/plot/all.dup.probe.coverage.pdf",
        i4="result/02_Map/target/plot/all.dup.target.coverage.mosaic.pdf"
    output:
        o1="result/02_Map/summary/summary_target_dedup.txt"
    shell:
        """
        cat <(cat {input.i1} | grep "sample_id" | head -n 1) <(cat {input.i1} | grep -v "sample_id") > {output.o1}
        """
    