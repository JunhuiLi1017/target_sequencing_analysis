rule venn_plot:
    input:
        mut_var="result/03_Variants/Mutect2/{sample}.mutect2.pass.snv.recode.vcf",
        pis_var="result/03_Variants/pisces/{sample}.pisces.pass.snv.recode.vcf",
        rep_var="result/04_overlap/00_input/{sample}.replow.pass.snv.txt"
    output:
        o1="result/04_overlap/00_input/{sample}.mutect2.diagram.txt",
        o2="result/04_overlap/00_input/{sample}.pisces.diagram.txt",
        o3="result/04_overlap/00_input/{sample}.replow.diagram.txt",
        o4="result/04_overlap/02_plot/{sample}.venn.diagram.png"
    shell:
        r'''
        cut -f 1,2,4,5 {input.mut_var} | sed "s/\t/_/g" | grep -v "^#" > {output.o1}
        cut -f 1,2,4,5 {input.pis_var} | sed "s/\t/_/g" > {output.o2}
        cut -f 1,2,4,5 {input.rep_var} | sed "s/\t/_/g" > {output.o3}
        if [ ! -s {output.o3} ]; then
            echo "empty" > {output.o3}
        fi
        Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/variants_venn_v1.0.R --mutect2 {output.o1} --pisces {output.o2} --replow {output.o3} --outfile {output.o4}
        '''

rule scatter_plot:
    input:
        mut_var="result/03_Variants/Mutect2/{sample}.mutect2.pass.snv.recode.vcf",
        pis_var="result/03_Variants/pisces/{sample}.pisces.pass.snv.recode.vcf",
        rep_var="result/04_overlap/00_input/{sample}.replow.pass.snv.txt"
    output:
        o1="result/04_overlap/00_input/{sample}.scatter.txt",
        o2="result/04_overlap/02_plot/{sample}.scatter.png"
    params:
        sample="{sample}"
    shell:
        r'''
        cut -f 1,2,4,5,10 {input.pis_var} | sed 's/:/\t/g' | cut -f 1-4,8,9 | grep -v "^#" | awk '{{print "pisces\t"$1"_"$2"_"$3"_"$4"\t"$5"\t"$6}}' > {output.o1}
        cut -f 1,2,4,5,10 {input.mut_var} | sed 's/:/\t/g' | cut -f 1-4,7,8 | grep -v "^#" | awk '{{print "mutect2\t"$1"_"$2"_"$3"_"$4"\t"$6"\t"$5}}' >> {output.o1}
        cut -f 1,2,4,5,8,10 {input.rep_var} | grep -v "^#" | awk '{{print "replow\t"$1"_"$2"_"$3"_"$4"\t"$5"\t"$6}}' >> {output.o1}
        Rscript /home/junhui.li11-umw/Pipeline/MosaicPipeline/bin/variants_scatter_v1.0.R --input {output.o1} --outfile {output.o2}
        '''

rule venn_plot1:
    input:
        i1=expand(["result/04_overlap/02_plot/{u.sample}.venn.diagram.png"], u=units.itertuples()),
        i2=expand(["result/04_overlap/00_input/{u.sample}.pisces.diagram.txt"], u=units.itertuples()),
        i3=expand(["result/04_overlap/02_plot/{u.sample}.scatter.png"], u=units.itertuples())
    output:
        o1="result/04_overlap/01_stat/all.summary.txt"
    shell:
        '''
        cat {input.i2} > {output.o1}
        '''

