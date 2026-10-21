# config file ---------------------------------------------------------------------------------------------------
configfile: '/scratch1/tanxinli/rnaseq/RNAseq/config.yaml'

################################################################################################
# config args #######################################################################################

# configs from arguements
dir_in=config["dir_in"]
dir_out=config["dir_out"]
samplelist=config["sample"]

# configs from config file
dir_script=config['dir_script']
ref_fa=config['ref_fa']
gtf=config['gtf']
bed_DUX4=config['bed_DUX4']
ref_star=config['ref_star']
ref_fusioncatcher=config['ref_fusioncatcher']
ref_RNApeg_fa=config['ref_RNApeg_fa']
ref_RNApeg_flat=config['ref_RNApeg_flat']

cores_star=config['cores_star']
cores_samtoolsSort=config['cores_samtoolsSort']
cores_fusioncatcher=config['cores_fusioncatcher']
cores_RNApeg=config['cores_RNApeg']

# human chromosomes
set_chr=['chr{}'.format(x) for x in list(range(1,23)) + ['X', 'Y']]


################################################################################################
# Pipeline ######################################################################################
rule all:
   input:
       expand(dir_out + "/log/{sample}.Done.txt", sample=samplelist)

rule SplitNCigarReads:
    input:
        bam=    dir_in + "/{sample}.MarkDup.sorted.aligned.bam",
        bai=    dir_in + "/{sample}.MarkDup.sorted.aligned.bai",
    output:
        bam_split=  temp(dir_out + "/SplitNCigarReads/{sample}.{chr}.split.bam"),
        bai_split=  temp(dir_out + "/SplitNCigarReads/{sample}.{chr}.split.bai")
    log:        dir_out + "/log/SplitNCigarReads/{sample}.{chr}.SplitNCigarReads.log"
    benchmark:  dir_out + "/log/SplitNCigarReads/{sample}.{chr}.SplitNCigarReads.bmk"
    params:
        ref=    ref_fa,
        chr=    "{chr}",
    shell:
        '''
        gatk SplitNCigarReads -R {params.ref} -I {input.bam} --intervals {params.chr} -O {output.bam_split} &> {log}
        '''

rule HaplotypeCaller:
    input:
        bam=    dir_out + "/SplitNCigarReads/{sample}.{chr}.split.bam",
        bai=    dir_out + "/SplitNCigarReads/{sample}.{chr}.split.bai",
    output:
        vcf=    dir_out + "/HaplotypeCaller/{sample}/{sample}.{chr}.HaplotypeCaller.vcf"
    log:        dir_out + "/log/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.log"
    benchmark:  dir_out + "/log/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.bmk"
    params:
        ref=    ref_fa,
        chr=    "{chr}",
    shell:
        '''
        gatk HaplotypeCaller --native-pair-hmm-threads {threads} \
        --dont-use-soft-clipped-bases \
        -R {params.ref} --intervals {params.chr} -I {input.bam} \
        -O {output.vcf} &> {log}
        '''

rule HaplotypeCaller_merge:
    input:
        vcf=lambda wildcards: expand(
            dir_out + "/HaplotypeCaller/{sample}/{sample}.{chr}.HaplotypeCaller.vcf",
            chr=set_chr,sample=wildcards.sample)
    output:
        vcf=    dir_out + "/Mutation/{sample}.HaplotypeCaller.vcf",
    log:        dir_out + "/log/HaplotypeCaller_merge/{sample}.HaplotypeCaller_merge.log"
    benchmark:  dir_out + "/log/HaplotypeCaller_merge/{sample}.HaplotypeCaller_merge.bmk"
    params:
        filename_prefix=    dir_out + "/HaplotypeCaller/{sample}/{sample}.*.HaplotypeCaller.vcf",
        idx_prefix=         dir_out + "/HaplotypeCaller/{sample}/{sample}.*.HaplotypeCaller.vcf.idx",
        filename_list=      dir_out + "/HaplotypeCaller/{sample}/vcf.list",
        dir_HP=             dir_out + "/HaplotypeCaller/",
    shell:
        '''
        ls {params.filename_prefix} > {params.filename_list}
        gatk MergeVcfs -I {params.filename_list} -O {output.vcf} &> {log}
       
        '''
rule Done:
    input:
        vcf=            dir_out + "/Mutation/{sample}.HaplotypeCaller.vcf"
    output:
        dir_out + "/log/{sample}.Done.txt"
    params:
        log_dir=dir_out + "/log/"
    shell:
        '''
        touch {output}
        '''