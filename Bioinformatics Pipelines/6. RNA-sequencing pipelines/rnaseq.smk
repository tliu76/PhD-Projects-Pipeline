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

rule Star:
    input:
        fq1=        dir_in + '/{sample}/{sample}_1.fq.gz',
        fq2=        dir_in + '/{sample}/{sample}_2.fq.gz',
    output:
        bam=        temp(dir_out + "/bam/{sample}.unsort.bam"),
        junction=   dir_out + "/bam/{sample}.star.junction",
    benchmark:      dir_out + "/log/{sample}.star.bmk",
    log:            dir_out + "/bam/{sample}.star.log",
    params:
        star_ref=       ref_star,
        rg=             "ID:{sample} PL:ILLUMINA PU:ILLUMINA LB:RNAseq SM:{sample}",
        dir_out=        dir_out + "/star/",
        bam_tmp=        dir_out + "/star/{sample}.aligned.out.bam",
        junction_tmp=   dir_out + "/star/{sample}.chimeric.out.junction"
    threads: cores_star
    shell:
        '''
        STAR \
        --twopassMode Basic \
        --runThreadN {threads} \
        --readFilesCommand zcat \
        --limitOutSAMoneReadBytes 90000000 \
        --outSAMtype BAM Unsorted \
        --outSAMstrandField intronMotif \
        --outSAMunmapped Within \
        --outSAMmapqUnique 60 \
        --outSAMmultNmax 1 \
        --outReadsUnmapped None \
        --outMultimapperOrder Random \
        --outSAMattributes NH HI AS nM NM MD \
        --outSAMattrRGline {params.rg} \
        --outFilterType BySJout \
        --outFilterMultimapNmax 100 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 4 \
        --alignInsertionFlush Right \
        --chimOutType Junctions \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 100 \
        --chimMultimapScoreRange 10 \
        --chimNonchimScoreDropMin 10 \
        --chimOutJunctionFormat 1 \
        --peOverlapNbasesMin 12 \
        --peOverlapMMp 0.1 \
        --outFileNamePrefix {params.dir_out} \
        --genomeDir {params.star_ref} \
        --genomeLoad NoSharedMemory \
        --readFilesIn {input.fq1} {input.fq2} \
        --sjdbOverhang 100 &> {log}

        mv {params.bam_tmp} {output.bam}
        mv {params.junction_tmp} {output.junction}
        rm -rf {params.dir_out}
        '''


rule samtools_sort:
    input:
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.unsort.bam",
    output:
        bam=        temp(dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.sorted.bam"),
        bai=        temp(dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.sorted.bam.bai"),
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/samtools_sort.bmk",
    params:
        dir=        dir_out + "/{sample}/TRANSCRIPTOME/bam/"
    threads: cores_samtoolsSort
    shell:
        '''
        samtools sort -m 1G -@ {threads} -O bam -T {params.dir} -o {output.bam} {input.bam}
        samtools index {output.bam} {output.bai} 
        '''

rule MarkDup:
    input:
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.sorted.bam",
        bai=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.sorted.bam.bai",
    output:
        bam=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai=        dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
        metrics=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.MarkDup.metrics.txt",
    log:            dir_out + "/{sample}/TRANSCRIPTOME/log/markDup.log"
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/markDup.bmk"
    params:
        bam_star=   dir_out + "/{sample}/TRANSCRIPTOME/star/{sample}.bam",
    shell:
        '''        
        gatk MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} &> {log}
        samtools index {output.bam} {output.bai}
        '''

rule HTSeq_count:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
    output:
        HTSeq_count=    dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.HTSeq",
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/HTSeq_count.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/HTSeq_count.bmk"
    params:
        gtf=    gtf
    shell:
        '''
        htseq-count -m intersection-strict -f bam -r pos -s no -a 10 {input.bam} {params.gtf} 1> {output.HTSeq_count} 2> {log}
        '''

rule HTSeqDUX4:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
        htseq=  dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.HTSeq",
    output:
        HTSeq=  dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.DUX4patched.HTSeq",
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/HTSeqDUX4.bmk"
    params:
        bed=    bed_DUX4,
        script= dir_script + "HTSeqDUX4patch.pl"
    shell:
        '''
        {params.script} -f {params.bed} -b {input.bam} -h {input.htseq} -o {output.HTSeq}
        '''

rule FusionCatcher:
    input:
        fq= lambda wildcards: expand(dir_in + '/{sample}.R{i}.fq.gz',i=[1,2],sample=wildcards.sample)
    output:
        dir=    directory(dir_out + "/{sample}/TRANSCRIPTOME/FusionCatcher")
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/FusionCatcher.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/FusionCatcher.bmk"
    params:
        ref=        ref_fusioncatcher,
        dir_out=    dir_out + "/{sample}/TRANSCRIPTOME/FusionCatcher",
        fq1=        dir_in + '/{sample}.R1.fq.gz',
        fq2=        dir_in + '/{sample}.R2.fq.gz',
    threads: cores_fusioncatcher
    shell:
        '''        
        fusioncatcher.py -p {threads} -d {params.ref} -i {params.fq1},{params.fq2} -o {output.dir} &>{log}
        '''

rule RNApeg:
    input:
        bam = dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai = dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
    output:
        junctions=  dir_out + "/{sample}/TRANSCRIPTOME/RNApeg/{sample}.bam.junctions.tab.shifted.tab"
    log:            dir_out + "/{sample}/TRANSCRIPTOME/log/RNApeg.log"
    benchmark:      dir_out + "/{sample}/TRANSCRIPTOME/log/RNApeg.bmk"
    params:
        ref1=    ref_RNApeg_fa,
        ref2=    ref_RNApeg_flat,
        dir_out= dir_out + "/{sample}/TRANSCRIPTOME/RNApeg",
        dir_bam= dir_out + "/{sample}/TRANSCRIPTOME/bam",
    threads: cores_RNApeg
    shell:
        '''
        RNApeg.sh -b {input.bam} -f {params.ref1} -r {params.ref2} -O {params.dir_out}  &> {log}
        '''

rule SplitNCigarReads:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam",
        bai=    dir_out + "/{sample}/TRANSCRIPTOME/bam/{sample}.bam.bai",
    output:
        bam_split=  temp(dir_out + "/{sample}/TRANSCRIPTOME/SplitNCigarReads/{sample}.{chr}.split.bam"),
        bai_split=  temp(dir_out + "/{sample}/TRANSCRIPTOME/SplitNCigarReads/{sample}.{chr}.split.bai")
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/SplitNCigarReads/{chr}.SplitNCigarReads.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/SplitNCigarReads/{chr}.SplitNCigarReads.bmk"
    params:
        ref=    ref_fa,
        chr=    "{chr}",
    shell:
        '''
        gatk SplitNCigarReads -R {params.ref} -I {input.bam} --intervals {params.chr} -O {output.bam_split} &> {log}
        '''

rule HaplotypeCaller:
    input:
        bam=    dir_out + "/{sample}/TRANSCRIPTOME/SplitNCigarReads/{sample}.{chr}.split.bam",
        bai=    dir_out + "/{sample}/TRANSCRIPTOME/SplitNCigarReads/{sample}.{chr}.split.bai",
    output:
        vcf=    dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.vcf"
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller/{chr}.HaplotypeCaller.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller/{chr}.HaplotypeCaller.bmk"
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
            dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.{chr}.HaplotypeCaller.vcf",
            chr=set_chr,sample=wildcards.sample)
    output:
        vcf=    dir_out + "/{sample}/TRANSCRIPTOME/Mutation/{sample}.HaplotypeCaller.vcf",
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller_merge.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/HaplotypeCaller_merge.bmk"
    params:
        filename_prefix=    dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.*.HaplotypeCaller.vcf",
        idx_prefix=         dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/{sample}.*.HaplotypeCaller.vcf.idx",
        filename_list=      dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller/vcf.list",
        dir_HP=             dir_out + "/{sample}/TRANSCRIPTOME/HaplotypeCaller",
    shell:
        '''
        ls {params.filename_prefix} > {params.filename_list}
        gatk MergeVcfs -I {params.filename_list} -O {output.vcf} &> {log}
        rm -rf {params.dir_HP}
        '''

rule RNAseqCNV:
    input:
        vcf=        dir_out + "/{sample}/TRANSCRIPTOME/Mutation/{sample}.HaplotypeCaller.vcf",
        counts=     dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.DUX4patched.HTSeq",
    output:
        cnv=    directory(dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV/{sample}"),
    log:        dir_out + "/{sample}/TRANSCRIPTOME/log/RNAseqCNV.log"
    benchmark:  dir_out + "/{sample}/TRANSCRIPTOME/log/RNAseqCNV.bmk"
    params:
        config_file=    dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV/config",
        outdir=         dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV",
        count_dir=      dir_out + "/{sample}/TRANSCRIPTOME/HTSeq",
        snv_dir=        dir_out + "/{sample}/TRANSCRIPTOME/Mutation",
        meta_file=      dir_out + "/{sample}/TRANSCRIPTOME/RNAseqCNV/metadata",
        sample=         "{sample}",
        htseq_file=     "{sample}.DUX4patched.HTSeq",
        vcf_file=       "{sample}.HaplotypeCaller.vcf",
        script=         dir_script + "RNAseqCNV.R"
    shell:
        '''
        cat <(echo 'out_dir = "{params.outdir}"') <(echo 'count_dir = "{params.count_dir}"') <(echo 'snv_dir = "{params.snv_dir}"')  > {params.config_file}
        cat <(echo '{params.sample} {params.htseq_file} {params.vcf_file}') > {params.meta_file}
        Rscript_GEP {params.script} {params.config_file} {params.meta_file} &> {log}
        '''

rule Done:
    input:
        HTSeq_count=    dir_out + "/{sample}/TRANSCRIPTOME/HTSeq/{sample}.DUX4patched.HTSeq",
        vcf=            dir_out + "/{sample}/TRANSCRIPTOME/Mutation/{sample}.HaplotypeCaller.vcf",
        RNApeg=         dir_out + "/{sample}/TRANSCRIPTOME/log/RNApeg.bmk",
        fusioncatcher=  dir_out + "/{sample}/TRANSCRIPTOME/log/FusionCatcher.bmk",
        RNAseqCNV=      dir_out + "/{sample}/TRANSCRIPTOME/log/RNAseqCNV.bmk",
    output:
        dir_out + "/{sample}/TRANSCRIPTOME/log/Done.txt"
    params:
        log_dir=dir_out + "/{sample}/TRANSCRIPTOME/log/"
    shell:
        '''
        touch {output}
        '''

