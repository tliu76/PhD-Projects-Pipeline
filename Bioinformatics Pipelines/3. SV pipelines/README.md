# Somatic structural variants and RAG-mediated deletions

This is the HPC pipeline for calling structural variations using three softwares from 12 matched tumor-normal sequencing bams from 6 childhood acute lymphoblastic leukemia cases. 

## SV calling

### Smoove (lumpy)
- [1.1 PreMergeSVPerSampleSmoove](SV_1.1_pre_merge_sv_per_sample_smoove.sh)
- [1.2 SmooveUnionSites](SV_1.2_smoove_unionsites)
- [1.3 SmooveGenotype](SV_1.3_smoove_genotype)

## Manta
- [2.1 PreMergeSVPerSampleMantaConfigure](SV_2.1_Pre_Merge_SV_Per_Sample_Manta_configure)
- [2.2 MantaWorkflow](SV_2.2_Manta_workflow)
- [2.3 MantaClean](SV_2.3_Manta_clean)

## Delly
- [3.1 DellyPerSampleCall](SV_3.1_Delly_per_sample_call)
- [3.2 DellyPrefilter](SV_3.2_Delly_prefilter)
- [3.3 DellyGeno](SV_3.3_Delly_geno)
- [3.4 DellyPostfilter](SV_3.4_Delly_postfilter)
- [3.5 DellyBcf2vcf](SV_3.5_Delly_bcf2vcf)

## Finding variants that were called by at least 2 callers
- [4.1 SURVIVOR](SV_4.1_SURVIVOR_new)

## Annotate variants
- [4.2 Annotate](SV_4.2_annotate)


## Motif enrichment at breakpoints
- (To be continued...)



