# Somatic SNV and indel calling (Mutect2)

Somatic short-variant discovery on matched tumor-normal WGS BAMs produced by [GATK data pre-processing](../gatk_data_preprocessing/README.md). The workflow follows the GATK [somatic short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels) best practices, scattered across intervals so that one WGS pair finishes in hours rather than days.

All scripts are SLURM array jobs written against GATK 4.2.6.1.

