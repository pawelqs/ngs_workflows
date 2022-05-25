#!/bin/bash

# This script is not used anywhere, just may be used to check how many mutations wil be kept for certain bcftools view options
# To run it:
# ls ~/amlro_run/results/mutect2/* | xargs -i workflow/scripts/sim_vcf_filtering.sh {}

tumor_samples=`bcftools query -l $1 | grep -v ".*_Tc.*" | tr "\n" "," | sed "s/,$//"`
echo $tumor_samples


bcftools view -s $tumor_samples -Ou $1 |
            bcftools view \
                -O v \
                -i 'INFO/CSQ !~"intron_variant" & INFO/CSQ !~"synonymous_variant" &
                    INFO/CSQ !~"intergenic_variant" & INFO/CSQ !~"non_coding" & 
                    MIN(INFO/MMQ) > 50 & AVG(FORMAT/DP) > 10' | 
        grep -v '^#' | awk '$10 !~ /\.\/\./ && $9 ~ /:AD:/' |
        wc -l

# AMLRO_10_Dx,AMLRO_10_Rx
# 1014
# AMLRO_11_Dx,AMLRO_11_Rx
# 699
# AMLRO_12_Dx,AMLRO_12_Rx
# 3900
# AMLRO_15_Dx,AMLRO_15_Rx
# 805
# AMLRO_1_Dx,AMLRO_1_Rx
# 3465
# AMLRO_2_Dx,AMLRO_2_Rx
# 1117
# AMLRO_3_Dx,AMLRO_3_Rx
# 935
# AMLRO_4_Dx,AMLRO_4_Rx
# 4278
# AMLRO_6_Dx,AMLRO_6_Rx
# 577
# AMLRO_8_Dx,AMLRO_8_Rx
# 1018
# AMLRO_9_Dx,AMLRO_9_Rx
# 614