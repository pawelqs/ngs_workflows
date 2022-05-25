rule facets_all:
    input: expand("results/facets/{sample}.csv", sample=tumor_samples)


rule snp_pile:
    input:
        tumor = lambda wildcards: bams[wildcards.sample],
        normal = lambda wildcards: bams[paired_normals[wildcards.sample]]
    output: "results/facets/{sample}.pileup.csv.gz"
    params:
        snps = config["FACETS_snp_vcf"]
    threads: 10
    shell: "snp-pileup -g -q10 -Q10 -P100 -r25,0 -d10000 {params.snps} {output} {input.normal} {input.tumor}"


rule facets:
    input: "results/facets/{sample}.pileup.csv.gz"
    output: "results/facets/{sample}.0csv"
    params:
        outdir = "results/facets/",
        nhet = 15,
        sample_name = "{sample}"
    conda: "../envs/facets_env.yaml"
    log: "outputs/run_facets_{sample}.txt"
    threads: 5
    script: "../scripts/run_facets.R"

    
rule facets_rep_chrom_names:
    input: "results/facets/{sample}.0csv"
    output: "results/facets/{sample}.csv"
    shell: "sed 's/^23/X/' {input} | sed 's/^24/Y/' > {output}"
