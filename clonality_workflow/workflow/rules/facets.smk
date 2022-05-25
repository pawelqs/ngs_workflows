rule facets_all:
    input: expand(outdir + "/facets/{sample}.csv", sample=tumor_samples)


rule snp_pile:
    input:
        tumor = lambda wildcards: bams[wildcards.sample],
        normal = lambda wildcards: bams[paired_normals[wildcards.sample]]
    output: outdir + "/facets/{sample}.pileup.csv.gz"
    params:
        snps = config["FACETS_snp_vcf"]
    threads: 10
    shell: "snp-pileup -g -q10 -Q10 -P100 -r25,0 -d10000 {params.snps} {output} {input.normal} {input.tumor}"


rule facets:
    input: outdir + "/facets/{sample}.pileup.csv.gz"
    output: outdir + "/facets/{sample}.0csv"
    params:
        outdir = outdir + "/facets/",
        nhet = 15,
        sample_name = "{sample}"
    conda: "../envs/facets_env.yaml"
    log: "outputs/run_facets_{sample}.txt"
    threads: 5
    script: "../scripts/run_facets.R"

    
rule facets_rep_chrom_names:
    input: outdir + "/facets/{sample}.0csv"
    output: outdir + "/facets/{sample}.csv"
    shell: "sed 's/^23/X/' {input} | sed 's/^24/Y/' > {output}"
