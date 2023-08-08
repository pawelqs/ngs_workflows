# It does not work with the ../envs/ascat_env.yml environment: ASCAT needs to be updated with
# devtools::install_github('VanLoo-lab/ascat/ASCAT') in the conda env

rule ascat_all:
    input: expand("ascat/{sample_id}.done", sample_id=tumor_samples)


sex_to_sex_chromosomes = {
    "female": "XX",
    "male": "XY"
}

ascat_genders = {
    key: sex_to_sex_chromosomes[value]
    for key, value in sex.items()
}


rule run_ascat:
    input: 
        tumor_bam = lambda wildcards: bams[wildcards.sample_id],
        normal_bam = lambda wildcards: bams[paired_normals[wildcards.sample_id]]
    output: "ascat/{sample_id}.done"
    params:
        tumor_name = "{sample_id}",
        normal_name = lambda wildcards: paired_normals[wildcards.sample_id],
        BED_file = config["targets_bed"] if "targets_bed" in config else None,
        allelecounter_path = config["allelecounter_path"],
        allele_prefix = config["ascat_allele_prefix"],
        loci_prefix = config["ascat_loci_prefix"],
        gc_file = config["ascat_gc_file"],
        rt_file = config["ascat_rt_file"],
        gender = lambda wildcards: ascat_genders[wildcards.sample_id],
        genome_version = config["genomeBuild"]
#    conda: "../envs/ascat_env.yml"
    conda: "ascat_env"
    threads: 5
    script: "../scripts/run_ASCAT.R"
