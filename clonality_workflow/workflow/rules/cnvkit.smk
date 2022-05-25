
def get_cnvkit_info(patient):
    patient_samples = samples[samples.patient_id == patient]
    tumor_samples = patient_samples[patient_samples.sample_type == "tumor"].sample_id.tolist()
    normal_samples = patient_samples[patient_samples.sample_type == "normal"].sample_id.tolist()
    tumor_bams = patient_samples[patient_samples.sample_type == "tumor"].bam.tolist()
    normal_bams = patient_samples[patient_samples.sample_type == "normal"].bam.tolist()
    sex = patient_samples["sex"].tolist()[1]
    res = {
        "tumor_samples": tumor_samples,
        "normal_samples": normal_samples,
        "tumor_bams": tumor_bams,
        "normal_bams": normal_bams,
        "sex": sex
    }
    return(res)

patient_cnvkit = {patient: get_cnvkit_info(patient) for patient in patients}
access_file = "cnvkit/excludes." + config["genomeBuild"] + ".bed"


rule cnvkit_all:
    input: 
        expand("cnvkit/{patient}/done", patient=patients)


rule cnvkit_access:
    input: config["genome_fa"]
    output: access_file
    conda: "../envs/cnvkit_env.yaml"
    shell:
        "cnvkit.py access {input} -o {output}"


rule cnvkit_pipeline:
    input:
        normal = lambda wildcards: patient_cnvkit[wildcards.patient]["normal_bams"],
        tumor = lambda wildcards: patient_cnvkit[wildcards.patient]["tumor_bams"],
        genome_excludes = access_file
    output:
        dir = directory("cnvkit/{patient}/"),
        done = "cnvkit/{patient}/done"
    params:
        targets = config["targets_bed"],
        genome_fa = config["genome_fa"],
        cnn_ref = "cnvkit/{patient}/reference.cnn"
    conda: "../envs/cnvkit_env.yaml"
    threads: 10
    shell:
        """
        mkdir -p {output.dir}
        srun cnvkit.py batch {input.tumor} \\
            --normal {input.normal} \\
            --targets {params.targets} \\
            --fasta {params.genome_fa} \\
            --access {input.genome_excludes} \\
            --output-reference {params.cnn_ref} \\
            --output-dir {output.dir} \\
            --diagram --scatter \\
            -p {threads}
        find {output.dir} -name '*.cns' ! -name '*bintest*' ! -name '*call*' -exec cp {{}} cnvkit/ \;
        touch {output.done}
        """
