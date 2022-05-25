
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


rule cnvkit_all:
    input: 
        expand("results/cnvkit/{patient}/done", patient=patients)


rule cnvkit_pipeline:
    input:
        normal = lambda wildcards: patient_cnvkit[wildcards.patient]["normal_bams"],
        tumor = lambda wildcards: patient_cnvkit[wildcards.patient]["tumor_bams"]
    output:
        dir = directory("results/cnvkit/{patient}/"),
        done = "results/cnvkit/{patient}/done"
    params:
        targets = "resources/agilent_sureselect_v6/Agilent_SureSelect_v6r2_S07604514_Covered_hg38.cnvkit.bed",
        cnn_ref = "results/cnvkit/{patient}/reference.cnn"
    conda: "../envs/cnvkit_env.yaml"
    threads: 10
    shell:
        """
        mkdir -p {output.dir}
        srun cnvkit.py batch {input.tumor} \\
            --normal {input.normal} \\
            --targets {params.targets} \\
            --fasta ~/resources/hg38.fa \\
            --access resources/access-excludes.hg38.bed \\
            --output-reference {params.cnn_ref} \\
            --output-dir {output.dir} \\
            --diagram --scatter \\
            -p {threads}
        find {output.dir} -name '*.cns' ! -name '*bintest*' ! -name '*call*' -exec cp {{}} results/cnvkit/ \;
        touch {output.done}
        """
