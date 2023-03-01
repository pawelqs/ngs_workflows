
rule ascat_sc_all:
    input: expand("ASCAT_sc/{patient}.Rda", patient=patients)

def get_ascat_sc_info(patient):
    patient_samples = samples[samples.patient_id == patient]
    tumor_samples = patient_samples[patient_samples.sample_type == "tumor"].sample_id.tolist()
    normal_samples = patient_samples[patient_samples.sample_type == "normal"].sample_id.tolist()
    tumor_bams = patient_samples[patient_samples.sample_type == "tumor"].bam.tolist()
    normal_bams = patient_samples[patient_samples.sample_type == "normal"].bam.tolist()
    all_bams = patient_samples.bam.tolist()
    sex = patient_samples["sex"].tolist()[1]
    res = {
        "tumor_samples": tumor_samples,
        "normal_samples": normal_samples,
        "tumor_bams": tumor_bams,
        "normal_bams": normal_bams,
        "all_bams": all_bams,
        "sex": sex
    }
    return(res)

patient_ascat_sc = {patient: get_ascat_sc_info(patient) for patient in patients}


rule run_ascat_sc:
    input:
        bams = lambda wildcards: patient_ascat_sc[wildcards.patient]["all_bams"],
    output:
        rda = "ASCAT_sc/{patient}.Rda"
    params:
        genome_fa = config["genome_fa"],
        workdir = "ASCAT_sc/",
        genome = config["genomeBuild"], # either hg19 or hg38 so far
        projectname = "{patient}",
        multipcf = True,
        sex = lambda wildcards: patient_ascat_sc[wildcards.patient]["sex"],
        chrstring_bam = config["chrstring"]
    conda: "ascat_sc"
    log: "outputs/run_ASCAT_sc_{patient}.txt"
    threads: 8
    script: "../scripts/run_ASCAT_sc.R"
