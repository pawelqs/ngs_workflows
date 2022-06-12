import subprocess
import numpy as np


def get_pyclone_info(patient):
    patient_samples = samples[samples.patient_id == patient]
    tumor_samples = patient_samples[patient_samples.sample_type == "tumor"].sample_id.tolist()
    normal_samples = patient_samples[patient_samples.sample_type == "normal"].sample_id.tolist()
    sex = patient_samples["sex"].tolist()[1]
    res = {
        "tumor_samples": tumor_samples,
        "normal_samples": normal_samples,
        "sex": sex
    }
    return(res)

patients_pyclone = {patient: get_pyclone_info(patient) for patient in patients}


rule pyclone_all:
    input:
        expand("pyclone_facets/{patient}.pyvi-output.tsv", patient=patients)


def pyclone_get_cnv_files(method, patient):
    sample_ids = patients_pyclone[patient]["tumor_samples"]
    # if method == "titan":
    #     return "titan_optimal/%s/%s.segs.tsv" % (sample_id, sample_id)
    if method == "facets":
        files = ["facets/%s.csv" % sample_id for sample_id in sample_ids]
    return files


def get_patient_sex(patient):
    sex = patients_pyclone[patient]["sex"]
    if pd.isnull(sex):
        return "unknown"
    else:
        return sex


rule pyclone_parse_tsv:
    input:
        vcf = "filtered_vcf/{patient}.vcf",
        cnv = lambda wildcards: pyclone_get_cnv_files(wildcards.cnv_method, wildcards.patient)
    output: "pyclone_{cnv_method}/{patient}.pyvi-input.tsv"
    params:
        sample_ids = lambda wildcards: patients_pyclone[wildcards.patient]["tumor_samples"],
        sex = lambda wildcards: get_patient_sex(wildcards.patient)
    conda: "clonalityParsers"
    script:
        "../scripts/prepare_pyclone-vi_input.R"


rule pyclone_fit:
    input: "pyclone_{cnv_method}/{patient}.pyvi-input.tsv"
    output: "pyclone_{cnv_method}/{patient}.h5"
    conda: "../envs/pyclonevi_env.yml"
    shell:
        "pyclone-vi fit -i {input} -o {output} -c 40 -d beta-binomial -r 50"


rule pyclone_write_results:
    input: "pyclone_{cnv_method}/{patient}.h5"
    output: "pyclone_{cnv_method}/{patient}.pyvi-output.tsv"
    conda: "../envs/pyclonevi_env.yml"
    shell:
        "pyclone-vi write-results-file -i {input} -o {output}"
