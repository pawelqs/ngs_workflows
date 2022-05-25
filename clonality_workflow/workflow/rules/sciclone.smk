
def get_sciclone_info(patient):
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

patient_sciclone = {patient: get_sciclone_info(patient) for patient in patients}


rule sciclone_all:
    input: expand(outdir + "/sciclone_facets/{patient}.sciclone.tsv", patient = patients)


rule run_sciclone_facets:
    input:
        vcf = outdir + "/filtered_vcf/{patient}.vcf",
        cnv = lambda wildcards: [outdir + "/facets/%s.csv" % sample
            for sample in patient_sciclone[wildcards.patient]["tumor_samples"]]
    output: outdir + "/sciclone_facets/{patient}.sciclone.tsv"
    params:
        sample_names = lambda wildcards: patient_sciclone[wildcards.patient]["tumor_samples"]
    conda: "../envs/sciclone_env.yml"
    shell:
        """
        Rscript workflow/scripts/run_sciclone.R \\
            --vcf {input.vcf} \\
            --cnv_files {input.cnv} \\
            --sample_names {params.sample_names} \\
            --cnv_type FACETS \\
            --prefix {wildcards.patient} \\
            --target_regions resources/agilent_sureselect_v6/Agilent_SureSelect_v6r2_S07604514_Covered_hg38.cnvkit.bed \\
            --outdir results/sciclone_facets
        """
