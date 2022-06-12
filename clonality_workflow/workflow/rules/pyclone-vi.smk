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



# def pyclone_write_samples_info(out_file, samples, tsv_files, purities):
#     with open(out_file, 'a') as f:
#         f.write(" ".join(samples) + "\n")
#         f.write(" ".join(tsv_files) + "\n")
#         f.write(" ".join([str(p) for p in purities]) + "\n")

# def pyclone_get_facets_purity(facets_out):
#     purity = pd.read_csv(facets_out).loc[1, "Purity"]
#     if np.isnan(purity):
#         purity = 1
#         print("Purity NaN, set to 1")
#     else:
#         print("Purity: %f" % purity)
#     return str(purity)

# rule pyclone_setup_facets:
#     input: 
#         purity_files = lambda wildcards: ["facets/%s.csv" % sample 
#             for sample in patients_pyclone[wildcards.patient]["tumor_samples"]],
#         tsv_files = lambda wildcards: ["pyclone_facets/%s/input/%s.tsv" % (wildcards.patient, sample) 
#             for sample in patients_pyclone[wildcards.patient]["tumor_samples"]]
#     output: "pyclone_facets/{patient}/input/samples_info.txt"
#     run:
#         samples = patients_pyclone[wildcards.patient]["tumor_samples"]
#         purities = [pyclone_get_facets_purity(str(facets_out)) for facets_out in input.purity_files]
#         pyclone_write_samples_info(str(output), samples, input.tsv_files, purities)


# def pyclone_get_titan_purity(params_file):
#     params = open(params_file).read().split("\n")
#     normal_cont_line = [line for line in params if "Normal contamination estimate" in line]
#     normal_contamination = float(normal_cont_line[0].split("\t")[1])
#     if normal_contamination > 0.9:
#         purity = 1
#         print("Estimated purity below 0.1, set to 1 as inaccurrate")
#     else:
#         purity = 1 - normal_contamination
#         print("Purity: %f" % purity)
#     return str(purity)

# rule pyclone_setup_titan:
#     input: 
#         purity_files = lambda wildcards: ["titan_optimal/%s/%s.params.txt" % (sample, sample) 
#             for sample in patients_pyclone[wildcards.patient]["tumor_samples"]],
#         tsv_files = lambda wildcards: ["pyclone_titan/%s/input/%s.tsv" % (wildcards.patient, sample) 
#             for sample in patients_pyclone[wildcards.patient]["tumor_samples"]]
#     output: "pyclone_titan/{patient}/input/samples_info.txt"
#     run:
#         samples = patients_pyclone[wildcards.patient]["tumor_samples"]
#         purities = [pyclone_get_titan_purity(str(pfile)) for pfile in input.purity_files]
#         pyclone_write_samples_info(str(output), samples, input.tsv_files, purities)


# rule pyclone_setup_analysis:
#     input: "pyclone_{cnv_method}/{patient}/input/samples_info.txt"
#     output: "pyclone_{cnv_method}/{patient}/config.yaml"
#     conda: "../envs/pyclone_env.yml"
#     shell: 
#         """
#         samples=`sed '1q;d' {input}`
#         cnv=`sed '2q;d' {input}`
#         purities=`sed '3q;d' {input}`
        
#         echo $samples
#         echo $cnv
#         echo $purities

#         PyClone setup_analysis \\
#             --samples $samples \\
#             --in_files $cnv \\
#             --tumour_contents $purities \\
#             --working_dir results/pyclone_{wildcards.cnv_method}/{wildcards.patient}
        
#         sed -i 's/.nan/null/g' {output}
#         """


# rule pyclone_run_analysis:
#     input: "pyclone_{cnv_method}/{patient}/config.yaml"
#     output: "pyclone_{cnv_method}/{patient}/trace/labels.tsv.bz2"
#     conda: "../envs/pyclone_env.yml"
#     shell: "PyClone run_analysis --config_file {input} --seed 4"



# rule pyclone_copy_results:
#     input: 
#         cluster_table = "pyclone_{cnv_method}/{patient}/table_loci.tsv",
#         loci_table = "pyclone_{cnv_method}/{patient}/table_cluster.tsv"
#     output: 
#         cluster_table = "pyclone_{cnv_method}/{patient}.table_loci.tsv",
#         loci_table = "pyclone_{cnv_method}/{patient}.table_cluster.tsv"
#     shell: 
#         """
#         cp {input.cluster_table} {output.cluster_table}
#         cp {input.loci_table} {output.loci_table} 
#         """
