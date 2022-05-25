import subprocess


def get_pwgs_info(patient):
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

patient_pwgs = {patient: get_pwgs_info(patient) for patient in patients}
file_types = ["mutass.zip", "muts.json.gz", "summ.json.gz"]


rule phylowgs_all:
    input: 
        expand(outdir + "/phylowgs_{cnv_method}/{patient}.Rds", cnv_method=["titan", "facets"], patient = patients)
        # expand(outdir + "/phylowgs_{cnv_method}/{patient}.{file_type}", 
        #        cnv_method=["titan", "facets"], patient = patients, file_type = file_types),
        # expand(outdir + "/phylowgs_{cnv_method}/{patient}/trees.zip", cnv_method=["titan", "facets"], patient=patients)
        # Only parse the files, stop before the multievolve:
        # expand(outdir + "/phylowgs_input/{patient}.{cnv_method}.params.json", patient=patients, cnv_method=["titan", "facets"])


rule prepare_avcf:
    input: outdir + "/filtered_vcf/{patient}.vcf"
    output: outdir + "/phylowgs_input/{patient}_{sample}.vcf"
    params:
        vcf_col = "{patient}_{sample}",
        temp_vcf = outdir + "/phylowgs_input/{patient}_{sample}.temp.vcf",
    shell: 
        """
        gatk SelectVariants \\
            --sample-name {params.vcf_col} \\
            -V {input} \\
            -O {params.temp_vcf}
        grep '^#' {params.temp_vcf} > {output}
        grep -v  '^#'  {params.temp_vcf} | awk '$10 !~ /\.\/\./ && $9 ~ /:AD:/' >> {output}
        rm {params.temp_vcf}
        """


rule parse_facets:
    input: outdir + "/facets/{sample}.csv"
    output: outdir + "/phylowgs_input/{sample}.cnvs.facets.txt"
    run:
        df = pd.read_csv(str(input))
        purity = df.loc[1, "Purity"]
        print("Purity: %f" % purity)
        subprocess.call([
            "python3", "workflow/scripts/parse_cnvs_facets_extension.py",
            "-f", "facets",
            "-c", str(purity),
            "--cnv-output", str(output),
            str(input)
        ])


rule parse_titan:
    input: 
        segments = outdir + "/titan_optimal/{sample}/{sample}.segs.tsv",
        params = outdir + "/titan_optimal/{sample}/{sample}.params.txt"
    output: outdir + "/phylowgs_input/{sample}.cnvs.titan.txt"
    run:
        params = open(str(input.params)).read().split("\n")
        normal_cont_line = [line for line in params if "Normal contamination estimate" in line]
        print(normal_cont_line)
        normal_contamination = float(normal_cont_line[0].split("\t")[1])
        purity = 1 - normal_contamination
        print("Purity: %f" % purity)
        subprocess.call([
            "python3", "workflow/scripts/parse_cnvs_facets_extension.py",
            "-f", "titan",
            "-c", str(purity),
            "--cnv-output", str(output),
            str(input.segments)
        ])


rule pwgs_create_inputs:
    input:
        vcfs = lambda wildcards: [outdir + "/phylowgs_input/%s.vcf" % sid 
                                  for sid in patient_pwgs[wildcards.patient]["tumor_samples"]],
        cnvs = lambda wildcards: [outdir + "/phylowgs_input/%s.cnvs.%s.txt" % (sample_id, wildcards.cnv_method) 
                                  for sample_id in patient_pwgs[wildcards.patient]["tumor_samples"]]
    output:
        ssms = outdir + "/phylowgs_input/{patient}.{cnv_method}.ssm_data.txt",
        cnvs = outdir + "/phylowgs_input/{patient}.{cnv_method}.cnv_data.txt",
        json = outdir + "/phylowgs_input/{patient}.{cnv_method}.params.json"
    params:
        samples = lambda wildcards: patient_pwgs[wildcards.patient]["tumor_samples"],
        sex = lambda wildcards: patient_pwgs[wildcards.patient]["sex"],           # "male" / "female"
        out_dir = outdir + "/{patient}/",
        cnvs = lambda wildcards, input: "".join(["--cnvs %s=%s " % (sample, cnv_file) 
                                        for sample, cnv_file 
                                        in zip(patient_pwgs[wildcards.patient]["tumor_samples"], input.cnvs)]),
        vcf_types = lambda wildcards, input: "".join(["--vcf-type %s=mutect_smchet " % sample 
                                             for sample 
                                             in patient_pwgs[wildcards.patient]["tumor_samples"]]),
        vcfs = lambda wildcards, input: "".join(["%s=%s " % (sample, vcf_file) 
                                        for sample, vcf_file 
                                        in zip(patient_pwgs[wildcards.patient]["tumor_samples"], input.vcfs)]),
        subset_size = 2000
    conda: "../envs/phylowgs_env.yml"
    shell: 
        """
        python2 workflow/scripts/create_phylowgs_inputs.py \\
            --sex {params.sex} --regions all  -s {params.subset_size} \\
            --output-cnvs {output.cnvs} --output-variants {output.ssms} --output-params {output.json} \\
            {params.cnvs} \\
            {params.vcf_types} \\
            {params.vcfs}
        """


rule pwgs_multievolve:
    input:
        ssms = outdir + "/phylowgs_input/{patient}.{cnv_method}.ssm_data.txt",
        cnvs = outdir + "/phylowgs_input/{patient}.{cnv_method}.cnv_data.txt"
    output: outdir + "/phylowgs_{cnv_method}/{patient}/trees.zip"
    threads: 20
    params:
        out_dir = outdir + "/phylowgs_{cnv_method}/{patient}"
    conda: "../envs/phylowgs_env.yml"
    shell: 
        """
        srun python2 ~/programs/phylowgs/multievolve.py \\
            --num-chains {threads} \\
            --ssms {input.ssms} --cnvs {input.cnvs} \\
            -O {params.out_dir}
        """


rule pwgs_write_results:
    input: outdir + "/phylowgs_{cnv_method}/{patient}/trees.zip"
    output:
        mutass = outdir + "/phylowgs_{cnv_method}/{patient}.mutass.zip",
        muts_gz = outdir + "/phylowgs_{cnv_method}/{patient}.muts.json.gz",
        summ_gz = outdir + "/phylowgs_{cnv_method}/{patient}.summ.json.gz",
        muts = outdir + "/phylowgs_{cnv_method}/{patient}.muts.json",
        summ = outdir + "/phylowgs_{cnv_method}/{patient}.summ.json"
    conda: "../envs/phylowgs_env.yml"
    shell: 
        """
        python2 ~/programs/phylowgs/write_results.py \\
            --include-ssm-names \\
            --include-multiprimary --max-multiprimary 1 \\
            {wildcards.patient} {input} {output.summ_gz} {output.muts_gz} {output.mutass}
        gunzip -c {output.summ_gz} > {output.summ}
        gunzip -c {output.muts_gz} > {output.muts}
        """


rule pwgs_extract_best_result:
    input:
        mutass = outdir + "/phylowgs_{cnv_method}/{patient}.mutass.zip",
        muts = outdir + "/phylowgs_{cnv_method}/{patient}.muts.json",
        summ = outdir + "/phylowgs_{cnv_method}/{patient}.summ.json"
    output:
        rds_file = outdir + "/phylowgs_{cnv_method}/{patient}.Rds",
        tree_file = outdir + "/phylowgs_{cnv_method}/{patient}.tree.tsv",
        clusters_file = outdir + "/phylowgs_{cnv_method}/{patient}.clusters.tsv",
        mutations_file = outdir + "/phylowgs_{cnv_method}/{patient}.mutations.tsv"
    conda: "../envs/genomicR_env.yml"
    params:
        out_dir = outdir + "/phylowgs_{cnv_method}",
        prefix = "{patient}"
    shell:
        """
        Rscript workflow/scripts/phylowgs_extract_best.R \\
            --summ_json_file {input.summ} \\
            --muts_json_file {input.muts} \\
            --mutass_zip_file {input.mutass} \\
            --out_dir {params.out_dir} \\
            --prefix {params.prefix}
        """
