

samples = pd.read_csv(config["manifest"], sep = "\t")
mutect_files = samples.loc[samples.mutect2 != ".", ["patient_id", "mutect2"]]
mutect_files = dict(zip(mutect_files.patient_id, mutect_files.mutect2))

# rule all:
#     input: expand("bcftools/{patient}.sex", patient=patients)


rule guess_gender:
    input: lambda wildcards: mutect_files[wildcards.patient]
    output: "bcftools/{patient}.sex"
    shell: 		
        """
        bcftools +guess-ploidy -e 0.1 -t GT {input} > {output}
        """
