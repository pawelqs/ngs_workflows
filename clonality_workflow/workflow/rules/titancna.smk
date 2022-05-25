# Source: https://github.com/gavinha/TitanCNA/tree/master/scripts/snakemake

import pandas as pd
import os.path


TITAN_CLUST = {1:[1], 2:[1,2], 3:[1,2,3], 4:[1,2,3,4], 5:[1,2,3,4,5], 6:[1,2,3,4,5,6], 7:[1,2,3,4,5,6,7], 8:[1,2,3,4,5,6,7,8], 9:[1,2,3,4,5,6,7,8,9], 10:[1,2,3,4,5,6,7,8,9,10]}
TITAN_PLOIDY = {2:[2], 3:[2,3], 4:[2,3,4]}


###########################################################
##################### TitanCNA part
###########################################################

# rule all:
    # input: 
        # expand("titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
               # tumor = tumor_samples, 
               # clustNum = TITAN_CLUST[config["TitanCNA_maxNumClonalClusters"]], 
               # ploidy=TITAN_PLOIDY[config["TitanCNA_maxPloidy"]]),
        # expand("titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt",
               # tumor = tumor_samples, 
               # clustNum = TITAN_CLUST[config["TitanCNA_maxNumClonalClusters"]], 
               # ploidy = TITAN_PLOIDY[config["TitanCNA_maxPloidy"]]),
        # expand("titan/hmm/titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt",
               # tumor = tumor_samples, 
               # clustNum = TITAN_CLUST[config["TitanCNA_maxNumClonalClusters"]], 
               # ploidy = TITAN_PLOIDY[config["TitanCNA_maxPloidy"]]),
        # "titan/hmm/optimalClusterSolution.txt",
        # "titan/hmm/optimalClusterSolution/"


### Give here 10 cpus!
rule runTitanCNA:
    input:
        alleleCounts = "titan/tumCounts/{tumor}.tumCounts.txt",
        corrDepth = "ichorCNA/{tumor}/{tumor}.correctedDepth.txt"        
    output:        
        titan = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
        param = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
        segTxt = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt",
        seg = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.seg"
    conda: "../envs/titancna_env.yml"
    params:
        outRoot = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}/",
        titanRscript = config["TitanCNA_rscript"],
        libdir = config["TitanCNA_libdir"],
        numCores = config["TitanCNA_numCores"],
        normal = config["TitanCNA_normalInit"],
        chrs = config["TitanCNA_chrs"],
        sex = lambda wildcards: sex[wildcards.tumor],
        genomeStyle = config["genomeStyle"],
        genomeBuild = config["genomeBuild"],
        cytobandFile = config["cytobandFile"],
        estimatePloidy = config["TitanCNA_estimatePloidy"],
        estimateClonality = config["TitanCNA_estimateClonality"],
        estimateNormal = config["TitanCNA_estimateNormal"],
        centromere = config["centromere"],
        alphaK = config["TitanCNA_alphaK"],
        #alphaR = config["TitanCNA_alphaR"],
        #alleleModel = config["TitanCNA_alleleModel"],
        txnExpLen = config["TitanCNA_txnExpLen"],
        plotYlim = config["TitanCNA_plotYlim"]
    log: "outputs/titan.hmm.titanCNA_ploidy{ploidy}.{tumor}_cluster{clustNum}.log"
    threads: 10
    shell:
        """
        Rscript {params.titanRscript} \\
            --hetFile {input.alleleCounts} \\
            --cnFile {input.corrDepth} \\
            --outFile {output.titan} \\
            --outSeg {output.segTxt} \\
            --outParam {output.param} \\
            --outIGV {output.seg} \\
            --outPlotDir {params.outRoot} \\
            --libdir {params.libdir} \\
            --id {wildcards.tumor} \\
            --numClusters {wildcards.clustNum} \\
            --numCores {threads} \\
            --normal_0 {params.normal} \\
            --ploidy_0 {wildcards.ploidy} \\
            --genomeStyle {params.genomeStyle} \\
            --genomeBuild {params.genomeBuild} \\
            --cytobandFile {params.cytobandFile} \\
            --chrs \"{params.chrs}\" \\
            --gender {params.sex} \\
            --estimateNormal {params.estimateNormal} \\
            --estimatePloidy {params.estimatePloidy} \\
            --estimateClonality {params.estimateClonality}  \\
            --centromere {params.centromere} \\
            --alphaK {params.alphaK} \\
            --txnExpLen {params.txnExpLen} \\
            --plotYlim \"{params.plotYlim}\" \\
            > {log} 2> {log}
           """


rule combineTitanAndIchorCNA:
    input:
        titanSeg = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.segs.txt", 
        titanBin = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.txt",
        titanParam = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.params.txt",
        ichorSeg = "ichorCNA/{tumor}/{tumor}.seg.txt",
        ichorBin = "ichorCNA/{tumor}/{tumor}.cna.seg",
        ichorParam = "ichorCNA/{tumor}/{tumor}.params.txt"
    output:
        segFile = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.seg.txt",
        binFile = "titan/{tumor}/hmm.titanCNA_ploidy{ploidy}/{tumor}_cluster{clustNum}.titan.ichor.cna.txt",
    conda: "../envs/titancna_env.yml"
    params:
        combineScript = config["TitanCNA_combineTitanIchorCNA"],
        libdir = config["TitanCNA_libdir"],
        centromere = config["centromere"],
        sex = lambda wildcards: sex[wildcards.tumor],
        mergeIchorHOMD = config["mergeIchorHOMD"]
    log: "outputs/titan.hmm.titanCNA_ploidy{ploidy}.{tumor}_cluster{clustNum}.combineTitanIchorCNA.log"
    shell:
        """
        Rscript {params.combineScript} \\
            --libdir {params.libdir} \\
            --titanSeg {input.titanSeg} \\
            --titanBin {input.titanBin} \\
            --titanParam {input.titanParam} \\
            --ichorSeg {input.ichorSeg} \\
            --ichorBin {input.ichorBin} \\
            --ichorParam {input.ichorParam} \\
            --mergeIchorHOMD {params.mergeIchorHOMD} \\
            --sex {params.sex} \\
            --outSegFile {output.segFile} \\
            --outBinFile {output.binFile} \\
            --centromere {params.centromere} \\
            > {log} 2> {log}
        """            


rule selectSolution:
    input:
        #ploidyDirs=expand("titan/hmm/titanCNA_ploidy{ploidy}/", ploidy=TITAN_PLOIDY[config["TitanCNA_maxPloidy"]]),
        resultFiles=expand("titan/{{tumor}}/hmm.titanCNA_ploidy{ploidy}/{{tumor}}_cluster{clustNum}.titan.txt", 
                           #tumor = tumor_samples, 
                           clustNum = TITAN_CLUST[config["TitanCNA_maxNumClonalClusters"]], 
                           ploidy = TITAN_PLOIDY[config["TitanCNA_maxPloidy"]])
    output: "titan/{tumor}/optimalClusterSolution.txt"
    conda: "../envs/titancna_env.yml"
    params:
        solutionRscript = config["TitanCNA_selectSolutionRscript"],
        threshold = config["TitanCNA_solutionThreshold"]
    log: "outputs/titan.{tumor}.selectSolution.log"
    shell:
        """
        ploidyRun2=results/titan/{wildcards.tumor}/hmm.titanCNA_ploidy2/
        if [ -d results/titan/{wildcards.tumor}/hmm.titanCNA_ploidy3/ ]; then
            ploidyRun3=results/titan/{wildcards.tumor}/hmm.titanCNA_ploidy3/
        else
            ploidyRun3=NULL
        fi
        if [ -d results/titan/{wildcards.tumor}/hmm.titanCNA_ploidy4/ ]; then
            ploidyRun4=results/titan/{wildcards.tumor}/hmm.titanCNA_ploidy4/
        else
            ploidyRun4=NULL
        fi
        Rscript {params.solutionRscript} \\
            --ploidyRun2 $ploidyRun2 \\
            --ploidyRun3 $ploidyRun3 \\
            --ploidyRun4 $ploidyRun4 \\
            --threshold {params.threshold} \\
            --outFile {output} \\
            > {log} 2> {log}
        """


rule copyOptSolution:
    input: "titan/{tumor}/optimalClusterSolution.txt"
    output: 
        segs = "titan_optimal/{tumor}/{tumor}.segs.txt",
        params = "titan_optimal/{tumor}/{tumor}.params.txt",
    params:
        out_dir = "titan_optimal/{tumor}"
    log: "outputs/titan.{tumor}.copyOptSolution.log"
    shell:
        """
        curDir=`pwd`
        for i in `cut -f11 {input} | grep -v "path"`;
        do
            echo -e "Copying $curDir/${{i}} to {params.out_dir}/"
            cp -r ${{curDir}}/${{i}}* {params.out_dir}/
        done
        
        cp {params.out_dir}/{wildcards.tumor}_cluster*.segs.txt {params.out_dir}/{wildcards.tumor}.segs.txt
        cp {params.out_dir}/{wildcards.tumor}_cluster*.params.txt {params.out_dir}/{wildcards.tumor}.params.txt
        cp {params.out_dir}/{wildcards.tumor}_cluster*.RData {params.out_dir}/{wildcards.tumor}.RData
        cp {params.out_dir}/{wildcards.tumor}_cluster*.seg {params.out_dir}/{wildcards.tumor}.seg
        cp {params.out_dir}/{wildcards.tumor}_cluster*.titan.txt {params.out_dir}/{wildcards.tumor}.titan.txt
        """


rule titan_clear_segfile:
    input: "titan_optimal/{sample}/{sample}.segs.txt"
    output: "titan_optimal/{sample}/{sample}.segs.tsv"
    run:
        cnv = pd.read_csv(str(input), sep = "\t", dtype = str)
        cnv = cnv.rename(columns={
            "Start_Position.bp.": "Start_Position(bp)",
            "End_Position.bp.": "End_Position(bp)",
            "Cellular_Prevalence": "Clonal_Frequency"
        })
        cnv = cnv[cnv["Start_Position(bp)"] != cnv["End_Position(bp)"]]
        cnv["Chromosome"] = cnv.Chromosome.str.replace("chr", "")
        print(cnv)
        cnv.to_csv(str(output), sep = "\t", index=False)


# ###########################################################
# ##################### getAlleleCounts part
# ###########################################################

# ## USERS MUST MODIFY THIS TO CHOSE THE CHROMOSOME NAMING
# # use this line if using NCBI naming
# #CHRS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X']
# # use this line if using UCSC naming
CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']

# rule tumCounts:
    # input: 
        # expand("titan/tumCounts/{tumor}/{tumor}.tumCounts.chr{chr}.txt", tumor=config["pairings"], chr=CHRS),        
        # expand("titan/tumCounts/{tumor}.tumCounts.txt", tumor=config["pairings"])


rule getHETsites:
    #input: lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
    input: lambda wildcards: bams[paired_normals[wildcards.tumor]]
    output: "titan/hetPosns/{tumor}/{tumor}.chr{chr}.vcf"
    params:
        refFasta = config["refFasta"],
        snpDB = config["snpVCF"],
        samtoolsCmd = config["samTools"],
        bcftoolsCmd = config["bcfTools"]
    log: "outputs/titan.hetPosns.{tumor}.chr{chr}.log"
    conda: "../envs/titancna_env.yml"
    shell:
        """
        {params.samtoolsCmd} mpileup -uv -I -f {params.refFasta} -r {wildcards.chr} -l {params.snpDB} {input} | \\
            {params.bcftoolsCmd} call -v -c - | \\
            grep -e '0/1' -e '#' \\
            > {output} 2> {log}
        """


rule getAlleleCountsByChr:
    input:
        hetSites = "titan/hetPosns/{tumor}/{tumor}.chr{chr}.vcf",
        #tumBam=lambda wildcards: config["samples"][wildcards.tumor]
        tumBam = lambda wildcards: bams[wildcards.tumor]
    output: "titan/tumCounts/{tumor}/{tumor}.tumCounts.chr{chr}.txt"
    params:
        countScript = config["pyCountScript"],
        #pyEnv=config["pyEnv"],
        #refFasta=config["refFasta"],
        mapQ = config["map_quality"],
        baseQ = config["base_quality"],
        vcfQ = config["vcf_quality"]
    log: "outputs/titan.tumCounts.{tumor}.chr{chr}.log"
    conda: "../envs/titancna_env.yml"
    shell: "python {params.countScript} {wildcards.chr} {input.hetSites} {input.tumBam} {params.baseQ} {params.mapQ} {params.vcfQ} > {output} 2> {log}"


rule catAlleleCountFiles:
    input: expand("titan/tumCounts/{{tumor}}/{{tumor}}.tumCounts.chr{chr}.txt", chr=CHRS)
    output: "titan/tumCounts/{tumor}.tumCounts.txt"
    log: "outputs/titan.tumCounts.{tumor}.cat.log"
    shell: "cat {input} | grep -v Chr > {output} 2> {log}"
        

###########################################################
##################### ichorCNA part
###########################################################

# rule correct_depth:
    # input:
        # expand("ichorCNA/{tumor}/{tumor}.cna.seg", tumor=tumor_samples),
        # expand("ichorCNA/{tumor}/{tumor}.seg.txt", tumor=tumor_samples),
        # expand("ichorCNA/{tumor}/{tumor}.params.txt", tumor=tumor_samples),
        # expand("ichorCNA/{tumor}/{tumor}.correctedDepth.txt", tumor=tumor_samples),
        # expand("readDepth/{sample}.bin{binSize}.wig", sample=samples.sample_id, binSize=str(config["binSize"]))


rule read_counter:
    input: 
        #bam = lambda wildcards: config["samples"][wildcards.samples]
        bam = lambda wildcards: bams[wildcards.sample],
        bai = lambda wildcards: "%s.bai" % bams[wildcards.sample]
    output: "readDepth/{sample}.bin{binSize}.wig"
    conda: "../envs/titancna_env.yml"        
    params:
        readCounter = config["readCounterScript"],
        binSize = config["binSize"],
        qual = "20",
        chrs = config["chrs"],
        out_dir = "readDepth"
    resources:
        mem = 4
    log: "outputs/readDepth.{sample}.bin{binSize}.log"
    shell: 
        """
        mkdir -p {params.out_dir}
        {params.readCounter} {input.bam} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}
        """


rule ichorCNA:
    input:
        # tum = "readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
        # norm = lambda wildcards: "readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
        tum = "readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
        norm = lambda wildcards: "readDepth/" + paired_normals[wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
    output:
        corrDepth = "ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        param = "ichorCNA/{tumor}/{tumor}.params.txt",
        cna = "ichorCNA/{tumor}/{tumor}.cna.seg",
        segTxt = "ichorCNA/{tumor}/{tumor}.seg.txt",
        #seg = "ichorCNA/{tumor}/{tumor}.seg",
        #rdata = "ichorCNA/{tumor}/{tumor}.RData"
    conda: "../envs/titancna_env.yml"
    params:
        outDir = "ichorCNA/{tumor}/",
        rscript = config["ichorCNA_rscript"],
        libdir = config["ichorCNA_libdir"],
        id = "{tumor}",
        ploidy = config["ichorCNA_ploidy"],
        normal = config["ichorCNA_normal"],
        genomeStyle = config["genomeStyle"],
        genomeBuild = config["genomeBuild"],
        gcwig = config["ichorCNA_gcWig"],
        mapwig = config["ichorCNA_mapWig"],
        estimateNormal = config["ichorCNA_estimateNormal"],
        estimatePloidy = config["ichorCNA_estimatePloidy"],
        estimateClonality = config["ichorCNA_estimateClonality"],
        scStates = config["ichorCNA_scStates"],
        maxCN = config["ichorCNA_maxCN"],
        includeHOMD = config["ichorCNA_includeHOMD"],
        chrs = config["ichorCNA_chrs"],
        #chrTrain = config["ichorCNA_chrTrain"],
        centromere = config["centromere"],
        exons = config["ichorCNA_exons"],
        txnE = config["ichorCNA_txnE"],
        txnStrength = config["ichorCNA_txnStrength"],
        fracReadsChrYMale = "0.001",
        plotFileType = config["ichorCNA_plotFileType"],
        plotYlim = config["ichorCNA_plotYlim"]
    resources:
        mem = 4
    log: "outputs/ichorCNA.{tumor}.log"    
    shell:
        """
        Rscript {params.rscript} \\
            --libdir {params.libdir} \\
            --id {params.id} \\
            --WIG {input.tum} \\
            --gcWig {params.gcwig} \\
            --mapWig {params.mapwig} \\
            --NORMWIG {input.norm} \\
            --ploidy \"{params.ploidy}\" \\
            --normal \"{params.normal}\" \\
            --maxCN {params.maxCN} \\
            --includeHOMD {params.includeHOMD} \\
            --genomeStyle {params.genomeStyle} \\
            --genomeBuild {params.genomeBuild} \\
            --chrs \"{params.chrs}\" \\
            --estimateNormal {params.estimateNormal} \\
            --estimatePloidy {params.estimatePloidy} \\
            --estimateScPrevalence {params.estimateClonality} \\
            --scStates \"{params.scStates}\" \\
            --centromere {params.centromere} \\
            --exons.bed {params.exons} \\
            --txnE {params.txnE} \\
            --txnStrength {params.txnStrength} \\
            --fracReadsInChrYForMale {params.fracReadsChrYMale} \\
            --plotFileType {params.plotFileType} \\
            --plotYLim \"{params.plotYlim}\" \\
            --outDir {params.outDir} \\
            > {log} 2> {log}
        """

