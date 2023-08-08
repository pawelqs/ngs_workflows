library(RColorBrewer)
library(ASCAT)
options(bitmapType='cairo')


ascat_dir = "ascat/"
dir.create("ascat/preproc", showWarnings = FALSE, recursive = TRUE)
tumor_bam = snakemake@input[["tumor_bam"]]
normal_bam = snakemake@input[["normal_bam"]]
tumor_name = snakemake@params[["tumor_name"]]
normal_name = snakemake@params[["normal_name"]]

BED_file = snakemake@params[["BED_file"]]
allelecounter_exe = snakemake@params[["allelecounter_path"]]
allele_prefix = snakemake@params[["allele_prefix"]]
loci_prefix = snakemake@params[["loci_prefix"]]
gc_file = snakemake@params[["gc_file"]]
rt_file = snakemake@params[["rt_file"]]
gender = snakemake@params[["gender"]]
purity = snakemake@params[["purity"]]
ploidy = snakemake@params[["ploidy"]]
genome_version = snakemake@params[["genome_version"]]
nthreads = snakemake@threads

tumourLogR_file = paste0(ascat_dir, "preproc/", tumor_name, ".LogR.txt")
tumourBAF_file = paste0(ascat_dir, "preproc/", tumor_name, ".BAF.txt")
normalLogR_file = paste0(ascat_dir, "preproc/", normal_name, ".LogRgermline_LogR.txt")
normalBAF_file = paste0(ascat_dir, "preproc/", normal_name, ".germline_LogR.txt")


#prepare from BAM files
ascat.prepareHTS(
  tumourseqfile = tumor_bam,
  normalseqfile = normal_bam,
  tumourname = tumor_name,
  normalname = normal_name,
  allelecounter_exe = allelecounter_exe,
  alleles.prefix = allele_prefix,
  loci.prefix = loci_prefix,
  gender = gender,
  genomeVersion = genome_version,
  nthreads = nthreads,
  # bed_file_arg = ,
  # chrom_names = c('21', '22'),
  tumourLogR_file = tumourLogR_file,
  tumourBAF_file = tumourBAF_file,
  normalLogR_file = normalLogR_file,
  normalBAF_file = normalBAF_file,
  seed = 42
)


#Load the data
ascat.bc = ascat.loadData(
  Tumor_LogR_file = tumourLogR_file,
  Tumor_BAF_file = tumourBAF_file,
  Germline_LogR_file = normalLogR_file,
  Germline_BAF_file = normalBAF_file,
  genomeVersion = genome_version,
  gender = gender
)

#Plot the raw data
prefix = paste0("ascat/", tumor_name)
ascat.plotRawData(ascat.bc, img.prefix = paste0(prefix, ".before_correction."))

#optional LogRCorrection
if(gc_file != "NULL") {
  # gc_input = paste0(normalizePath("$gc_input"), "/", "$gc_input", ".txt")
  
  if(rt_file != "NULL"){
    # rt_input = paste0(normalizePath("$rt_input"), "/", "$rt_input", ".txt")
    ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_file, replictimingfile = rt_file)
    #Plot raw data after correction
    ascat.plotRawData(ascat.bc,
      img.dir = "ascat/",
      img.prefix = paste0(tumor_name, ".after_correction_gc_rt.")
    )
  }
  else {
    ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = gc_file, replictimingfile = rt_file)
    #Plot raw data after correction
    ascat.plotRawData(ascat.bc, 
      img.dir = "ascat/",
      img.prefix = paste0(tumor_name, ".after_correction_gc.")
    )
  }
}
saveRDS(ascat.bc, paste0(tumor_name, ".Rds"))
#Segment the data
ascat.bc = ascat.aspcf(ascat.bc, seed=42)

#Plot the segmented data
ascat.plotSegmentedData(
  ascat.bc,
  img.dir = "ascat/",
  img.prefix = paste0(tumor_name, "."),
)

#Run ASCAT to fit every tumor to a model, inferring ploidy, normal cell contamination, and discrete copy numbers
#If psi and rho are manually set:
if (!is.null(purity) && !is.null(ploidy)){
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=purity, psi_manual=ploidy)
} else if(!is.null(purity) && is.null(ploidy)){
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1, rho_manual=purity)
} else if(!is.null(ploidy) && is.null(purity)){
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1, psi_manual=ploidy)
} else {
  ascat.output <- ascat.runAscat(ascat.bc, gamma=1)
}

#Extract metrics from ASCAT profiles
QC = ascat.metrics(ascat.bc,ascat.output)

#Write out segmented regions (including regions with one copy of each allele)
write.table(ascat.output[["segments"]], file=paste0(prefix, ".segments.txt"), sep="\t", quote=F, row.names=F)

#Write out CNVs in bed format
cnvs=ascat.output[["segments"]][2:6]
write.table(cnvs, file=paste0(prefix, ".cnvs.txt"), sep="\t", quote=F, row.names=F, col.names=T)

#Write out purity and ploidy info
summary <- tryCatch({
  matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]), ncol=2, byrow=TRUE)}, error = function(err) {
    # error handler picks up where error was generated
    print(paste("Could not find optimal solution:  ",err))
    return(matrix(c(0,0),nrow=1,ncol=2,byrow = TRUE))
  }
)
colnames(summary) <- c("AberrantCellFraction","Ploidy")
write.table(summary, file=paste0(prefix, ".purityploidy.txt"), sep="\t", quote=F, row.names=F, col.names=T)

write.table(QC, file=paste0(prefix, ".metrics.txt"), sep="\t", quote=F, row.names=F)

# # version export
# # f <- file("versions.yml","w")
# # alleleCounter_version = system(paste("alleleCounter --version"), intern = T)
# # ascat_version = sessionInfo()\$otherPkgs\$ASCAT\$Version
# # writeLines(paste0('"', "$task.process", '"', ":"), f)
# # writeLines(paste("    alleleCounter:", alleleCounter_version), f)
# # writeLines(paste("    ascat:", ascat_version), f)
# # close(f)

f <- file(paste0("ascat/", tumor_name, ".done"), "w")
writeLines("Done", f)
close(f)


sessionInfo()