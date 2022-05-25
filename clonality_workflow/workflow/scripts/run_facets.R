### https://github.com/aleighbrown/facets_snakemake/blob/master/run_facets.R

#load the libraries
library(facets)
library(data.table)

#argument should be the snppileup produced in the previous step and the desired output folder

set.seed(1234)

#critical value to call a change
CVAL = 150
#number of hetrozygous snps required to call minor copy number
NHET = snakemake@params[["nhet"]]
#the gzfile produced by snppileup
inputFilename = snakemake@input[[1]]
#the place where things should be written
outputFilename = snakemake@output[[1]]
#the directory where you want to write things
outputDir = snakemake@params[["outdir"]]
sampleName = snakemake@params[["sample_name"]]
#make sure the directory exists already, this will throw a warning if the output folder exists already
dir.create(outputDir)
print(inputFilename)
print(paste0("Reading file: ", inputFilename))
#read the Snp matrix, feed a file path of the snp pileup 
rcmat = readSnpMatrix(inputFilename)

#perform the preprocessing
print(paste("preprocessing file",inputFilename))
xx = preProcSample(rcmat, gbuild = "hg38")
#process the sample
print(paste("Processing file",inputFilename))

oo = procSample(xx,cval=CVAL)
print(oo$flags)
#we're using the emcncf2 function to allow the data to provide subclonal copy number calls as well
print(paste("Running fit",inputFilename))
fit = emcncf(oo, min.nhet = NHET)
#the data is store in cnfc
fit_table = as.data.table(fit$cncf)
fit_table[,"Purity" := fit$purity]
fit_table[,"Ploidy" := fit$ploidy]

fwrite(fit_table, outputFilename)
print(paste(inputFilename, "written"))

print("Drawing diagnostic plots")

print("Plotting Sample")
png(paste0(outputDir, sampleName, "_", CVAL, NHET, "_cnv.png"), units="px", width=1600, height=1600, res=300)
sname <- sprintf('%s; ploidy= %.2f; purity= %.2f', sampleName, fit$ploidy, fit$purity)
plotSample(x = oo, emfit = fit, sname = sname)

while(!is.null(dev.list())){dev.off()}
print("Ploting Spider")
png(paste0(outputDir, sampleName, "_", CVAL, NHET, "_spider.png"), width = 3.25, height = 3.25, units = "in", res=1200, pointsize = 4)
    par(mar = c(5, 5, 2, 2), xaxs = "i", yaxs = "i",cex.axis = 2,cex.lab  = 2)
logRlogORspider(oo$out,oo$dipLogR)

while(!is.null(dev.list())){dev.off()}


print(paste("All done with sample:",inputFilename))