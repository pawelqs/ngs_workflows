library(ASCAT)

ascat.prepareHTS(
  tumourseqfile = "../data/MK_OPUS_WES/G35_L1.dedup.recal.bam ",
  normalseqfile = "../data/MK_OPUS_WES/G35_C.dedup.recal.bam ",
  tumourname = "G35_L1",
  normalname = "G35_C",
  BED_file = "/scratch/scratch-hdd/pkus/resources/agilent/Agilent_SureSelect_v6r2_S07604514_Covered_hg38.bed",
  allelecounter_exe = "alleleCounter",
  alleles.prefix = "/scratch/scratch-hdd/pkus/resources/ASCAT/hg38/alleles/G1000_alleles_hg38_chr",
  loci.prefix = "/scratch/scratch-hdd/pkus/resources/ASCAT/hg38/loci/G1000_loci_hg38_chr",
  gender = "XX",
  genomeVersion = "hg38",
  nthreads = 8,
  tumourLogR_file = "G35_L1.LogR.txt",
  tumourBAF_file = "G35_L1.BAF.txt",
  normalLogR_file = "G35.germline_LogR.txt",
  normalBAF_file = "G35.germline_BAF.txt")