testthat::test_that("vcftodb works", {
  
  db_path <- system.file("extdata","testdata", package = "SomaVarDBTools")
  prefix <- "test"
  vcf_list <- list.files(system.file("extdata","testdata", package = "SomaVarDBTools"), pattern = ".gz$", full.names = TRUE)

  future::plan("multicore")
  options(future.globals.maxSize = 10000*1024^2)
  options(mc.cores = 3)
  for (vcf in vcf_list){
    print(vcf)
    #SomaVarDBTools::buildDB_seqone(prefix = prefix, vcf_name = vcf, db_path = db_path)
    #SomaVarDBTools::buildDB_sarek_haplotypecaller(prefix = prefix, vcf_name = vcf, db_path = db_path)
    SomaVarDBTools::buildDB_sarek_mutect2(prefix = prefix, vcf_name = vcf, db_path = db_path)
  }

  SomaVarDBTools::addMDtodb(db_path = db_path, prefix = prefix, API_key = "l6ln_tsCbdZk8ICVpy6z86vhhb2rFUkFFTfV3wq7L-4", workers = 2)
  
  SomaVarDBTools::update_db_metadata(db_path = db_path, prefix = prefix)
  
  SomaVarDBTools::addSampleAttributes(db_path = db_path, prefix = prefix, matches_file = system.file("extdata","testdata/samples_run_match.tsv", package = "SomaVarDBTools"))

  SomaVarDBTools::compute_frequency(db_path = db_path, prefix = prefix ,attribute = "constant_group")
  
  SomaVarDBTools::CNVtodb_seqone(db_path = db_path, prefix = prefix, cnvfile_path = system.file("extdata","testdata/cnv.tsv", package = "SomaVarDBTools"))

  SomaVarDBTools::QCtodb_seqone(db_path = db_path, prefix = prefix, qcfile_path = system.file("extdata","testdata/QCs.csv", package = "SomaVarDBTools"))
  
  SomaVarDBTools::initcomment(db_path = db_path, prefix = prefix)
})
