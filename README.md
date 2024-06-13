# :clipboard: Requirements

R 4.2.0

# :hammer_and_wrench: Installation

`devtools::install_github("clbenoit/SomaVarDBTools")`

# :dna: Manage your genomic variations database

- `SomaVarDBTools::buildDB_seqone(prefix = prefix, db_path = db_path, vcf_name = "pathtovcf")` Will import genomic variations found on the specified vcf (seqone format) to the selected database. [More information abouth the vcf specifications](inst/extdata/testdata/README.md)

- `SomaVarDBTools::addMDtodb(db_path = db_path, prefix = prefix, API_key = "MD_API_key", workers = 2)` Will annotate variants in the selected database. A [Mobidetails](https://mobidetails.iurc.montp.inserm.fr/MD) account is required
  
- `SomaVarDBTools::addSampleAttributes(db_path = db_path, prefix = prefix, matches_file = system.file("extdata","testdata/samples_run_match.tsv", package = "SomaVarDBTools"))` Add samples attributes to the selected databases. [These attributes should be defined in a tsv file ](inst/extdata/testdata/samples_run_match.tsv)

- `SomaVarDBTools::compute_frequency(db_path = db_path, prefix = prefix ,attribute = "run")` Compute genomic variant frequencies accoding to a specific samples attribute. Previous Example : According to sample sequencing run belonging.

- `SomaVarDBTools::CNVtodb_seqone(db_path = db_path, prefix = prefix, cnvfile_path = system.file("extdata","testdata/cnv.tsv", package = "SomaVarDBTools"))` Will import copy number variations found on the specified cnv file [(seqone format)](inst/extdata/testdata/cnv.tsv) to the selected database. 

- `SomaVarDBTools::QCtodb_seqone(db_path = db_path, prefix = prefix, qcfile_path = system.file("extdata","testdata/QCs.csv", package = "SomaVarDBTools"))` Will import copy number variations found on the specified [QCs file](inst/extdata/testdata/QCs.tsv) to the selected database. 

# :soon: Incoming features

- Feed database with vcf produced by [sarek nf core pipeline](https://nf-co.re/sarek/3.2.3)
- Possibility for the user to save filters and parameters presets and choose a default one 

# :warning: Troubleshouting

