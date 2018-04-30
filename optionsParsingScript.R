option_list <- list(
  make_option(c("-o", "--output"), type="character", default="simulatedData/analyticCase/", 
              help="Name of the output RData file (can contain path)"),
  make_option(c("-s", "--seed"), type="integer", default=0, 
              help="Seed for RNG")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
seed <- opt$seed
SAVEPATH <- opt$output