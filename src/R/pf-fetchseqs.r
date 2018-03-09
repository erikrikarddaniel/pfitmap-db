#!/usr/bin/env Rscript

# pf-fetchseqs.r
#
# Fetches missing sequences. 
#
# The script needs an sqlite database with a tblout table. Each unique accession in the accno
# column will be fetched if not already found in a sequences table. The latter will be created
# if not already present with accno and sequence columns. Optionally, the script reads an amino
# acid fasta file given on the command line.
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

# Get arguments
option_list = list(
  make_option(
    c('--sqlitedb'), type='character', 
    help='Name of sqlite file'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] file.faa", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg(sprintf("Opening %s", opt$options$sqlitedb))

logmsg("Done")
