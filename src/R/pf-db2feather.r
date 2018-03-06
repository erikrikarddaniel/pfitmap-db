#!/usr/bin/env Rscript

# pf-db2feather.r
#
# Outputs selected pfitmap tables from an sqlite3 database file to individual
# feather files.
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
#suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(feather))

SCRIPT_VERSION = "0.9.0"

# Get arguments
option_list = list(
  make_option(
    c('--prefix'), default='', help='Set a prefix for generated output files'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] file.sqlite3",
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

if ( opt$options$prefix != '' ) opt$options$prefix = sprintf("%s.", opt$options$prefix)

logmsg(sprintf("pf-db2feather.r %s, connecting to %s", SCRIPT_VERSION, opt$args[1]))

db <- DBI::dbConnect(RSQLite::SQLite(), opt$args[1])

list('accessions', 'dbsources', 'hmm_profiles', 'domains', 'dupfree_proteins', 'hmm_profiles', 'taxa') %>%
  walk(
    function(t) {
      logmsg(sprintf("Writing %s", t))
      db %>% tbl(t) %>% collect() %>%
        write_feather(sprintf("%s%s.feather", opt$options$prefix, t))
    }
  )

db %>% DBI::dbDisconnect()

logmsg("Done")
