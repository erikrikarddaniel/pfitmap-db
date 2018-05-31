#!/usr/bin/env Rscript

# pf-db2feather.r
#
# Outputs selected pfitmap tables from an sqlite3 database file to individual
# feather files.
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))

SCRIPT_VERSION = "1.0.3"

# Options for testing: opt <- list(options = list(verbose = TRUE, prefix='testing'), args = 'pf-classify.02.sqlite3')
# Get arguments
option_list = list(
  make_option(
    c('--dbs'), default='', help='Comma-separated list of databases to subset data to'
  ),
  make_option(
    c('--prefix'), default='', help='Set a prefix for generated output files'
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  ),
  make_option(
    c("-V", "--version"), action="store_true", default=FALSE, 
    help="Print program version and exit"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] file.sqlite3",
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION[1], stdout())
  quit('no')
}

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

dbs <- if(opt$options$dbs != '') str_split(opt$options$dbs, ',')[[1]] else  db %>% tbl('accessions') %>% distinct(db) %>% filter(!is.na(db)) %>% collect() %>% pull(db)

# All tables that are subset by db, must be taken care of separately
accessions <- db %>% tbl('accessions') %>% filter(db %in% dbs)
tt <- accessions %>% collect()
logmsg(sprintf("Writing accessions, %d rows", nrow(tt)))
tt %>% write_feather(sprintf("%saccessions.feather", opt$options$prefix))

tt <- db %>% tbl('taxa') %>% semi_join(accessions %>% distinct(taxon), by = 'taxon') %>% collect()
logmsg(sprintf("Writing taxa, %d rows", nrow(tt)))
tt %>% write_feather(sprintf("%staxa.feather", opt$options$prefix))

tt <- db %>% tbl('dbsources') %>% collect()
logmsg(sprintf("Writing dbsources, %d rows", nrow(tt)))
tt %>% write_feather(sprintf("%sdbsources.feather", opt$options$prefix))

tt <- db %>% tbl('hmm_profiles') %>% collect()
logmsg(sprintf("Writing hmm_profiles, %d rows", nrow(tt)))
tt %>% write_feather(sprintf("%shmm_profiles.feather", opt$options$prefix))

intersect(
  c('domains', 'dupfree_proteins', 'proteins', 'sequences', 'tblout', 'domtblout'),
  db %>% DBI::dbListTables()
) %>% walk(
    function(t) {
      tt <- db %>% tbl(t) %>% semi_join(accessions %>% transmute(accno = accto) %>% distinct(), by = 'accno') %>% collect()
      logmsg(sprintf("Writing %s, %d rows", t, nrow(tt)))
      tt %>% write_feather(sprintf("%s%s.feather", opt$options$prefix, t))
    }
  )

db %>% DBI::dbDisconnect()

logmsg("Done")
