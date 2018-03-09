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
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(purrr))

# Arguments for testing: opt <- list(options = list(sqlitedb = 'pf-fetchseqs.00.sqlite3', verbose = TRUE))

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

db <- DBI::dbConnect(RSQLite::SQLite(), dbname = opt$options$sqlitedb)

accessions <- db %>% tbl('tblout') %>% distinct(accno) %>% collect()

# Do we have a sequences table or not?
if ( 'sequences' %in% (db %>% DBI::dbListTables()) ) {
  sequences <- db %>% table('sequences') %>% collect()
} else {
  sequences <- tibble(accno = character(), sequence = character())
}

fetch_seq <- function(accno, filename) {
  system(sprintf("efetch -db protein -id %s -format fasta >> %s", accno, filename))
}

# Fetch sequences that are not yet in the sequences table
tmpfn <- tempfile()
logmsg(sprintf("Fetching fasta formated sequences to %s", tmpfn))
accessions %>% anti_join(sequences, by = 'accno') %>% pull(accno) %>%
  walk(fetch_seq, tmpfn)

logmsg("Done fetching, reading fasta file and updating/creating sequences table")

newseqs <- readAAStringSet(tmpfn)
sequences <- sequences %>% dplyr::union(
  tibble(accno = sub(' .*', '', names(newseqs)), sequence = as.character(newseqs))
)
logmsg(sprintf("Inserting new table with %d sequences", sequences %>% nrow()))

db %>% copy_to(sequences, 'sequences', temporary = FALSE, overwrite = TRUE)

db %>% DBI::dbDisconnect()

logmsg("Done")