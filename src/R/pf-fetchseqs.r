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
suppressPackageStartupMessages(library(readr))

# Arguments for testing: opt <- list(options = list(sqlitedb = 'pf-fetchseqs.02.original.sqlite3', verbose = TRUE))
SCRIPT_VERSION = "0.9"

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
logmsg(sprintf("pd-fetchseqs.r version %s, opening %s", SCRIPT_VERSION, opt$options$sqlitedb))

db <- DBI::dbConnect(RSQLite::SQLite(), dbname = opt$options$sqlitedb)

accessions <- db %>% tbl('tblout') %>% distinct(accno) %>% collect()
logmsg(sprintf("Read accessions, %d rows", accessions %>% nrow()))

# Do we have a sequences table or not?
if ( 'sequences' %in% (db %>% DBI::dbListTables()) ) {
  sequences <- db %>% tbl('sequences') %>% collect()
} else {
  sequences <- tibble(accno = character(), sequence = character())
}
logmsg(sprintf("Read sequences, %d rows", sequences %>% nrow()))

# Reads a tsv or fasta file with sequences and returns as a tibble
handle_input <- function(fn) {
  if ( grepl('\\.faa|\\.fasta$', fn) ) {
    logmsg(sprintf("Reading %s as fasta", fn))
    f <- readAAStringSet(fn)
    d <- tibble(accno = names(f), sequence = as.character(f))
  } else if ( grepl('\\.tsv(\\.gz)?$', fn) ) {
    logmsg(sprintf("Reading %s as tsv", fn))
    d <- read_tsv(fn, col_names = c('accno', 'sequence'), col_types = cols(.default = col_character())) %>%
      filter(accno != 'accno')
  }
  return(d)
}

# Loop over files on the command line
new_sequences <- tibble(accno = character(), sequence = character())
for ( f in opt$args ) {
  new_sequences <- new_sequences %>% union(handle_input(f))
}

# Add sequences from new_sequences that are not already in sequences but are in accessions
sequences <- sequences %>%
  union(
    new_sequences %>% anti_join(sequences, by = 'accno') %>%
      semi_join(accessions, by = 'accno')
  )
logmsg(sprintf("Added sequences from command line files, now %d sequences", sequences %>% nrow()))

# Function that fetches the sequence for an accno and appends to a file
fetch_seq <- function(accno, filename) {
  system(sprintf("efetch -db protein -id %s -format fasta >> %s 2>/dev/null", accno, filename))
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

remaining <- db %>% tbl('tblout') %>% distinct(accno) %>%
  anti_join(db %>% tbl('sequences') %>% distinct(accno), by = 'accno') %>%
  collect()

logmsg(sprintf("After insertion %d accessions remain without sequence", remaining %>% nrow()))

db %>% DBI::dbDisconnect()

logmsg("Done")
