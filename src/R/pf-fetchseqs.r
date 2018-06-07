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

# Arguments for testing: opt <- list(options = list(sqlitedb = 'pf-fetchseqs.07.original.sqlite3', fetch = TRUE, verbose = TRUE, sourcedbs = 'refseq,pdb', faalevel='pfamily', faadir='.'))
SCRIPT_VERSION = "1.3.4"

# Get arguments
option_list = list(
  make_option(
    c("--fetch"), action="store_true", default=TRUE,
    help="Run fetching"
  ),
  make_option(
    c("--fetchedseqs"), type='character',
    help="Save all sequences to this file. Implies that *sequences will not be saved to the database*."
  ),
  make_option(
    c("--faadir"), type='character', default='.',
    help="Save output files resulting from faas over profiles at '--faalevel' to this directory, default %default."
  ),
  make_option(
    c("--faahmmcov"), type='double', default='0.0',
    help="Only save sequences to faa files if the HMMER hit covers at least this fraction of the hmm, default %default. In the db, proteins.hmmlen will be compared with hmm_profiles.plen."
  ),
  make_option(
    c("--faalevel"), type='character', 
    help="If set to 'psuperfamily', 'pfamily' or 'pclass', this option will make the program write individual faa files with sequences for each entry at the particular hierarchy level of profiles. Files will be written to the directory specied with '--faadir'."
  ),
  make_option(
    c("--skipfetch"), action="store_false", dest='fetch',
    help="Skip fetching, inserting only what's found in files on the command line"
  ),
  make_option(
    c("--only_prefetch"), action="store_true", default=FALSE, 
    help="Only run until not present accession numbers are written to file, default %default"
  ),
  make_option(
    c('--prefetch_accnos'), type='character', 
    help='File name to write accession numbers not present with sequences *before* attempting fetch'
  ),
  make_option(
    c('--postfetch_accnos'), type='character', 
    help='File name to write accession numbers not present with sequences *after* attempting fetch'
  ),
  make_option(
    c('--sourcedbs'), type='character', 
    help='Comma-separated list of databases to consider (e.g. "refseq", "pdb").'
  ),
  make_option(
    c('--sqlitedb'), type='character', 
    help='Name of sqlite file'
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
    usage = "%prog [options] file0.faa|file0.tsv|file0.tsv.gz ... filen.faa|filen.tsv|filen.tsv.gz", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION[1], stdout())
  quit('no')
}

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(stringr))

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

if ( length(opt$options$sourcedbs) ) {
  accessions <- accessions %>%
    semi_join(
      db %>% tbl('accessions') %>% collect() %>% filter(db %in% str_split(opt$options$sourcedbs, ',')[[1]]) %>%
        distinct(accto),
      by = c('accno' = 'accto')
    )
  logmsg(sprintf("Subset accessions to those belonging to %s, %d remaining", opt$options$sourcedbs, accessions %>% nrow()))
}

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
  } else if ( grepl('\\.feather$', fn) ) {
    logmsg(sprintf("Reading %s as feather", fn))
    d <- read_feather(fn)
  } else {
    logmsg(sprintf("Don't know how to handle %s", fn))
  }
  return(d)
}

# Loop over files on the command line
new_sequences <- tibble(accno = character(), sequence = character())
for ( f in opt$args ) {
  new_sequences <- new_sequences %>% union(handle_input(f) %>% anti_join(new_sequences, by = 'accno'))
}

# Add sequences from new_sequences that are not already in sequences but are in accessions
sequences <- sequences %>%
  union(
    new_sequences %>% anti_join(sequences, by = 'accno') %>%
      semi_join(accessions, by = 'accno')
  )
logmsg(sprintf("Added sequences from command line files, now %d sequences", sequences %>% nrow()))

acctofetch <- accessions %>% anti_join(sequences, by = 'accno') 
if ( length(opt$options$prefetch_accnos) > 0 ) {
  acctofetch %>% arrange(accno) %>% write_tsv(opt$options$prefetch_accnos)
  logmsg(sprintf("Saved accessions to %s", opt$options$prefetch_accnos))
}

if ( opt$options$only_prefetch ) {
  logmsg("Not continuing with fetch due to only_prefetch flag set")
  quit('no')
}

# Fetch sequences that are not yet in the sequences table
faafn <- tempfile(pattern = 'pf-fetchseqs.', tmpdir = '/tmp', fileext = '.faa')
errfn <- tempfile(pattern = 'pf-fetchseqs.', tmpdir = '/tmp', fileext = '.err')

# Function that fetches the sequence for an accno and appends to a file
fetch_seq <- function(accno, faafile, errfile) {
  system(sprintf("efetch -db protein -id %s -format fasta >> %s 2>>%s", accno, faafile, errfile))
}

if ( opt$options$fetch & acctofetch %>% nrow() > 0 ) {
  logmsg(sprintf("Fetching %d fasta formated sequences to %s, stderr to %s", acctofetch %>% nrow(), faafn, errfn))
  acctofetch %>% pull(accno) %>% walk(fetch_seq, faafn, errfn)

  logmsg("Done fetching, reading fasta file and updating/creating sequences table")

  newseqs <- readAAStringSet(faafn)
  sequences <- sequences %>% dplyr::union(
    tibble(accno = sub(' .*', '', names(newseqs)), sequence = as.character(newseqs))
  )
} else {
  logmsg("Skipping fetch")
}

if ( length(opt$options$fetchedseqs) > 0 ) {
  logmsg(sprintf("Saving table with %d sequences to %s", sequences %>% nrow(), opt$options$fetchedseqs))
  if ( grepl('\\.feather', opt$options$fetchedseqs) ) {
    sequences %>% arrange(accno) %>% write_feather(opt$options$fetchedseqs)
  } else {
    sequences %>% arrange(accno) %>% write_tsv(opt$options$fetchedseqs)
  }
} else {
  logmsg(sprintf("Inserting new table with %d sequences", sequences %>% nrow()))
  db %>% copy_to(sequences, 'sequences', temporary = FALSE, overwrite = TRUE)

  remaining <- db %>% tbl('tblout') %>% distinct(accno) %>%
    anti_join(db %>% tbl('sequences') %>% distinct(accno), by = 'accno') %>%
    collect()

  logmsg(sprintf("After insertion %d accessions remain without sequence", remaining %>% nrow()))
}
if ( length(opt$options$faalevel) > 0 ) {
  logmsg(sprintf("Writing faa files, one per %s, to %s", opt$options$faalevel, opt$options$faadir))
  if ( ! dir.exists(opt$options$faadir) ) dir.create(opt$options$faadir)
  for ( 
    ptaxon in db %>% tbl('hmm_profiles') %>% 
      transmute(pt = !! rlang::sym(opt$options$faalevel)) %>% 
      distinct() %>% filter(!is.na(pt)) %>% collect() %>% pull(pt)
  ) {
    f <- sprintf("%s/%s.faa", opt$options$faadir, ptaxon)
    s <- db %>% tbl('hmm_profiles') %>% filter(!! rlang::sym(opt$options$faalevel) == ptaxon) %>%
      inner_join(db %>% tbl('tblout'), by = 'profile') %>%
      inner_join(db %>% tbl('accessions') %>% transmute(accno = accto, taxon, db), by = 'accno') %>%
      inner_join(db %>% tbl('taxa'), by = 'taxon') %>%
      inner_join(db %>% tbl('proteins'), by = c('accno', 'profile')) %>%
      filter(hmmlen/as.integer(plen) >= opt$options$faahmmcov) %>%
      distinct(accno, tdomain, tphylum, tclass, psuperfamily, pfamily, pclass, psubclass, pgroup, taxon, db) %>% collect() %>%
      inner_join(sequences, by = 'accno') %>%
      mutate(name = sprintf("%s_%s_%s_%s_%s_%s_%s_%s_%s@%s_%s", tdomain, tphylum, tclass, gsub(' ', '_', taxon), psuperfamily, pfamily, pclass, psubclass, pgroup, db, accno)) %>%
      arrange(accno)

    ss <- AAStringSet(s$sequence)
    names(ss) <- s$name
    writeXStringSet(ss, f)
    logmsg(sprintf("\t%s: %d sequences saved to %s", ptaxon, nrow(s), f))
  }
} 

if ( length(opt$options$postfetch_accnos) > 0 ) remaining %>% arrange(accno) %>% write_tsv(opt$options$postfetch_accnos)

db %>% DBI::dbDisconnect()

logmsg("Done")
