#!/usr/bin/env Rscript

# pf-fasta-unique-taxon-protein.r
#
# Reads a fasta file, the hmm_profiles, accessions, protein and taxa pfitmap
# feather tables (see option --featherprefix, and outputs a fasta file with one
# sequence per taxon and protein rank (options --trank and --prank
# respectively; currently only species and subclass).
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(feather))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))

SCRIPT_VERSION = "0.9.2"

# Testing arguments: opt <- list('options' = list('featherprefix' = 'pf-db2feather.01', 'prank' = 'psubclass', 'trank' = 'tspecies'), args = c('pf-fasta-unique-taxon-protein.00.faa'))
# Get arguments
option_list = list(
  make_option(
    c('--featherprefix'), default='', help='Prefix to feather files. The script assumes that <prefix>.n.feather exists for n in "accessions", "hmm_profiles", "proteins" and "rank". No default.'
  ),
  make_option(
    c('--prank'), default='psubclass', help='Protein rank to group on, i.e. output one accession per this and the chosen taxon rank (--trank), default %default. NOTE: Currently no other ranks are implemented.'
  ),
  make_option(
    c('--trank'), default='tspecies', help='Taxon rank to group on, i.e. output one accession per this and the chosen protein rank (--prank), default %default. NOTE: Currently no other ranks are implemented.'
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
    usage = '%prog [options] sequences.faa\n\n\tSubsets a fasta file only contain one sequence per protein and taxon rank.\n\n\tThe fasta file description lines must contain accessions like "@refseq_WP_123".',
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION, stdout())
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

logmsg(sprintf("Reading feather files with %s prefix", opt$options$featherprefix))

accessions   <- read_feather(sprintf("%s.accessions.feather",   opt$options$featherprefix))
hmm_profiles <- read_feather(sprintf("%s.hmm_profiles.feather", opt$options$featherprefix))
proteins     <- read_feather(sprintf("%s.proteins.feather",     opt$options$featherprefix))
taxa         <- read_feather(sprintf("%s.taxa.feather",         opt$options$featherprefix))

logmsg(sprintf("Reading fasta file %s", opt$args[1]))
f <- readAAStringSet(opt$args[1])
sequences <- tibble(name = names(f), seq = as.character(f)) %>%
  mutate(accno = gsub('.*@[a-z]+_([A-Za-z][A-Za-z_0-9.]+).*', '\\1', name))

logmsg(sprintf("Calculating one accession per %s and %s", opt$options$prank, opt$options$trank))
keep <- hmm_profiles %>%
  inner_join(proteins %>% select(profile, accno), by = 'profile') %>%
  inner_join(accessions %>% select(db, accno, accto, taxon), by = 'accno') %>%
  inner_join(taxa, by = 'taxon') %>%
  mutate(db = factor(db, levels = unique(c('refseq', order(unique(accessions$db)))), ordered = TRUE)) %>%
  arrange(db, accto) %>%
  group_by(psubclass, tspecies) %>% mutate(r = rank(accto)) %>% ungroup() %>% filter(r == 1) %>%
  transmute(accno = accto) %>% distinct() %>%
  arrange(accno) %>%
  inner_join(sequences, by = 'accno') %>%
  transmute(w = sprintf(">%s\n%s", name, seq))

logmsg(sprintf("Writing %d sequences to stdout", nrow(keep)))
write(keep$w, file = stdout())

logmsg("Done")
