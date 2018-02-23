#!/usr/bin/env Rscript

# hmmsearch2classification
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

SCRIPT_VERSION = "1.1.2"

# Get arguments
option_list = list(
  make_option(
    c('--dbsource', default='', help='dbsource:name:version')
  ),
  make_option(
    c('--profilehierarchies', default='', help='tsv file with profile hiearchies')
  ),
  make_option(
    c('--singletable', default='', help='Write data in a single tsv format to this filename.')
  ),
  make_option(
    c('--sqlitedb', default='', help='Write data in a SQLite database with this filename.')
  ),
  make_option(
    c('--taxflat', default='', help='Name of NCBI taxon table in "taxflat" format (see https://github.com/erikrikarddaniel/taxdata2taxflat).')
  ),
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  )
)
opt = parse_args(
  OptionParser(
    usage = "%prog [options] file0.tblout .. filen.tblout file0.domtblout .. filen.domtblout", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( length(grep('sqlitedb', names(opt$options), value = TRUE)) > 0 ) {
  suppressPackageStartupMessages(library(dbplyr))
}

# Args list for testing:
# opt = list(args = c('hmmsearch2classification.00.d/GRX.ncbi_nr.test.domtblout', 'hmmsearch2classification.00.d/GRX.ncbi_nr.test.tblout', 'hmmsearch2classification.00.d/NrdAe.tblout','hmmsearch2classification.00.d/NrdAg.tblout','hmmsearch2classification.00.d/NrdAh.tblout','hmmsearch2classification.00.d/NrdAi.tblout','hmmsearch2classification.00.d/NrdAk.tblout','hmmsearch2classification.00.d/NrdAm.tblout','hmmsearch2classification.00.d/NrdAn.tblout','hmmsearch2classification.00.d/NrdAq.tblout','hmmsearch2classification.00.d/NrdA.tblout','hmmsearch2classification.00.d/NrdAz3.tblout','hmmsearch2classification.00.d/NrdAz4.tblout','hmmsearch2classification.00.d/NrdAz.tblout','hmmsearch2classification.00.d/NrdAe.domtblout','hmmsearch2classification.00.d/NrdAg.domtblout','hmmsearch2classification.00.d/NrdAh.domtblout','hmmsearch2classification.00.d/NrdAi.domtblout','hmmsearch2classification.00.d/NrdAk.domtblout','hmmsearch2classification.00.d/NrdAm.domtblout','hmmsearch2classification.00.d/NrdAn.domtblout','hmmsearch2classification.00.d/NrdAq.domtblout','hmmsearch2classification.00.d/NrdA.domtblout','hmmsearch2classification.00.d/NrdAz3.domtblout','hmmsearch2classification.00.d/NrdAz4.domtblout','hmmsearch2classification.00.d/NrdAz.domtblout'), options=list(verbose=T, singletable='test.out.tsv', profilehierarchies='hmmsearch2classification.00.phier.tsv', taxflat='hmmsearch2classification.taxflat.tsv', sqlitedb='testdb.sqlite3', dbsource='NCBI:NR:20180212'))

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg(sprintf("hmmsearch2classification.r version %s: Starting classification", SCRIPT_VERSION))

dbsource = strsplit(opt$options$dbsource, ':')[[1]]

logmsg(sprintf("Reading profile hierarchies from %s", opt$options$profilehierarchies))
hmm_profiles <- read_tsv(opt$options$profilehierarchies, col_types=cols(.default=col_character()))

logmsg(sprintf("Reading taxflat from %s", opt$options$taxflat))
taxflat <- read_tsv(opt$options$taxflat, col_types=cols(.default=col_character(), ncbi_taxon_id=col_integer())) %>%
  transmute(
    ncbi_taxon_id, taxon, trank = rank,
    tdomain       = superkingdom, tkingdom = kingdom,
    tphylum       = phylum,       tclass   = class,
    torder        = order,        tfamily  = family,
    tgenus        = genus,        tspecies = species
  )

# Delete duplicate taxon, rank combinations belonging in Eukaryota
taxflat <- taxflat %>%
  anti_join(
    taxflat %>% group_by(taxon, trank) %>% summarise(n = n()) %>% ungroup() %>% filter(n > 1) %>%
      inner_join(taxflat %>% filter(tdomain == 'Eukaryota'), by = c('taxon', 'trank')),
    by = c('ncbi_taxon_id')
  )

# We will populate two tables, one with the full results, one with accessions
tblout = tibble(
  accno = character(), profile = character(),
  evalue = double(), score = double(), bias = double()
)
accessions = tibble(accno = character(), accto = character())

# Read all the tblout files
for ( tbloutfile in grep('\\.tblout', opt$args, value=TRUE) ) {
  logmsg(sprintf("Reading %s", tbloutfile))
  t =  read_fwf(
    tbloutfile, fwf_cols(content = c(1, NA)), 
    col_types = cols(content = col_character()), 
    comment='#'
  ) %>% 
    separate(
      content, 
      c('accno', 't0', 'profile', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'), 
      '\\s+', 
      extra='merge',
      convert = T
    )
  tblout = union(tblout, t %>% select(accno, profile, evalue, score, bias))
  accessions = union(accessions, t %>% transmute(accno, accto = sprintf("%s %s", accno, rest)))
}

# Split the accto field
accessions = accessions %>% separate_rows(accto, sep = '\x01') %>% 
  mutate(
    taxon = ifelse(
      grepl('[^[]\\[(.*)\\]', accto), 
      sub('.*[^[]\\[(.*)\\].*', '\\1', accto),
      'unknown'
    ),
    accto = sub(' .*', '', accto)
  )

domtblout = tibble(
  accno = character(), tlen = integer(), profile = character(), qlen = integer(), i = integer(), n = integer(), 
  dom_c_evalue = double(), dom_i_evalue = double(), dom_score = double(), dom_bias = double(),
  hmm_from = integer(), hmm_to = integer(), ali_from = integer(), ali_to = integer(), 
  env_from = integer(), env_to = integer()
)

# Read all the domtblout files
for ( domtbloutfile in grep('\\.domtblout', opt$args, value=TRUE) ) {
  logmsg(sprintf("Reading %s", domtbloutfile))
  t = read_fwf(
    domtbloutfile, fwf_cols(content = c(1, NA)), 
    col_types = cols(content = col_character()), 
    comment='#'
  ) %>% 
    separate(
      content, 
      c(
        'accno', 't0', 'tlen', 'profile', 't1', 'qlen',  'evalue', 'score', 'bias', 'i', 'n', 
        'dom_c_evalue', 'dom_i_evalue', 'dom_score', 'dom_bias', 
        'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'acc', 'rest'
      ),
      '\\s+', 
      extra='merge',
      convert = T
    )
  
  domtblout = union(
    domtblout,
    t %>% select(
      accno, tlen, profile, qlen, i, n, dom_c_evalue, dom_i_evalue, dom_score, dom_bias,
      hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
    )
  )
}

# Calculate lengths:
# This is long-winded because we have three lengths (hmm, ali and env) which
# have to be calculated separately after getting rid of overlapping rows in
# domtblout.

# Define a temporary table that will be filled with lengths, minimum from and maximum
# to values.
lengths <- tibble(
  accno = character(), profile = character(), type = character(), val = integer()
)

# Function that joins all with n > 1 with the next row and not occurring in the second
# table in the call.
do_nextjoin <- function(dt, no) {
  dt %>% 
    filter(n > 1) %>%
    left_join(
      domt %>% transmute(accno, profile, i = i - 1, next_row = TRUE, next_from = from, next_to = to), 
      by = c('accno', 'profile', 'i')
    ) %>%
    replace_na(list('next_row' = FALSE)) %>%
    return()
}

# FOR EACH FROM-TO PAIR:
for ( fs in list(
  c('hmm_from', 'hmm_to', 'hmmlen'),
  c('ali_from', 'ali_to', 'alilen'),
  c('env_from', 'env_to', 'envlen')
)) {
  logmsg(sprintf("Calculating %s", fs[3]))

  # 1. Set from and to to current pair, delete shorter rows with the same start and 
  # calculate i (rownumber) and n (total domains) for each combination of accno and profile
  domt <- domtblout %>% transmute(accno, profile, from = .data[[fs[1]]], to = .data[[fs[2]]]) 
  domt <- domt %>% distinct() %>%
    semi_join(
      domt %>% group_by(accno, profile, from) %>% filter(to == max(to)) %>% ungroup,
      by = c('accno', 'profile', 'from', 'to')
    ) %>% 
    arrange(accno, profile, from, to) %>%
    group_by(accno, profile) %>% mutate(i = row_number(from), n = n()) %>% ungroup()

  # This is where non-overlapping rows will be stored
  nooverlaps <- tibble(
    accno = character(), profile = character(),
    from = integer(), to = integer()
  )

  while ( domt %>% nrow() > 0 ) {
    logmsg(sprintf("Working on overlaps, nrow: %d", domt %>% nrow()))
    
    # 2. Move rows with n == 1 to nooverlaps
    nooverlaps <- nooverlaps %>%
      union(domt %>% filter(n == 1) %>% select(accno, profile, from, to))
    domt <- domt %>% filter(n > 1)
    
    nextjoin <- do_nextjoin(domt, nooverlaps)

    # As a debug aid, save the domt and nextjoin data sets
    write_tsv(domt, 'domt.tsv.gz')
    write_tsv(nextjoin, 'nextjoin.tsv.gz')
  
    # 3. Move rows that do not overlap with the next row
    nooverlaps <- nooverlaps %>%
      union(nextjoin %>% filter(to < next_from) %>% select(accno, profile, from, to))
    nextjoin <- nextjoin %>%
      anti_join(
        nextjoin %>% filter(to < next_from) %>% select(accno, profile, from, to),
        by = c('accno', 'profile', 'from', 'to')
      )
    
    # 4. Delete rows which are contained in the earlier row
    nextjoin <- nextjoin %>%
      anti_join(
        nextjoin %>% filter(from < next_from, to > next_to) %>%
          transmute(accno, profile, from = next_from, to = next_to),
        by = c('accno', 'profile', 'from', 'to')
      )
    
    # 5. Set a new "to" for those that overlap with the next row
    nextjoin <- nextjoin %>%
      mutate(to = ifelse(! is.na(next_from) & to >= next_from & next_to > to, next_to, to))
    
    # 6. Delete rows that are the last in their group of overlaps, they now have
    #   the same "to" as the previous row.
    nextjoin <- nextjoin %>%
      anti_join(
        nextjoin %>% select(accno, profile, from, to) %>% 
          inner_join(nextjoin %>% select(accno, profile, from, to), by = c('accno', 'profile', 'to')) %>% 
          filter(from.x != from.y) %>% 
          group_by(accno, profile, to) %>% summarise(from = max(from.x)) %>% ungroup(),
        by = c('accno', 'profile', 'from', 'to')
      )

    # 5. Calculate a new domt from nextjoin
    domt <- nextjoin %>% distinct(accno, profile, from, to) %>%
      group_by(accno, profile) %>% mutate(i = rank(from), n = n()) %>% ungroup()
  }
  
  lengths <- lengths %>%
    union(
      nooverlaps %>% mutate(len = to - from + 1) %>%
        group_by(accno, profile) %>% summarise(len = sum(len), from = min(from), to = max(to)) %>% ungroup() %>%
        gather(type, val, len, from, to) %>%
        mutate(
          type = case_when(
            type == 'len' ~ fs[3],
            type == 'from' ~ fs[1],
            type == 'to'   ~ fs[2]
          )
        )
    )
}

# Join in the above results with tlen and qlen from domtblout
align_lengths = domtblout %>% distinct(accno, profile, tlen, qlen) %>%
  inner_join(lengths %>% spread(type, val, fill = 0), by = c('accno', 'profile'))

logmsg("Calculated lengths, inferring source databases from accession numbers")

# Infer databases from the structure of accession numbers
accessions = accessions %>%
  mutate(db = ifelse(grepl('^.._', accto), 'refseq', NA)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[0-9A-Z]{4,4}_[0-9A-Z]$', accto)), 'pdb', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^P[0-9]+\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[O,P,Q][0-9][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[ADEKOJMNP][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'genbank', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[C][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'embl', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[BFGIL][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'dbj', db))

logmsg("Inferred databases, calculating best scoring profile for each accession")

# Create proteins with entries from tblout not matching hmm_profile entries with rank == 'domain'.
# Calculate best scoring profile for each accession
proteins <- tblout %>% 
  anti_join(hmm_profiles %>% filter(prank == 'domain'), by = 'profile') %>%
  group_by(accno) %>% top_n(1, score) %>% ungroup() %>%
  select(accno, profile, score, evalue)

logmsg("Calculated best scoring profiles, creating domains")

# Create table of domains as those that match domains specified in hmm_profiles
domains <- domtblout %>%
  semi_join(hmm_profiles %>% filter(prank == 'domain'), by = 'profile') %>%
  transmute(
   accno, profile, i, n,
   dom_c_evalue, dom_i_evalue, dom_score,
   hmm_from, hmm_to,
   ali_from, ali_to,
   env_from, env_to
)

# Join in lengths
logmsg(sprintf("Joining in lengths from domtblout, nrows before: %d", proteins %>% nrow()))
proteins <- proteins %>% inner_join(align_lengths, by = c('accno', 'profile'))

logmsg("Joined in lengths, writing data")

# If we were called with the singletable option, prepare data suitable for that
if ( length(grep('singletable', names(opt$options), value = TRUE)) > 0 ) {
  logmsg("Writing single table format")

  # Join proteins with accessions and drop profile to get a single table output
  logmsg(sprintf("Joining in all accession numbers and dropping profile column, nrows before: %d", proteins %>% nrow()))
  singletable <- proteins %>% 
    left_join(hmm_profiles, by='profile') %>%
    inner_join(accessions, by='accno') %>%
    mutate(accno = accto) 

  # If we have a taxflat NCBI taxonomy, read and join
  logmsg(sprintf("Adding NCBI taxon ids from taxflat, nrows before: %d", singletable %>% nrow()))
  singletable = singletable %>% 
    left_join(
      taxflat %>% select(taxon, ncbi_taxon_id),
      by='taxon'
    )

  logmsg(sprintf("Writing single table %s, nrows: %d", opt$options$singletable, singletable %>% nrow()))
  write_tsv(
    singletable %>% 
      select(db, accno, taxon, score, evalue, profile, psuperfamily:pgroup, ncbi_taxon_id, tlen, qlen, alilen, hmmlen,envlen) %>%
      arrange(accno, profile),
    opt$options$singletable
  )
}

# If the user specified a filename for a SQLite database, write that here
if ( length(grep('sqlitedb', names(opt$options), value = TRUE)) > 0 ) {
  logmsg(sprintf("Creating/opening SQLite database %s", opt$options$sqlitedb))
  con = DBI::dbConnect(RSQLite::SQLite(), opt$options$sqlitedb, create = TRUE)

  con %>% copy_to(
    tibble(source = dbsource[1], name = dbsource[2], version = dbsource[3]), 
    'dbsources', temporary = FALSE, overwrite = TRUE
  )

  # The accto field in accession should be turned into a list for each
  # combination of accno, db and taxon to ensure organisms do not show up as
  # having more than one exactly identical sequence, which they do with the new
  # redundant RefSeq entries (WP_ accessions).
  logmsg('Copying to "accessions", creating indices')
  accessions <- accessions %>%
    arrange(db, taxon, accno, accto) %>%
    group_by(db, taxon, accno) %>%
    summarise(accto = paste(accto, collapse = ',')) %>%
    ungroup()
  
  con %>% copy_to(accessions, 'accessions', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE UNIQUE INDEX "accessions.i00" ON "accessions"("db", "accno", "taxon");')
  con %>% DBI::dbExecute('CREATE INDEX "accessions.i01" ON "accessions"("accno");')
  con %>% DBI::dbExecute('CREATE INDEX "accessions.i02" ON "accessions"("taxon");')
  
  logmsg('Copying to "proteins", creating indices')
  con %>% copy_to(proteins, 'proteins', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE INDEX "proteins.i00" ON "proteins"("accno");')
  con %>% DBI::dbExecute('CREATE INDEX "proteins.i01" ON "proteins"("profile");')

  logmsg('Copying to "domains", creating indices')
  con %>% copy_to(domains, 'domains', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE INDEX "domains.i00" ON "domains"("accno");')
  con %>% DBI::dbExecute('CREATE INDEX "domains.i01" ON "domains"("profile");')

  logmsg('Copying to "hmm_profiles", creating indices')
  con %>% copy_to(hmm_profiles, 'hmm_profiles', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE UNIQUE INDEX "hmm_profiles.i00" ON "hmm_profiles"("profile");')

  # If we have a taxflat NCBI taxonomy, read and join
  logmsg(sprintf("Adding NCBI taxon ids from taxflat"))
  con %>% copy_to(
    taxflat %>% inner_join(accessions %>% distinct(taxon), by='taxon'),
    'taxa', temporary = FALSE, overwrite = TRUE
  )

  logmsg('Creating indices on "taxa"')
  con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i00" ON "taxa"("taxon", "trank");')
  con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i01" ON "taxa"("ncbi_taxon_id");')

  logmsg('Saving tblout and domtblout to database')
  con %>% copy_to(tblout    %>% arrange(accno, profile),    'tblout',    temporary = FALSE, overwrite = TRUE)
  con %>% copy_to(domtblout %>% arrange(accno, profile, i), 'domtblout', temporary = FALSE, overwrite = TRUE)
}

logmsg("Done")
