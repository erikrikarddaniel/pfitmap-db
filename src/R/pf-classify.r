#!/usr/bin/env Rscript

# pf-classify
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(feather))

SCRIPT_VERSION = "1.9.6"

options(warn = 1)

# Get arguments
option_list <- list(
  make_option(
    c('--dbsource'), default='', help='Database source in dbsource:name:version format'
  ),
  make_option(
    c('--featherprefix'), default = '', help = 'Turn on generation of feather output files, and use this sting to prefix file names with.',
  ),
  make_option(
    c('--fuzzy_factor'), type = 'integer', default = 1, action = 'store', help = 'Factor to make lengths fuzzy for reduction of possible duplicates, default %default.'
  ),
  make_option(
    c('--gtdbmetadata'), default = '', help = 'A concatenation of the bacterial (bac120_metadata.tsv) and archaeal (ar122_metadata.tsv) metadata files from GTDB. Make sure there is *only one header* line.',
  ),
  make_option(
    c('--hmm_mincov'), type = 'double', default = 0.0, action = 'store', help = 'Minimum coverage of hmm profile to include in output, default %default.'
  ),
  make_option(
  c('--profilehierarchies'), default='', help='A tsv file with profile hiearchies, including a header. Required fields: profile and plen, recommended psuperfamily, pfamily pclass, psubclass, pgroup, prank and version.'
  ),
  make_option(
    c('--seqfaa'), default='', help='Fasta file with amino acid sequences, same as the one used as database when hmmsearching. Will populate a sequence table if --sqlitedb is set.',
  ),
  make_option(
    c('--singletable'), default='', help='Write data in a single tsv format to this filename.'
  ),
  make_option(
    c('--sqlitedb'), default='', help='Write data in a SQLite database with this filename.'
  ),
  make_option(
    c('--taxflat'), default='', help='Name of NCBI taxon table in "taxflat" format (see https://github.com/erikrikarddaniel/taxdata2taxflat).'
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
opt <- parse_args(
  OptionParser(
    usage = "%prog [options] file0.tblout .. filen.tblout file0.domtblout .. filen.domtblout", 
    option_list = option_list
  ), 
  positional_arguments = TRUE
)

if ( opt$options$version ) {
  write(SCRIPT_VERSION[1], stdout())
  quit('no')
}

if ( length(grep('sqlitedb', names(opt$options), value = TRUE)) > 0 ) {
  suppressPackageStartupMessages(library(dbplyr))
}

# Args list for testing:
# NCBI: opt = list(args = Sys.glob('pf-classify.00.d/*tblout'), options=list(verbose=T, singletable='test.out.tsv', hmm_mincov=0.9, profilehierarchies='pf-classify.00.phier.tsv', taxflat='pf-classify.taxflat.tsv', sqlitedb='testdb.sqlite3', dbsource='NCBI:NR:20180212', fuzzy_factor=30))
# GTDB: opt = list(args = Sys.glob('pf-classify.gtdb.02.d/*tblout'), options=list(verbose=T, singletable='test.out.tsv', hmm_mincov=0.9, profilehierarchies='pf-classify.gtdb.02.phier.tsv', taxflat='pf-classify.taxflat.tsv', sqlitedb='testdb.sqlite3', dbsource='GTDB:GTDB:r86', fuzzy_factor=30, gtdbannotindex='pf-classify.gtdb.02.d/gtdb_prokka_index.tsv.gz', gtdbmetadata='pf-classify.gtdb.02.d/gtdb_metadata.tsv', gtdbtaxonomy='pf-classify.gtdb.02.d/gtdb_taxonomy.tsv', seqfaa='pf-classify.gtdb.03.d/genomes.faa'))
DEBUG   = 0
INFO    = 1
WARNING = 2
LOG_LEVELS = list(
  DEBUG   = list(n = 0, msg = 'DEBUG'),
  INFO    = list(n = 1, msg = 'INFO'),
  WARNING = list(n = 2, msg = 'WARNING'),
  ERROR   = list(n = 3, msg = 'ERROR')
)
logmsg    = function(msg, llevel='INFO') {
  if ( opt$options$verbose | LOG_LEVELS[[llevel]][["n"]] >= LOG_LEVELS[["WARNING"]][["n"]] ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg(sprintf("pf-classify.r version %s: Starting classification", SCRIPT_VERSION))
if ( length(grep('featherprefix', names(opt$options), value = TRUE)) > 0 & str_length(opt$options$featherprefix) > 0 ) logmsg(sprintf("Will write feather files prefixed with %s", opt$options$featherprefix))

# Make sure the dbsource parameter is given and in the proper format
if ( opt$options[['dbsource']] != '' ) {
  dbsource <- strsplit(opt$options$dbsource, ':')[[1]]
} else {
  logmsg(sprintf("--dbsource is required"), 'ERROR')
  quit('no', status = 2)
}

if ( length(dbsource) != 3 ) {
  logmsg(sprintf("--dbsource needs to contain three components separated with ':', you specified '%s'", opt$options$dbsource), 'ERROR')
  quit('no', status = 2)
}

# Check if the GTDB metadata parameter is given
gtdb <- ifelse(opt$options$gtdbmetadata > '', TRUE, FALSE)

logmsg(sprintf("Reading profile hierarchies from %s", opt$options$profilehierarchies))
hmm_profiles <- fread(opt$options$profilehierarchies)

# Check that it's unique on profile
if ( nrow(hmm_profiles) != nrow(unique(hmm_profiles[, .(profile)])) ) {
  logmsg(sprintf("The profile column in the hmm_profiles table (%s) needs to be unique", opt$options$profilehierarchies), 'ERROR')
  quit('no', status = 2)
}

# Read the taxonomy file, in GTDB or NCBI format
if ( gtdb ) {
  logmsg(sprintf('Reading GTDB metadata from %s', opt$options$gtdbmetadata))
  # Don't know how to separate() with data.table, read as tibble then convert
  gtdbmetadata <- read_tsv(
    opt$options$gtdbmetadata,
    col_types = cols(.default = col_character())
  ) %>%
    mutate(
      thier = str_remove_all(gtdb_taxonomy, '[a-z]__'), 
    ) %>%
    separate(thier, c('tdomain', 'tphylum', 'tclass', 'torder', 'tfamily', 'tgenus', 'tspecies'), sep = ';') %>%
    as.data.table()
  gtdbtaxonomy <- lazy_dt(gtdbmetadata) %>%
    mutate(
      accno0 = str_remove(accession, '^RS_') %>% str_remove('^GB_') %>% str_remove('\\.[0-9]'), 
      accno1 = ncbi_genbank_assembly_accession %>% str_remove('\\.[0-9]'),
      trank  = 'species',
      ncbi_taxon_id = ncbi_species_taxid
    ) %>%
    select(accno0, accno1, tdomain:tspecies, trank, ncbi_taxon_id) %>%
    as.data.table()
} else {
  logmsg(sprintf("Reading NCBI taxonomy from %s", opt$options$taxflat))
  taxflat <- fread(opt$options$taxflat) %>%
    lazy_dt() %>%
    transmute(
      ncbi_taxon_id, taxon, trank = rank,
      tdomain       = superkingdom, tkingdom = kingdom,
      tphylum       = phylum,       tclass   = class,
      torder        = order,        tfamily  = family,
      tgenus        = genus,        tspecies = species
    ) %>%
    as.data.table()

  # Delete duplicate taxon, rank combinations belonging in Eukaryota
  taxflat <- lazy_dt(taxflat) %>%
    anti_join(
      lazy_dt(taxflat) %>% group_by(taxon, trank) %>% summarise(n = n(), .groups = 'drop_last') %>% ungroup() %>% filter(n > 1) %>%
        inner_join(lazy_dt(taxflat) %>% filter(tdomain == 'Eukaryota'), by = c('taxon', 'trank')),
      by = c('ncbi_taxon_id')
    ) %>%
    as.data.table()
}

# We will populate two tables, one with the full results, one with accessions
tblout <- data.table(
  accno = character(), profile = character(),
  evalue = double(), score = double(), bias = double()
)
accessions <- data.table(accno = character(), accto = character())

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
    ) %>%
    as.data.table()
  tblout <- funion(tblout, t[, .(accno, profile, evalue, score, bias)])
  accessions <- funion(
    accessions, 
    t[, {  accno <- accno; accto <- sprintf("%s %s", accno, rest); .(accno = accno, accto = accto) }]
  )
}

# Split the accto field (data.table magic I don't understand, see https://stackoverflow.com/questions/13773770/split-comma-separated-strings-in-a-column-into-separate-rows/31514711#31514711)
accessions <- accessions[, strsplit(accto, "\x01", fixed=TRUE), by = .(accno, accto)][,.(accno, accto = V1)] %>% 
  lazy_dt() %>%
  mutate(
    taxon = ifelse(
      grepl('[^[]\\[(.*)\\]', accto), 
      sub('.*[^[]\\[(.*)\\].*', '\\1', accto),
      'unknown'
    ),
    accto = sub(' .*', '', accto)
  ) %>%
  as.data.table()

domtblout <- data.table(
  accno = character(), tlen = integer(), profile = character(), qlen = integer(), i = integer(), n = integer(), 
  dom_c_evalue = double(), dom_i_evalue = double(), dom_score = double(), dom_bias = double(),
  hmm_from = integer(), hmm_to = integer(), ali_from = integer(), ali_to = integer(), 
  env_from = integer(), env_to = integer()
)

# Read all the domtblout files
for ( domtbloutfile in grep('\\.domtblout', opt$args, value=TRUE) ) {
  logmsg(sprintf("Reading %s", domtbloutfile))
  t <- read_fwf(
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
    ) %>%
    as.data.table()
  
  domtblout <- funion(
    domtblout,
    t[, .(
      accno, tlen, profile, qlen, i, n, dom_c_evalue, dom_i_evalue, dom_score, dom_bias,
      hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
    )]
  )
}

# Calculate lengths:
# This is long-winded because we have three lengths (hmm, ali and env) which
# have to be calculated separately after getting rid of overlapping rows in
# domtblout.

# Define a temporary table that will be filled with lengths, minimum from and maximum
# to values.
lengths <- data.table(
  accno = character(), profile = character(), type = character(), val = integer()
)

# Function that joins all with n > 1 with the next row and not occurring in the second
# table in the call.
do_nextjoin <- function(dt) {
  # Filter dt on n > 1, then join with itself with transmuted column names
  suppressWarnings( # To avoid warnings about, I believe, empty data.tables
    t <- lazy_dt(dt) %>% 
      filter(n > 1) %>%
      left_join(
        lazy_dt(dt) %>% transmute(accno, profile, i = i - 1, next_row = TRUE, next_from = from, next_to = to), 
        by = c('accno', 'profile', 'i')
      ) %>%
      as.data.table()
  )
  # Set any NA next_row entries to FALSE (doesn't work through dtplyr)
  t$next_row = ifelse(is.na(t$next_row), FALSE, TRUE)
  return(t)
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
  logmsg("Creating domt table")
  domt <- domtblout[, .(accno = accno, profile = profile, from = get(fs[1]), to = get(fs[2]))]
  logmsg("Joining with itself, calculating row numbers and more")
  domt <- lazy_dt(domt) %>% distinct() %>%
    semi_join(
      lazy_dt(domt) %>% group_by(accno, profile, from) %>% filter(to == max(to)) %>% ungroup(),
      by = c('accno', 'profile', 'from', 'to')
    ) %>% 
    arrange(accno, profile, from, to) %>%
    group_by(accno, profile) %>% mutate(i = row_number(from), n = n()) %>% ungroup() %>%
    as.data.table()
  logmsg("domt done")

  # This is where non-overlapping rows will be stored
  nooverlaps <- data.table(
    accno = character(), profile = character(),
    from = integer(), to = integer()
  )

  while ( nrow(domt) > 0 ) {
    logmsg(sprintf("Working on overlaps, nrow: %d", domt %>% nrow()))
    
    # 2. Move rows with n == 1 to nooverlaps
    nooverlaps <- funion(nooverlaps, domt[n == 1, .(accno, profile, from, to)])
    domt <- domt[n > 1]
    
    nextjoin <- do_nextjoin(domt)

    # As a debug aid, save the domt and nextjoin data sets
    #write_tsv(domt, 'domt.tsv.gz')
    #write_tsv(nextjoin, 'nextjoin.tsv.gz')
  
    # 3. Move rows that do not overlap with the next row
    nooverlaps <- funion(nooverlaps, nextjoin[to < next_from, .(accno, profile, from, to)])
    nextjoin <- lazy_dt(nextjoin) %>%
      anti_join(
        lazy_dt(nextjoin) %>% filter(to < next_from) %>% select(accno, profile, from, to),
        by = c('accno', 'profile', 'from', 'to')
      ) %>%
      as.data.table()
    
    # 4. Delete rows which are contained in the earlier row
    nextjoin <- lazy_dt(nextjoin) %>%
      anti_join(
        lazy_dt(nextjoin) %>% filter(from < next_from, to > next_to) %>%
          transmute(accno, profile, from = next_from, to = next_to),
        by = c('accno', 'profile', 'from', 'to')
      ) %>%
      as.data.table()
    
    # 5. Set a new "to" for those that overlap with the next row
    nextjoin <- lazy_dt(nextjoin) %>%
      mutate(to = ifelse(! is.na(next_from) & to >= next_from & next_to > to, next_to, to)) %>%
      as.data.table()
    
    # 6. Delete rows that are the last in their group of overlaps, they now have
    #   the same "to" as the previous row.
    suppressWarnings( # The max(from.x) causes warnings: "no non-missing arguments to max; returning -Inf"; I can't find any errors
      nextjoin <- lazy_dt(nextjoin) %>%
        anti_join(
          lazy_dt(nextjoin) %>% select(accno, profile, from, to) %>% 
            inner_join(lazy_dt(nextjoin) %>% select(accno, profile, from, to), by = c('accno', 'profile', 'to')) %>% 
            filter(from.x != from.y) %>% 
            group_by(accno, profile, to) %>% summarise(from = max(from.x), .groups = 'drop_last') %>% ungroup(),
          by = c('accno', 'profile', 'from', 'to')
        ) %>%
        as.data.table()
    )

    # 7. Calculate a new domt from nextjoin
    if ( nrow(nextjoin) > 0 ) {
      domt <- lazy_dt(nextjoin) %>% distinct(accno, profile, from, to) %>%
        group_by(accno, profile) %>% mutate(i = rank(from), n = n()) %>% ungroup() %>%
        as.data.table()
    }
  }
  
  # I'm resorting to tibbles here because of the gather() buried relatively deep in the pipeline
  lengths <- as_tibble(lengths) %>%
    union(
      as_tibble(nooverlaps) %>% mutate(len = to - from + 1) %>%
        group_by(accno, profile) %>% summarise(len = sum(len), from = min(from), to = max(to), .groups = 'drop_last') %>% ungroup() %>%
        gather(type, val, len, from, to) %>% # Change to pivot_longer()!
        mutate(
          type = case_when(
            type == 'len' ~ fs[3],
            type == 'from' ~ fs[1],
            type == 'to'   ~ fs[2]
          )
        )
    ) %>%
    as.data.table()
}

# Join in the above results with tlen and qlen from domtblout
align_lengths <- lazy_dt(domtblout) %>% 
  distinct(accno, profile, tlen, qlen) %>%
  inner_join(as_tibble(lengths) %>% spread(type, val, fill = 0), by = c('accno', 'profile')) %>% as.data.table() %>% lazy_dt() %>%
  as.data.table()

logmsg("Calculated lengths, inferring source databases from accession numbers")

if ( gtdb ) {
  accessions$db <- 'gtdb'
  accessions <- lazy_dt(accessions) %>%
    transmute(db, genome_accno = str_remove(taxon, '\\..*'), accno) %>%
    as.data.table()
} else {
  # Infer databases from the structure of accession numbers
  accessions <- lazy_dt(accessions) %>%
    mutate(db = ifelse(grepl('^.._', accto), 'refseq', NA)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[0-9A-Z]{4,4}_[0-9A-Z]$', accto)), 'pdb', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^P[0-9]+\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[O,P,Q][0-9][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[ADEKOJMNP][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'genbank', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[C][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'embl', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[BFGIL][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'dbj', db)) %>%
    as.data.table()
}

logmsg("Inferred databases, calculating best scoring profile for each accession")

# Create proteins with entries from tblout not matching hmm_profile entries with rank == 'domain'.
# Calculate best scoring profile for each accession
proteins <- lazy_dt(tblout) %>% 
  anti_join(lazy_dt(hmm_profiles) %>% filter(prank == 'domain'), by = 'profile') %>%
  group_by(accno) %>% filter(score == max(score)) %>% ungroup() %>%
  select(accno, profile, score, evalue) %>%
  as.data.table()
proteins <- setorder(setDT(proteins), -score)[, head(.SD, 1), keyby = accno]

logmsg("Calculated best scoring profiles, creating domains")

# Create table of domains as those that match domains specified in hmm_profiles
domains <- lazy_dt(domtblout) %>%
  transmute(
    accno, profile, i, n, dom_c_evalue, dom_i_evalue, dom_score,
    hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
  ) %>%
  as.data.table()

# Join in lengths
logmsg(sprintf("Joining in lengths from domtblout, nrows before: %d", proteins %>% nrow()))
proteins <- lazy_dt(proteins) %>% inner_join(align_lengths, by = c('accno', 'profile')) %>% as.data.table()

logmsg("Joined in lengths, writing data")

# Subset data with hmm_mincov parameter
logmsg(sprintf("Subsetting output to proteins and domains covering at least %f of the hmm profile", opt$options$hmm_mincov))

# 1. proteins table
#logmsg(sprintf("proteins before: %d", proteins %>% nrow()), 'DEBUG')
p <- lazy_dt(proteins) %>% 
  inner_join(lazy_dt(hmm_profiles) %>% select(profile, plen), by = 'profile') %>%
  filter(hmmlen/plen >= opt$options$hmm_mincov) %>%
  select(-plen) %>%
  as.data.table()
#logmsg(sprintf("proteins after: %d", proteins %>% nrow()), 'DEBUG')
p <- lazy_dt(p) %>%
  union(
    lazy_dt(proteins) %>% anti_join(lazy_dt(p), by = 'accno') %>%
      semi_join(lazy_dt(accessions) %>% filter(db == 'pdb'), by = c('accno' = 'accno'))
  ) %>%
  as.data.table()
proteins <- p
rm(p)

# 1.b add pdb entries that are not present due to not passing the hmm_mincov criterion

# 2. domains
domains <- lazy_dt(domains) %>%
  inner_join(lazy_dt(hmm_profiles) %>% select(profile, plen), by = 'profile') %>%
  filter((hmm_to - hmm_from + 1)/plen >= opt$options$hmm_mincov) %>%
  select(-plen) %>%
  as.data.table()

# 3. accessions
accessions <- lazy_dt(accessions) %>% 
  semi_join(
    union(lazy_dt(proteins) %>% select(accno), lazy_dt(domains) %>% select(accno)) %>% distinct(accno), 
    by = 'accno'
  ) %>%
  as.data.table()

# 4. tblout
tblout <- lazy_dt(tblout) %>% semi_join(lazy_dt(accessions) %>% distinct(accno), by = 'accno') %>% as.data.table()

# 5. domtblout
domtblout <- lazy_dt(domtblout) %>% semi_join(lazy_dt(accessions) %>% distinct(accno), by = 'accno') %>% as.data.table()

# Sequences in fasta file?
if ( opt$options$seqfaa != '' ) {
  logmsg(sprintf('Reading %s', opt$options$seqfaa))
  cmd <- c(
    sprintf("%s %s", ifelse(grepl('\\.gz$', opt$options$seqfaa), 'zcat', 'cat'), opt$options$seqfaa),
    "sed '/^>/s/ .*//'",
    "awk '/^>/ { printf(\"\\n%s\\n\",$$0); next; } { printf(\"%s\",$$0);} END { printf(\"\\n\");}'",
    "grep -v '^$'", 
    "sed 's/^>//'", 
    "paste - -"
  )
  sequences <- fread(cmd = paste(cmd, collapse = '|'), col.names = c('accno', 'sequence'), header = FALSE) %>% lazy_dt() %>%
    semi_join(lazy_dt(accessions), by = 'accno') %>% 
    as.data.table()
  logmsg(sprintf("Read %d sequences", nrow(sequences)))
}

# If we were called with the singletable option, prepare data suitable for that
if ( opt$options$singletable > '' ) {
  logmsg("Writing single table format")

  # Join proteins with accessions and drop profile to get a single table output
  logmsg(sprintf("Joining in all accession numbers, nrows before: %d", proteins %>% nrow()))
  singletable <- as_tibble(proteins) %>% 
    left_join(as_tibble(hmm_profiles), by='profile') %>% # Why is this a *left* join?
    inner_join(as_tibble(accessions), by='accno')
  
  if ( ! gtdb ) singletable <- singletable %>% mutate(accno = accto) 

  # Join in taxonomies, either GTDB or taxflat
  if ( gtdb ) {
    singletable <- union(
      singletable %>% 
        inner_join(as_tibble(gtdbtaxonomy), by = c('genome_accno' = 'accno0')) %>% select(-accno1),
      singletable %>% 
        anti_join(as_tibble(gtdbtaxonomy),  by = c('genome_accno' = 'accno0')) %>% 
        left_join(as_tibble(gtdbtaxonomy), by = c('genome_accno' = 'accno1')) %>% select(-accno0)
    )
    if ( singletable %>% filter(is.na(tspecies)) %>% nrow() > 0 ) {
      logmsg(
        sprintf("*** Accessions without GTDB species assignment: %s ***", singletable %>% filter(is.na(tspecies)) %>% pull(genome_accno) %>% paste(collapse = ', ')),
        "WARNING"
      )
    }
    logmsg(sprintf("Writing single table %s, nrows: %d", opt$options$singletable, singletable %>% nrow()))
    write_tsv(
      singletable %>% 
        select(db, accno, score, evalue, profile, psuperfamily:pgroup, genome_accno, tdomain:tspecies, tlen, qlen, alilen, hmmlen,envlen) %>%
        arrange(accno, profile),
      opt$options$singletable
      )
  } else {
    logmsg(sprintf("Adding NCBI taxon ids from taxflat, nrows before: %d", singletable %>% nrow()))
    singletable <- singletable %>% 
      left_join(
        as_tibble(taxflat) %>% select(taxon, ncbi_taxon_id),
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
}

if ( ! gtdb ) {
  # The accto field in accession should be turned into a list for each
  # combination of accno, db and taxon to ensure organisms do not show up as
  # having more than one exactly identical sequence, which they do with the new
  # redundant RefSeq entries (WP_ accessions).
  logmsg('Copying to "accessions", creating indices')
  # I'm converting to tibble here, as I don't know how to to the separate_rows in
  # data.table and it's not in dtplyr.
  accessions <- as_tibble(accessions) %>%
    arrange(db, taxon, accno, accto) %>%
    group_by(db, taxon, accno) %>%
    summarise(accto = paste(accto, collapse = ','), .groups = 'drop_last') %>%
    ungroup() %>%
    separate_rows(accto, sep = ',') %>%
    distinct()
}
if ( gtdb ) {
  taxa <- union(
    as_tibble(gtdbtaxonomy) %>% semi_join(as_tibble(accessions), by = c('accno0' = 'genome_accno')) %>% 
      select(-accno1) %>% rename(genome_accno = accno0),
    as_tibble(gtdbtaxonomy) %>% anti_join(as_tibble(accessions), by = c('accno0' = 'genome_accno')) %>% 
      mutate(genome_accno = ifelse(accno1 == 'none', accno0, accno1)) %>%
      select(-accno0, -accno1)
  )
} else {
  taxa <- lazy_dt(taxflat) %>% semi_join(lazy_dt(accessions) %>% distinct(taxon), by='taxon') %>% as.data.table()
}

# If the user specified a filename for a SQLite database, write that here
if ( length(grep('sqlitedb', names(opt$options), value = TRUE)) > 0 & str_length(opt$options$sqlitedb) > 0 ) {
  logmsg(sprintf("Creating/opening SQLite database %s", opt$options$sqlitedb))
  con <- DBI::dbConnect(RSQLite::SQLite(), opt$options$sqlitedb, create = TRUE)

  con %>% copy_to(
    tibble(source = dbsource[1], name = dbsource[2], version = dbsource[3]), 
    'dbsources', temporary = FALSE, overwrite = TRUE
  )
  
  con %>% copy_to(accessions, 'accessions', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE INDEX "accessions.i00" ON "accessions"("accno");')
  if ( gtdb ) {
    con %>% DBI::dbExecute('CREATE INDEX "accessions.i01" ON "accessions"("db", "accno", "genome_accno");')
    con %>% DBI::dbExecute('CREATE INDEX "accessions.i02" ON "accessions"("genome_accno");')
  } else {
    con %>% DBI::dbExecute('CREATE INDEX "accessions.i01" ON "accessions"("db", "accto", "taxon");')
    con %>% DBI::dbExecute('CREATE INDEX "accessions.i02" ON "accessions"("taxon");')
  }
  
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

  logmsg('Copying to "taxa", creating indices')
  con %>% copy_to(taxa, 'taxa', temporary = FALSE, overwrite = TRUE)
  if ( gtdb ) {
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i00" ON "taxa"("genome_accno");')
  } else {
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i00" ON "taxa"("taxon", "trank");')
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i01" ON "taxa"("ncbi_taxon_id");')
  }

  logmsg('Saving tblout and domtblout to database')
  con %>% copy_to(tblout    %>% arrange(accno, profile),    'tblout',    temporary = FALSE, overwrite = TRUE)
  con %>% copy_to(domtblout %>% arrange(accno, profile, i), 'domtblout', temporary = FALSE, overwrite = TRUE)

  if ( ! gtdb ) {
  logmsg(sprintf('Creating dupfree_proteins, using %d as fuzzy factor', opt$options$fuzzy_factor))
    dp <- as_tibble(proteins) %>% inner_join(as_tibble(accessions) %>% transmute(accno = accto, db, taxon), by = 'accno') %>%
      mutate(
        alilen = as.integer(round(round(alilen / opt$options$fuzzy_factor) * opt$options$fuzzy_factor)),
        envlen = as.integer(round(round(envlen / opt$options$fuzzy_factor) * opt$options$fuzzy_factor)),
        hmmlen = as.integer(round(round(hmmlen / opt$options$fuzzy_factor) * opt$options$fuzzy_factor))
      ) %>%
      group_by(db, taxon, profile, alilen, hmmlen, envlen) %>% mutate(r = rank(accno)) %>% ungroup()

    con %>% copy_to(
      dp %>% filter(r < 2) %>% select(-db, -taxon),
      'dupfree_proteins',
      temporary = FALSE, overwrite = TRUE
    )
  }

  if ( opt$options$seqfaa != '' ) {
    logmsg(sprintf('Saving sequences table', opt$options$seqfaa))
    con %>% copy_to(sequences, 'sequences',
      temporary = FALSE, overwrite = TRUE
    )
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "sequences.i00" ON "sequences"("accno");')
  }

  logmsg('Disconnecting from sqlite3 db')
  con %>% DBI::dbDisconnect()
}

if ( length(grep('featherprefix', names(opt$options), value = TRUE)) > 0 & str_length(opt$options$featherprefix) > 0 ) {
  logmsg("Writing dbsources to feather file")
  write_feather(
    tibble(source = dbsource[1], name = dbsource[2], version = dbsource[3]),
    sprintf("%s.dbsources.feather", opt$options$featherprefix)
  )

  logmsg("Writing accessions to feather file")
  write_feather(accessions, sprintf("%s.accessions.feather", opt$options$featherprefix))

  logmsg("Writing proteins to feather file")
  write_feather(proteins, sprintf("%s.proteins.feather", opt$options$featherprefix))

  logmsg("Writing domains to feather file")
  write_feather(domains, sprintf("%s.domains.feather", opt$options$featherprefix))

  logmsg("Writing hmm_profiles to feather file")
  write_feather(hmm_profiles, sprintf("%s.hmm_profiles.feather", opt$options$featherprefix))

  logmsg("Writing taxa to feather file")
  write_feather(taxa, sprintf("%s.taxa.feather", opt$options$featherprefix))

  logmsg("Writing tblout to feather file")
  write_feather(tblout, sprintf("%s.tblout.feather", opt$options$featherprefix))

  logmsg("Writing domtblout to feather file")
  write_feather(domtblout, sprintf("%s.domtblout.feather", opt$options$featherprefix))

  if ( opt$options$seqfaa != '' ) {
    logmsg("Writing sequences to feather file")
    write_feather(sequences, sprintf("%s.sequences.feather", opt$options$featherprefix))
  }
}

logmsg("Done")
