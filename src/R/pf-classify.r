#!/usr/bin/env Rscript

# pf-classify
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(feather))

SCRIPT_VERSION = "1.9.14"
ROWS_PER_SEQUENCE_TSV = 1e7

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
    c('--missing'), default='', help='Name of file to write missing genome information to, default not set. Only works in GTDB mode!'
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
# GTDB: opt = list(args = Sys.glob('pf-classify.gtdb.02.d/*tblout'), options=list(verbose=T, singletable='test.out.tsv', hmm_mincov=0.9, profilehierarchies='pf-classify.gtdb.02.phier.tsv', taxflat='pf-classify.taxflat.tsv', sqlitedb='testdb.sqlite3', dbsource='GTDB:GTDB:r86', fuzzy_factor=30, gtdbmetadata='pf-classify.gtdb.02.d/gtdb_metadata.tsv', seqfaa='pf-classify.gtdb.03.d/genomes.faa'))
# GTDB 06: opt = list(args = Sys.glob('pf-classify.gtdb.06.d/*tblout'), options=list(verbose=T, featherprefix='/tmp/pf-classify-testing', hmm_mincov=0.8, profilehierarchies='pf-classify.gtdb.06.phier.tsv', dbsource='GTDB:GTDB:r86', fuzzy_factor=30, gtdbmetadata='pf-classify.gtdb.06.d/gtdb_metadata.tsv', seqfaa='pf-classify.gtdb.06.d/genomes.faa'))
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
  if ( opt$options$verbose | LOG_LEVELS[[llevel]][["n"]] >= LOG_LEVELS[["INFO"]][["n"]] ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg(sprintf("pf-classify.r version %s: Starting classification", SCRIPT_VERSION))

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

logmsg(sprintf("Reading profile hierarchies from %s", opt$options$profilehierarchies), 'DEBUG')
hmm_profiles <- data.table::fread(opt$options$profilehierarchies)

# Check that it's unique on profile
if ( nrow(hmm_profiles) != nrow(unique(hmm_profiles[, .(profile)])) ) {
  logmsg(sprintf("The profile column in the hmm_profiles table (%s) needs to be unique", opt$options$profilehierarchies), 'ERROR')
  quit('no', status = 2)
}

# Read the taxonomy file, in GTDB or NCBI format
if ( gtdb ) {
  logmsg(sprintf('Reading GTDB metadata from %s', opt$options$gtdbmetadata), 'DEBUG')
  gtdbmetadata <- read_tsv(
    opt$options$gtdbmetadata,
    col_types = cols(.default = col_character())
  ) %>%
    mutate(
      thier = str_remove_all(gtdb_taxonomy, '[a-z]__'), 
    ) %>%
    separate(thier, c('tdomain', 'tphylum', 'tclass', 'torder', 'tfamily', 'tgenus', 'tspecies'), sep = ';')
  #print(gtdbmetadata %>% transmute(accno0 = str_remove(accession, '^RS_') %>% str_remove('^GB_') %>% str_remove('\\.[0-9]'), accession)
  gtdbtaxonomy <- gtdbmetadata %>%
    mutate(
      accno0 = str_remove(accession, '^RS_') %>% str_remove('^GB_') %>% str_remove('\\.[0-9]'), 
      accno1 = ncbi_genbank_assembly_accession %>% str_remove('\\.[0-9]'),
      trank  = 'species',
      ncbi_taxon_id = ncbi_species_taxid
    ) %>%
    select(
      accno0, accno1, tdomain:tspecies, trank, ncbi_taxon_id, checkm_completeness, checkm_contamination, checkm_strain_heterogeneity,
      contig_count, genome_size, gtdb_genome_representative, gtdb_representative, l50_contigs, l50_scaffolds, longest_contig,
      longest_scaffold, mean_contig_length, mean_scaffold_length, n50_contigs, n50_scaffolds, ncbi_bioproject, ncbi_biosample,
      ncbi_genbank_assembly_accession, protein_count, scaffold_count
    ) 
} else {
  logmsg(sprintf("Reading NCBI taxonomy from %s", opt$options$taxflat), 'DEBUG')
  taxflat <- data.table::fread(opt$options$taxflat) %>%
    as_tibble() %>%
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
      taxflat %>% group_by(taxon, trank) %>% summarise(n = n(), .groups = 'drop_last') %>% ungroup() %>% filter(n > 1) %>%
        inner_join(taxflat %>% filter(tdomain == 'Eukaryota'), by = c('taxon', 'trank')),
      by = c('ncbi_taxon_id')
    ) 
}

# We will populate two tables, one with the full results, one with accessions
tblout <- tibble(
  accno = character(), profile = character(),
  evalue = double(), score = double(), bias = double()
)
accessions <- tibble(accno = character(), accto = character())

# Read all the tblout files
for ( tbloutfile in grep('\\.tblout', opt$args, value=TRUE) ) {
  if ( file.info(tbloutfile)$size == 0 ) {
    logmsg(sprintf("Skipping %s -- empty", tbloutfile), 'DEBUG')
  } else {
    logmsg(sprintf("Reading %s", tbloutfile), 'DEBUG')
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
    tblout <- union(tblout, t %>% select(accno, profile, evalue, score, bias))
    accessions <- union(
      accessions, 
      t %>% transmute(accno = accno, accto = sprintf("%s %s", accno, rest))
    )
  }
}

# Split the accto field into multiple rows
accessions <- accessions %>%
  separate_rows(accto, sep = "\x01") %>%
  mutate(
    taxon = ifelse(
      grepl('[^[]\\[(.*)\\]', accto), 
      sub('.*[^[]\\[(.*)\\].*', '\\1', accto),
      'unknown'
    ),
    accto = sub(' .*', '', accto)
  )

domtblout <- tibble(
  accno = character(), tlen = integer(), profile = character(), qlen = integer(), i = integer(), n = integer(), 
  dom_c_evalue = double(), dom_i_evalue = double(), dom_score = double(), dom_bias = double(),
  hmm_from = integer(), hmm_to = integer(), ali_from = integer(), ali_to = integer(), 
  env_from = integer(), env_to = integer()
)

# Read all the domtblout files
for ( domtbloutfile in grep('\\.domtblout', opt$args, value=TRUE) ) {
  if ( file.info(domtbloutfile)$size == 0 ) {
    logmsg(sprintf("Skipping %s -- empty", domtbloutfile), 'DEBUG')
  } else {
    logmsg(sprintf("Reading %s", domtbloutfile), 'DEBUG')
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
      )
    
    domtblout <- union(
      domtblout,
      t %>% select(
        accno, tlen, profile, qlen, i, n, dom_c_evalue, dom_i_evalue, dom_score, dom_bias,
        hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
      )
    )
  }
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
do_nextjoin <- function(dt) {
  # Filter dt on n > 1, then join with itself with transmuted column names
  t <- dt %>% 
    filter(n > 1) %>%
    left_join(
      dt %>% transmute(accno, profile, i = i - 1, next_row = TRUE, next_from = from, next_to = to), 
      by = c('accno', 'profile', 'i')
    ) %>%
    mutate(next_row = ifelse(is.na(next_row), FALSE, TRUE))
  # Set any NA next_row entries to FALSE (doesn't work through dtplyr)
  #t$next_row = ifelse(is.na(t$next_row), FALSE, TRUE)
  return(t)
}

logmsg("Calculating lengths -- this takes a long time!")
# FOR EACH FROM-TO PAIR:
for ( fs in list(
  c('hmm_from', 'hmm_to', 'hmmlen'),
  c('ali_from', 'ali_to', 'alilen'),
  c('env_from', 'env_to', 'envlen')
)) {
  logmsg(sprintf("Calculating %s", fs[3]), 'DEBUG')

  # 1. Set from and to to current pair, delete shorter rows with the same start and 
  # calculate i (rownumber) and n (total domains) for each combination of accno and profile
  #logmsg("\tCreating domt table", 'DEBUG')
  domt <- domtblout %>%
    transmute(accno = accno, profile = profile, from = get(fs[1]), to = get(fs[2])) %>%
    distinct()
  
  # Join with itself to get rid of shorter alignments starting at the same position as another
  # alignment; create row numbers for the from column and count the number of rows in each accno, 
  # profile group.
  #logmsg("\tSubset and create row number", 'DEBUG')
  domt <- domt %>%
    semi_join(
      domt %>% group_by(accno, profile, from) %>% filter(to == max(to)) %>% ungroup(),
      by = c('accno', 'profile', 'from', 'to')
    ) %>% 
    group_by(accno, profile) %>% mutate(i = row_number(from), n = n()) %>% ungroup()

  # This is where non-overlapping rows will be stored
  nooverlaps <- tibble(
    accno = character(), profile = character(),
    from = integer(), to = integer()
  )

  while ( nrow(domt) > 0 ) {
    #tic('dt')
    logmsg(sprintf("Working on overlaps, nrow: %d", domt %>% nrow()), 'DEBUG')
    
    # 2. Move rows with n == 1 to nooverlaps
    nooverlaps <- union(nooverlaps, domt %>% filter(n == 1) %>% select(accno, profile, from, to))
    domt <- domt %>% filter(n > 1)
    
    nextjoin <- do_nextjoin(domt)

    # As a debug aid, save the domt and nextjoin data sets
    #write_tsv(domt, 'domt.tsv.gz')
    #write_tsv(nextjoin, 'nextjoin.tsv.gz')
  
    # 3. Move rows that do not overlap with the next row
    nooverlaps <- union(nooverlaps, nextjoin %>% filter(to < next_from) %>% select(accno, profile, from, to))
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
    suppressWarnings( # The max(from.x) causes warnings: "no non-missing arguments to max; returning -Inf"; I can't find any errors
      nextjoin <- nextjoin %>%
        anti_join(
          nextjoin %>% select(accno, profile, from, to) %>% 
            inner_join(nextjoin %>% select(accno, profile, from, to), by = c('accno', 'profile', 'to')) %>% 
            filter(from.x != from.y) %>% 
            group_by(accno, profile, to) %>% summarise(from = max(from.x), .groups = 'drop_last') %>% ungroup(),
          by = c('accno', 'profile', 'from', 'to')
        )
    )

    # 7. Calculate a new domt from nextjoin
    if ( nrow(nextjoin) > 0 ) {
      domt <- nextjoin %>% distinct(accno, profile, from, to) %>%
        group_by(accno, profile) %>% mutate(i = rank(from), n = n()) %>% ungroup()
    }
    #toc()
  }
  
  # Add lengths to growing table
  lengths <- lengths %>%
    union(
      nooverlaps %>% mutate(len = to - from + 1) %>%
        group_by(accno, profile) %>% summarise(len = sum(len), from = min(from), to = max(to), .groups = 'drop_last') %>% ungroup() %>%
        gather(type, val, len, from, to) %>% # Change to pivot_longer()!
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
align_lengths <- domtblout %>% 
  distinct(accno, profile, tlen, qlen) %>%
  inner_join(lengths %>% spread(type, val, fill = 0), by = c('accno', 'profile'))

### dtfile <- "domtblout_after_calc_lens.feather"
### alfile <- "align_lengths_after_cal_lens.feather"
### logmsg(sprintf("Writing DEBUG tables: %s, %s", dtfile, alfile))
### write_feather(domtblout, dtfile)
### write_feather(align_lengths, alfile)

logmsg("Calculated lengths, inferring source databases from accession numbers", 'DEBUG')

if ( gtdb ) {
  accessions$db <- 'gtdb'
  accessions <- accessions %>%
    transmute(db, genome_accno = str_remove(taxon, '\\..*'), accno)
} else {
  # Infer databases from the structure of accession numbers
  accessions <- accessions %>%
    mutate(db = ifelse(grepl('^.._', accto), 'refseq', NA)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[0-9A-Z]{4,4}_[0-9A-Z]$', accto)), 'pdb', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^P[0-9]+\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[O,P,Q][0-9][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[ADEKOJMNP][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'genbank', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[C][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'embl', db)) %>%
    mutate(db = ifelse((is.na(db) & grepl('^[BFGIL][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'dbj', db))
}

logmsg("Inferred databases, calculating best scoring profile for each accession", 'DEBUG')

# Create proteins with entries from tblout not matching hmm_profile entries with rank == 'domain'.
# Calculate best scoring profile for each accession
proteins <- tblout %>% 
  anti_join(hmm_profiles %>% filter(prank == 'domain'), by = 'profile') %>%
  group_by(accno) %>% filter(score == max(score)) %>% ungroup() %>%
  select(accno, profile, score, evalue)
#proteins <- setorder(setDT(proteins), -score)[, head(.SD, 1), keyby = accno]

write_feather(proteins, "proteins_1st.feather")

logmsg("Calculated best scoring profiles, creating domains", 'DEBUG')

# Create table of domains as those that match domains specified in hmm_profiles
domains <- tblout %>%
  semi_join(hmm_profiles %>% filter(prank == 'domain'), by = 'profile') %>%
  select(accno, profile, score, evalue)

write_feather(domains, "domains_1st.feather")

# Join in lengths
logmsg(sprintf("Joining in lengths from domtblout, nrows before: %d", proteins %>% nrow()), 'DEBUG')
proteins <- proteins %>% inner_join(align_lengths, by = c('accno', 'profile'))
domains  <- domains %>% inner_join(align_lengths, by = c('accno', 'profile'))

write_feather(proteins, "proteins_2nd.feather")
write_feather(domains, "domains_2nd.feather")

logmsg("Joined in lengths, writing data", 'DEBUG')

# Subset data with hmm_mincov parameter
logmsg(sprintf("Subsetting output to proteins and domains covering at least %f of the hmm profile", opt$options$hmm_mincov), 'DEBUG')

# 1. proteins table
#logmsg(sprintf("proteins before: %d", proteins %>% nrow()), 'DEBUG')
p <- proteins %>% 
  inner_join(hmm_profiles %>% select(profile, plen), by = 'profile') %>%
  filter(hmmlen/plen >= opt$options$hmm_mincov) %>%
  select(-plen)
#logmsg(sprintf("proteins after: %d", proteins %>% nrow()), 'DEBUG')
p <- p %>%
  union(
    proteins %>% anti_join(p, by = 'accno') %>%
      semi_join(accessions %>% filter(db == 'pdb'), by = c('accno' = 'accno'))
  )
proteins <- p
rm(p)

write_feather(proteins, "proteins_3rd.feather")

# 1.b add pdb entries that are not present due to not passing the hmm_mincov criterion

# 2. domains
domains <- domains %>%
  inner_join(hmm_profiles %>% select(profile, plen), by = 'profile') %>%
  filter((hmm_to - hmm_from + 1)/plen >= opt$options$hmm_mincov) %>%
  select(-plen)

# 3. accessions
accessions <- accessions %>% 
  semi_join(
    union(proteins %>% select(accno), domains %>% select(accno)) %>% distinct(accno), 
    by = 'accno'
  )
accnov <- unique(accessions$accno)

# 4. tblout
tblout <- tblout %>% filter(accno %in% accnov)

# 5. domtblout
domtblout <- domtblout %>% filter(accno %in% accnov)

# Sequences in fasta file?
if ( opt$options$seqfaa != '' ) {
  logmsg(sprintf('Reading %s, splitting into %d sequences per batch', opt$options$seqfaa, ROWS_PER_SEQUENCE_TSV))
  cmd <- c(
    sprintf("%s %s", ifelse(grepl('\\.gz$', opt$options$seqfaa), 'zcat', 'cat'), opt$options$seqfaa),
    "sed '/^>/s/ .*//'",
    "gawk '/^>/ { printf(\"\\n%s\\n\",$$0); next; } { printf(\"%s\",$$0);} END { printf(\"\\n\");}'",
    "grep -v '^$'", 
    "sed 's/^>//'", 
    "paste - -",
    sprintf("split -a 3 --additional-suffix=.tsv -l %d - %s/sequences.", ROWS_PER_SEQUENCE_TSV, tempdir())
  )
  logmsg("Creating multiple tsv files with sequences", 'DEBUG')
  system(paste(cmd, collapse = '|'))

  sequences <- tibble(accno = character(), sequence = character())

  for ( f in Sys.glob(sprintf("%s/sequences.*.tsv", tempdir())) ) {
    logmsg(sprintf("Reading %s", f), 'DEBUG')
    sequences <- sequences %>%
      union(
        data.table::fread(f, col.names = c('accno', 'sequence'), header = FALSE) %>%
          as_tibble() %>%
          semi_join(accessions, by = 'accno')
      )
    logmsg(sprintf("\tCurrently %d sequences", nrow(sequences)), 'DEBUG')
  }
  logmsg(sprintf("Read %d sequences", nrow(sequences)), 'DEBUG')
}

# If we were called with the singletable option, prepare data suitable for that
if ( opt$options$singletable > '' ) {
  logmsg(sprintf("Writing single table format to %s", opt$options$singletable))

  # Join proteins with accessions and drop profile to get a single table output
  logmsg(sprintf("Joining in all accession numbers, nrows before: %d", proteins %>% nrow()), 'DEBUG')
  singletable <- proteins %>% 
    left_join(hmm_profiles, by='profile') %>% # Why is this a *left* join?
    inner_join(accessions, by='accno')
  
  if ( ! gtdb ) singletable <- singletable %>% mutate(accno = accto) 

  # Join in taxonomies, either GTDB or taxflat
  if ( gtdb ) {
    singletable <- union(
      singletable %>% 
        inner_join(gtdbtaxonomy, by = c('genome_accno' = 'accno0')) %>% select(-accno1),
      singletable %>% 
        anti_join(gtdbtaxonomy,  by = c('genome_accno' = 'accno0')) %>% 
        left_join(gtdbtaxonomy, by = c('genome_accno' = 'accno1')) %>% select(-accno0)
    )
    if ( singletable %>% filter(is.na(tspecies)) %>% nrow() > 0 ) {
      logmsg(
        sprintf("*** Accessions without GTDB species assignment: %s ***", singletable %>% filter(is.na(tspecies)) %>% pull(genome_accno) %>% paste(collapse = ', ')),
        "WARNING"
      )
    }
    logmsg(sprintf("Writing single table %s, nrows: %d", opt$options$singletable, singletable %>% nrow()), 'DEBUG')
    write_tsv(
      singletable %>% 
        select(db, accno, score, evalue, profile, psuperfamily:pgroup, genome_accno, tdomain:tspecies, tlen, qlen, alilen, hmmlen,envlen) %>%
        arrange(accno, profile),
      opt$options$singletable
    )
  } else {
    logmsg(sprintf("Adding NCBI taxon ids from taxflat, nrows before: %d", singletable %>% nrow()), 'DEBUG')
    singletable <- singletable %>% 
      left_join(
        taxflat %>% select(taxon, ncbi_taxon_id),
        by='taxon'
      )
    logmsg(sprintf("Writing single table %s, nrows: %d", opt$options$singletable, singletable %>% nrow()), 'DEBUG')
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
  logmsg('Copying to "accessions", creating indices', 'DEBUG')
  # I'm converting to tibble here, as I don't know how to do the separate_rows in
  # data.table and it's not in dtplyr.
  accessions <- accessions %>%
    arrange(db, taxon, accno, accto) %>%
    group_by(db, taxon, accno) %>%
    summarise(accto = paste(accto, collapse = ','), .groups = 'drop_last') %>%
    ungroup() %>%
    separate_rows(accto, sep = ',') %>%
    distinct()
}
if ( gtdb ) {
  taxa <- union(
    gtdbtaxonomy %>% semi_join(accessions, by = c('accno0' = 'genome_accno')) %>% 
      select(-accno1) %>% rename(genome_accno = accno0),
    gtdbtaxonomy %>% anti_join(accessions, by = c('accno0' = 'genome_accno')) %>% 
      mutate(genome_accno = ifelse(accno1 == 'none', accno0, accno1)) %>%
      select(-accno0, -accno1)
  )
} else {
  taxa <- taxflat %>% semi_join(accessions %>% distinct(taxon), by='taxon')
}

# Write information about missing genomes, if asked for
if ( length(grep('missing', names(opt$options), value = TRUE)) > 0 & str_length(opt$options$missing) > 0 ) {
  logmsg(sprintf("Writing information about missing genomes to '%s'", opt$options$missing), 'DEBUG')
  write(
    sprintf(
      "Genomes missing from GTDB metadata:\n%s",
      paste(accessions %>% anti_join(taxa, by = 'genome_accno') %>% arrange(genome_accno) %>% pull(genome_accno), collapse = '\n')
    ),
    opt$options$missing
  )
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
  
  logmsg('Copying to "proteins", creating indices', 'DEBUG')
  con %>% copy_to(proteins, 'proteins', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE INDEX "proteins.i00" ON "proteins"("accno");')
  con %>% DBI::dbExecute('CREATE INDEX "proteins.i01" ON "proteins"("profile");')

  logmsg('Copying to "domains", creating indices', 'DEBUG')
  con %>% copy_to(domains, 'domains', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE INDEX "domains.i00" ON "domains"("accno");')
  con %>% DBI::dbExecute('CREATE INDEX "domains.i01" ON "domains"("profile");')

  logmsg('Copying to "hmm_profiles", creating indices', 'DEBUG')
  con %>% copy_to(hmm_profiles, 'hmm_profiles', temporary = FALSE, overwrite = TRUE)
  con %>% DBI::dbExecute('CREATE UNIQUE INDEX "hmm_profiles.i00" ON "hmm_profiles"("profile");')

  logmsg('Copying to "taxa", creating indices', 'DEBUG')
  con %>% copy_to(taxa, 'taxa', temporary = FALSE, overwrite = TRUE)
  if ( gtdb ) {
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i00" ON "taxa"("genome_accno");')
  } else {
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i00" ON "taxa"("taxon", "trank");')
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "taxa.i01" ON "taxa"("ncbi_taxon_id");')
  }

  logmsg('Saving tblout and domtblout to database', 'DEBUG')
  con %>% copy_to(tblout    %>% arrange(accno, profile),    'tblout',    temporary = FALSE, overwrite = TRUE)
  con %>% copy_to(domtblout %>% arrange(accno, profile, i), 'domtblout', temporary = FALSE, overwrite = TRUE)

  if ( ! gtdb ) {
  logmsg(sprintf('Creating dupfree_proteins, using %d as fuzzy factor', opt$options$fuzzy_factor), 'DEBUG')
    dp <- proteins %>% inner_join(accessions %>% transmute(accno = accto, db, taxon), by = 'accno') %>%
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
    logmsg(sprintf('Saving sequences table', opt$options$seqfaa), 'DEBUG')
    con %>% copy_to(sequences, 'sequences',
      temporary = FALSE, overwrite = TRUE
    )
    con %>% DBI::dbExecute('CREATE UNIQUE INDEX "sequences.i00" ON "sequences"("accno");')
  }

  logmsg('Disconnecting from sqlite3 db', 'DEBUG')
  con %>% DBI::dbDisconnect()
}

if ( length(grep('featherprefix', names(opt$options), value = TRUE)) > 0 & str_length(opt$options$featherprefix) > 0 ) {
  logmsg(sprintf("Writing data to feather tables prefixed with %s", opt$options$featherprefix))
  logmsg("Writing dbsources to feather file", 'DEBUG')
  write_feather(
    tibble(source = dbsource[1], name = dbsource[2], version = dbsource[3]),
    sprintf("%s.dbsources.feather", opt$options$featherprefix)
  )

  logmsg("Writing accessions to feather file", 'DEBUG')
  write_feather(accessions %>% arrange(accno), sprintf("%s.accessions.feather", opt$options$featherprefix))

  logmsg("Writing proteins to feather file", 'DEBUG')
  write_feather(proteins %>% arrange(accno), sprintf("%s.proteins.feather", opt$options$featherprefix))

  logmsg("Writing domains to feather file", 'DEBUG')
  write_feather(domains %>% arrange(accno), sprintf("%s.domains.feather", opt$options$featherprefix))

  logmsg("Writing hmm_profiles to feather file", 'DEBUG')
  write_feather(hmm_profiles %>% arrange(profile), sprintf("%s.hmm_profiles.feather", opt$options$featherprefix))

  logmsg("Writing taxa to feather file", 'DEBUG')
  write_feather(taxa, sprintf("%s.taxa.feather", opt$options$featherprefix))

  logmsg("Writing tblout to feather file", 'DEBUG')
  write_feather(tblout, sprintf("%s.tblout.feather", opt$options$featherprefix))

  logmsg("Writing domtblout to feather file", 'DEBUG')
  write_feather(domtblout %>% arrange(accno), sprintf("%s.domtblout.feather", opt$options$featherprefix))

  if ( opt$options$seqfaa != '' ) {
    logmsg("Writing sequences to feather file", 'DEBUG')
    write_feather(sequences %>% arrange(accno), sprintf("%s.sequences.feather", opt$options$featherprefix))
  }
}

logmsg("Done")
