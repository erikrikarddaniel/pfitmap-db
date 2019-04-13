DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

SQLITE_SELECT_CLASSIFICATION = sqlite3 $@.sqlite3 "select '--> accessions <--'; select accno, accto, taxon, db from accessions order by accno, accto, taxon, db; select '--> proteins <--'; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from proteins order by accno, profile; select '--> dupfree_proteins <--'; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from dupfree_proteins order by accno, profile; select '--> hmm_profiles <--'; select profile, psuperfamily, pfamily, pclass, psubclass, pgroup, prank, version, plen from hmm_profiles order by profile; select '--> taxa <--'; select ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, taxon, trank from taxa order by ncbi_taxon_id; select '--> dbsources <--'; SELECT source, name, version FROM dbsources; select '--> domains <--'; SELECT accno, profile, i, n, dom_c_evalue, dom_i_evalue, dom_score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to FROM domains ORDER BY accno, profile, i; select '--> tblout <--'; SELECT * FROM tblout ORDER BY accno, profile; select '--> domtblout <--'; SELECT * FROM domtblout ORDER BY accno, profile, i;" > $@.out
SQLITE_SELECT_FETCHSEQS = sqlite3 $@.sqlite3 "SELECT accno, sequence FROM sequences ORDER BY accno;" > $@.out
#SQLITE_SELECT_CLASSIFICATION = sqlite3 $@.sqlite3 "select accno, accto, taxon, db from accessions order by accno, accto, taxon, db; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from proteins order by accno, profile; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from dupfree_proteins order by accno, profile; select accno_from, accno_to from dupfree_accessions order by accno, accto; select profile, psuperfamily, pfamily, pclass, psubclass, pgroup, prank, version, plen from hmm_profiles order by profile; select ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, taxon, trank from taxa order by ncbi_taxon_id; SELECT source, name, version FROM dbsources; SELECT accno, profile, i, n, dom_c_evalue, dom_i_evalue, dom_score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to FROM domains ORDER BY accno, profile, i; SELECT * FROM tblout ORDER BY accno, profile; SELECT * FROM domtblout ORDER BY accno, profile, i;" > $@.out

DB2FEATHER_OUT = for f in $@.*.feather; do Rscript --default-packages=dplyr,feather -e "print(paste0('$$f', ': ', read_feather('$$f') %>% nrow(), ' rows'))"; done > $@.out

all: pf-classify pf-db2feather pf-fetchseqs

clean:
	rm -f *.out

all.sqlite: pf-classify.02 pf-classify.03 pf-classify.04

pf-classify: pf-classify.00 pf-classify.01 pf-classify.02 pf-classify.03 pf-classify.04 pf-classify.05 pf-classify.06 pf-classify.07 pf-classify.08

pf-db2feather: pf-db2feather.00 pf-db2feather.01 pf-db2feather.02 pf-db2feather.03

pf-fetchseqs: pf-fetchseqs.00 pf-fetchseqs.01 pf-fetchseqs.02 pf-fetchseqs.03 pf-fetchseqs.04 pf-fetchseqs.05 pf-fetchseqs.06 pf-fetchseqs.07 pf-fetchseqs.08

pf-fasta-unique-taxon-protein: pf-fasta-unique-taxon-protein.00 pf-fasta-unique-taxon-protein.01

Nrd.test.tar.gz:
	scp rnrdb:projects/profiles/profiles/Nrd.test.tar.gz .

pf-classify.00:
	../R/pf-classify.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(CHECK)

pf-classify.01:
	../R/pf-classify.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(CHECK)

pf-classify.02:
	@rm -f $@.sqlite3
	../R/pf-classify.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

pf-classify.03:
	@rm -f $@.sqlite3
	../R/pf-classify.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

# Test the fuzzy_factor to produce a taxon reduced protein list
pf-classify.04:
	@rm -f $@.sqlite3
	../R/pf-classify.r --verbose --fuzzy_factor=30 --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

# Check that the script fails when no dbsource is specified
pf-classify.05:
	( ../R/pf-classify.r --verbose --dbsource x --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout 2>&1 | grep -o 'needs to contain three components' > $@.out )
	@$(CHECK)

pf-classify.06:
	( ../R/pf-classify.r --verbose --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout 2>&1 | grep -o 'dbsource is required' > $@.out )
	@$(CHECK)

pf-classify.07:
	@rm -f $@.sqlite3
	../R/pf-classify.r --hmm_mincov=0.95 --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

# Thought I had problems with the 4CON structure...
pf-classify.08:
	@rm -f $@.sqlite3
	../R/pf-classify.r --hmm_mincov=0.95 --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=pf-classify.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

pf-db2feather.00: 
	@rm -f $@*.feather
	@../R/pf-classify.r --dbsource NCBI:NR:20180109 --profilehierarchies=pf-classify.02.phier.tsv --sqlitedb=pf-classify.02.sqlite3 --taxflat=pf-classify.taxflat.tsv pf-classify.02.d/*.tblout pf-classify.02.d/*.domtblout
	@../R/pf-db2feather.r --verbose --prefix=$@ pf-classify.02.sqlite3
	@$(DB2FEATHER_OUT)
	@$(CHECK)

pf-db2feather.01: 
	@rm -f $@*.feather
	@../R/pf-classify.r --dbsource NCBI:NR:20180109 --profilehierarchies=pf-classify.03.phier.tsv --sqlitedb=pf-classify.03.sqlite3 --taxflat=pf-classify.taxflat.tsv pf-classify.03.d/*.tblout pf-classify.03.d/*.domtblout
	@../R/pf-db2feather.r --verbose --prefix=$@ pf-classify.03.sqlite3
	@$(DB2FEATHER_OUT)
	@$(CHECK)

pf-db2feather.02:
	@rm -f $@*.feather
	@../R/pf-db2feather.r --verbose --prefix=$@ $@.sqlite3
	@$(DB2FEATHER_OUT)
	@$(CHECK)

pf-db2feather.03:
	@rm -f $@*.feather
	@../R/pf-db2feather.r --verbose --dbs=refseq,pdb --prefix=$@ $@.sqlite3
	@$(DB2FEATHER_OUT)
	@$(CHECK)

pf-fetchseqs.00:
	@cp $@.original.sqlite3 $@.sqlite3
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3
	@$(SQLITE_SELECT_FETCHSEQS)
	@$(CHECK)
	
pf-fetchseqs.01:
	@cp $@.original.sqlite3 $@.sqlite3
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3
	@$(SQLITE_SELECT_FETCHSEQS)
	@$(CHECK)
	
pf-fetchseqs.02:
	@cp $@.original.sqlite3 $@.sqlite3
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3 $@.tsv $@.faa $@.fasta
	@$(SQLITE_SELECT_FETCHSEQS)
	@$(CHECK)
	
pf-fetchseqs.03:
	@cp $@.original.sqlite3 $@.sqlite3
	@../R/pf-fetchseqs.r --verbose --sourcedbs=refseq,pdb --sqlitedb=$@.sqlite3 $@.tsv $@.faa $@.fasta
	@$(SQLITE_SELECT_FETCHSEQS)
	@$(CHECK)

pf-fetchseqs.04:
	@cp $@.original.sqlite3 $@.sqlite3
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3 --prefetch_accnos=$@.prefetch.accnos --postfetch_accnos=$@.postfetch.accnos
	@$(SQLITE_SELECT_FETCHSEQS)
	@cat $@.prefetch.accnos >> $@.out
	@cat $@.postfetch.accnos >> $@.out
	@$(CHECK)

pf-fetchseqs.05:
	@cp $@.original.sqlite3 $@.sqlite3
	@../R/pf-fetchseqs.r --verbose --skipfetch --sqlitedb=$@.sqlite3
	@$(SQLITE_SELECT_FETCHSEQS)
	@$(CHECK)

pf-fetchseqs.06:
	@cp $@.original.sqlite3 $@.sqlite3
	@rm -f $@.out
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3 --fetchedseqs=$@.tsv
	@mv $@.tsv $@.out
	@$(CHECK)

pf-fetchseqs.07:
	@cp $@.original.sqlite3 $@.sqlite3
	@rm -f $@.out
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3 --sourcedbs=refseq --fetchedseqs=$@.feather --faalevel=psuperfamily --faadir=$@.psuperfamilies
	@cat $@.psuperfamilies/*.faa > $@.out
	@ls $@.feather >> $@.out
	@$(CHECK)

pf-fetchseqs.08:
	@cp $@.original.sqlite3 $@.sqlite3
	@rm -f $@.out
	@../R/pf-fetchseqs.r --verbose --sqlitedb=$@.sqlite3 --sourcedbs=refseq --fetchedseqs=$@.feather --faalevel=psuperfamily --faadir=$@.psuperfamilies --faahmmcov=1.0
	@cat $@.psuperfamilies/*.faa > $@.out
	@ls $@.feather >> $@.out
	@$(CHECK)

pf-fasta-unique-taxon-protein.00: pf-db2feather.01
	../R/pf-fasta-unique-taxon-protein.r --verbose --featherprefix=$< --trank=species --prank=subclass $@.faa > $@.out
	@$(CHECK)

pf-fasta-unique-taxon-protein.01: 
	../R/pf-fasta-unique-taxon-protein.r --verbose --featherprefix=$@ --trank=species --prank=subclass $@.faa > $@.out
	@$(CHECK)
