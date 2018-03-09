DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

SQLITE_SELECT_CLASSIFICATION = sqlite3 $@.sqlite3 "select '--> accessions <--'; select accno, accto, taxon, db from accessions order by accno, accto, taxon, db; select '--> proteins <--'; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from proteins order by accno, profile; select '--> dupfree_proteins <--'; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from dupfree_proteins order by accno, profile; select '--> hmm_profiles <--'; select profile, psuperfamily, pfamily, pclass, psubclass, pgroup, prank, version, plen from hmm_profiles order by profile; select '--> taxa <--'; select ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, taxon, trank from taxa order by ncbi_taxon_id; select '--> dbsources <--'; SELECT source, name, version FROM dbsources; select '--> domains <--'; SELECT accno, profile, i, n, dom_c_evalue, dom_i_evalue, dom_score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to FROM domains ORDER BY accno, profile, i; select '--> tblout <--'; SELECT * FROM tblout ORDER BY accno, profile; select '--> domtblout <--'; SELECT * FROM domtblout ORDER BY accno, profile, i;" > $@.out
SQLITE_SELECT_FETCHSEQS = sqlite3 $@.sqlite3 "SELECT accno, sequence FROM sequences ORDER BY accno;" > $@.out
#SQLITE_SELECT_CLASSIFICATION = sqlite3 $@.sqlite3 "select accno, accto, taxon, db from accessions order by accno, accto, taxon, db; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from proteins order by accno, profile; select accno, profile, score, evalue, tlen, qlen, hmmlen, alilen, envlen, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to from dupfree_proteins order by accno, profile; select accno_from, accno_to from dupfree_accessions order by accno, accto; select profile, psuperfamily, pfamily, pclass, psubclass, pgroup, prank, version, plen from hmm_profiles order by profile; select ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, taxon, trank from taxa order by ncbi_taxon_id; SELECT source, name, version FROM dbsources; SELECT accno, profile, i, n, dom_c_evalue, dom_i_evalue, dom_score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to FROM domains ORDER BY accno, profile, i; SELECT * FROM tblout ORDER BY accno, profile; SELECT * FROM domtblout ORDER BY accno, profile, i;" > $@.out

all: hmmsearch2classification pf-db2feather pf-fetchseqs

all.sqlite: hmmsearch2classification.02 hmmsearch2classification.03 hmmsearch2classification.04

hmmsearch2classification: hmmsearch2classification.00 hmmsearch2classification.01 hmmsearch2classification.02 hmmsearch2classification.03 hmmsearch2classification.04

pf-db2feather: pf-db2feather.00 pf-db2feather.01

pf-fetchseqs: pf-fetchseqs.00 pf-fetchseqs.01

Nrd.test.tar.gz:
	scp rnrdb:projects/profiles/profiles/Nrd.test.tar.gz .

hmmsearch2classification.00:
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=hmmsearch2classification.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(CHECK)

hmmsearch2classification.01:
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --singletable=$@.out --taxflat=hmmsearch2classification.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(CHECK)

hmmsearch2classification.02:
	@rm -f $@.sqlite3
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=hmmsearch2classification.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

hmmsearch2classification.03:
	@rm -f $@.sqlite3
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=hmmsearch2classification.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

hmmsearch2classification.04:
	@rm -f $@.sqlite3
	../R/hmmsearch2classification.r --verbose --fuzzy_factor=30 --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=hmmsearch2classification.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)

pf-db2feather.00: 
	@../R/hmmsearch2classification.r --dbsource NCBI:NR:20180109 --profilehierarchies=hmmsearch2classification.02.phier.tsv --sqlitedb=hmmsearch2classification.02.sqlite3 --taxflat=hmmsearch2classification.taxflat.tsv hmmsearch2classification.02.d/*.tblout hmmsearch2classification.02.d/*.domtblout
	@../R/pf-db2feather.r --verbose --prefix=$@ hmmsearch2classification.02.sqlite3
	@ls $@.*.feather > $@.out
	@$(CHECK)

pf-db2feather.01: 
	@../R/hmmsearch2classification.r --dbsource NCBI:NR:20180109 --profilehierarchies=hmmsearch2classification.03.phier.tsv --sqlitedb=hmmsearch2classification.03.sqlite3 --taxflat=hmmsearch2classification.taxflat.tsv hmmsearch2classification.03.d/*.tblout hmmsearch2classification.03.d/*.domtblout
	@../R/pf-db2feather.r --verbose --prefix=$@ hmmsearch2classification.03.sqlite3
	@ls $@.*.feather > $@.out
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
	
# Data for this tc available by request from daniel.lundin@dbb.su.se
hmmsearch2classification.large_tc:
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT_CLASSIFICATION)
	@$(CHECK)
