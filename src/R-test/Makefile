DIFF = grep -v "^[ ]*$$" $@.out | diff $@.expect -

CHECK = if [ ` $(DIFF) | wc -l ` -eq 0 ]; then echo "*** $@ OK ***"; else echo "*** $@ FAILED ***"; echo "Here's a diff: "; $(DIFF); exit 1; fi

SQLITE_SELECT = sqlite3 $@.sqlite3 "select accno, accto, taxon, db from accessions order by accno, accto, taxon, db; select accno, profile, score, evalue, tlen, qlen, alilen from proteins order by accno, profile; select profile, psuperfamily, pfamily, pclass, psubclass, pgroup, prank, version, plen from hmm_profiles order by profile; select ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, taxon, trank from taxa order by ncbi_taxon_id; SELECT source, name, version FROM dbsources; SELECT accno, profile, i, n, dom_c_evalue, dom_i_evalue, dom_score, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to FROM domains ORDER BY accno, i; SELECT * FROM tblout ORDER BY accno, profile; SELECT * FROM domtblout ORDER BY accno, profile, i;" > $@.out

all: hmmsearch2classification.00 hmmsearch2classification.01 hmmsearch2classification.02 hmmsearch2classification.03

all.sqlite: hmmsearch2classification.02 hmmsearch2classification.03

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
	@$(SQLITE_SELECT)
	@$(CHECK)

hmmsearch2classification.03:
	@rm -f $@.sqlite3
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=hmmsearch2classification.taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT)
	@$(CHECK)

# Data for this tc available by request from daniel.lundin@dbb.su.se
hmmsearch2classification.large_tc:
	../R/hmmsearch2classification.r --verbose --dbsource NCBI:NR:20180109 --profilehierarchies=$@.phier.tsv --sqlitedb=$@.sqlite3 --taxflat=taxflat.tsv $@.d/*.tblout $@.d/*.domtblout
	@$(SQLITE_SELECT)
	@$(CHECK)
