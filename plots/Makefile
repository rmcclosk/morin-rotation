all: $(patsubst mafs/%.snvs.maf, VAF_distributions/%.pdf, $(wildcard mafs/*snvs.maf))

VAF_distributions/%.pdf: mafs/%.snvs.maf
	./VAF_distribution.R \
		-m $^ \
		-p $(shell grep `basename $@ | cut -d '.' -f 1` purity.tsv | cut -f 2) \
		-o $@ \
		-t $(shell basename $@ | cut -d '.' -f 1)
