NORM_DIR=coverage-normal
TUM_DIR=coverage-tumor
SEG_DIR=segments
ADMIX_DIR=admixture
LOG_DIR=logs

READ_LENGTH=100
CHR_LIST="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X"
NCPUS=4
MEM=1G

all: $(patsubst $(TUM_DIR)/%.sample_interval_summary.bz2, \
				$(SEG_DIR)/%.seg.bz2, \
                $(wildcard $(TUM_DIR)/*.sample_interval_summary.bz2))

$(SEG_DIR)/%.seg.bz2: $(SEG_DIR)/%.seg
	bzip2 $^

$(SEG_DIR)/%.seg: $(SEG_DIR)/%.segment.lrr.txt.bz2
	echo -e "Chromosome\tStart\tEnd\tNum_Probes\tSegment_Mean" > $@
	bzcat $^ | awk '{print substr($$1, 4)"\t"$$2"\t"$$3"\t"($$3-$$2)/1000"\t"$$4}' >> $@

$(SEG_DIR)/%.segment.copynumber.txt.bz2: $(SEG_DIR)/%.sh
	qsub -V -sync yes -cwd -S /bin/sh -N ExomeCNV -pe ncpus $(NCPUS) \
	     -o $(LOG_DIR)/$(notdir $(basename $^)).stdout \
		 -e $(LOG_DIR)/$(notdir $(basename $^)).stderr \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^

$(SEG_DIR)/%.sh: $(NORM_DIR)/%.sample_interval_summary.bz2 $(TUM_DIR)/%.sample_interval_summary.bz2 $(ADMIX_DIR)/%.txt
	@echo    "#!/bin/sh" > $@
	@echo -n "ExomeCNV.R " >> $@
	@echo -n "-a $(word 1, $^) " >> $@
	@echo -n "-b $(word 2, $^) " >> $@
	@echo -n "-c $(READ_LENGTH) " >> $@
	@echo -n "-d $(shell cat $(word 3, $^)) " >> $@
	@echo -n "-e $(CHR_LIST) " >> $@
	@echo -n "-p $(basename $@).segment.copynumber.txt.bz2 " >> $@
	@echo -n "-q $(basename $@).segment.lrr.txt.bz2 " >> $@
	@echo    "-s $(NCPUS)" >> $@
