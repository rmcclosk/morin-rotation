SEG_DIR=segments
ABS_DIR=absolute
MAF_DIR=mafs
RUN=qsub -V -sync yes -cwd -S /bin/sh -o /dev/null -e /dev/null -pe ncpus $(NCPUS) -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) -N $(NAME) $(COMMAND)
PLATFORM=Illumina_WES
PRIMARY_DISEASE=colorectal_cancer
SIGMA_P=0.01
MAX_SIGMA_H=0.03
MIN_PLOIDY=1
MAX_PLOIDY=6
MAX_AS_SEG_COUNT=1500
MAX_NON_CLONAL=0.1
MAX_NEG_GENOME=0.05
MIN_MUT_AF=0.05

ifdef MAF_DIR
	ABS_DEPS = $(SEG_DIR)/%.seg.bz2 $(MAF_DIR)/%.maf.bz2
else
	ABS_DEPS = $(SEG_DIR)/%.seg.bz2
endif

all: summary.PP-calls_tab.txt

summary.PP-calls_tab.txt: MEM := 1G
summary.PP-calls_tab.txt: NCPUS := 1
summary.PP-calls_tab.txt: summary.sh
	@$(eval COMMAND=$^)
	@$(eval NAME=summary)
	$(RUN)

ALL_RDATA = $(patsubst $(SEG_DIR)/%.seg.bz2,  \
                       $(ABS_DIR)/%.ABSOLUTE.RData, \
                       $(wildcard $(SEG_DIR)/*.seg.bz2))
summary.sh:	$(ALL_RDATA)
	@echo    "#!/bin/sh" > $@
	@echo -n "summarize.R " >> $@
	@echo -n "-p $(ABS_DIR) -n summary -d $(shell pwd) " >> $@

$(ABS_DIR)/%.ABSOLUTE.RData: MEM := 1G
$(ABS_DIR)/%.ABSOLUTE.RData: NCPUS := 1
$(ABS_DIR)/%.ABSOLUTE.RData: $(ABS_DIR)/%.sh
	@$(eval COMMAND=$^)
	@$(eval NAME=absolute-$(notdir $(basename $^)))
	$(RUN)

$(ABS_DIR)/%.sh: $(ABS_DEPS)
	@echo    "#!/bin/sh" > $@
	@echo -n "ABSOLUTE.R " >> $@
	@echo -n "-s $(word 1, $^) -n $(PLATFORM) -o $(shell pwd)/$(ABS_DIR) -d $(PRIMARY_DISEASE) " >> $@
	@echo -n "-p $(SIGMA_P) -j $(MAX_SIGMA_H) -i $(MIN_PLOIDY) " >> $@
	@echo -n "-a $(MAX_PLOIDY) -c $(MAX_AS_SEG_COUNT) " >> $@
	@echo -n "-l $(MAX_NON_CLONAL) -x $(MAX_NEG_GENOME) " >> $@
ifdef MAF_DIR
	@echo -n "-m $(word 2, $^) -f $(MIN_MUT_AF) " >> $@
endif
