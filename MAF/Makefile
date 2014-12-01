MAF_DIR=mafs
BAM_DIR=bams
BYSAMPLE_DIR=by-sample
LOG_DIR=logs
REF=/extscratch/morinlab/shared/rmccloskey/GRCh37-lite.fa

all: $(patsubst $(MAF_DIR)/%, $(BYSAMPLE_DIR)/%, $(wildcard $(MAF_DIR)/*))

RUNSCRIPT=echo "`date`: starting $(notdir $^)"; \
	mkdir -p $(dir $@); \
	mkdir -p $(LOG_DIR); \
	qsub -V -sync yes -cwd -S /bin/sh -N x$(notdir $^) \
		 -o $(LOG_DIR)/$(notdir $^).stdout \
		 -e $(LOG_DIR)/$(notdir $^).stderr \
		 -pe ncpus $(NCPUS) \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^ && \
	echo "`date`: finished $(notdir $^)"

$(BYSAMPLE_DIR)/%.maf: NCPUS := 1
$(BYSAMPLE_DIR)/%.maf: MEM := 1G
$(BYSAMPLE_DIR)/%.maf: $(BYSAMPLE_DIR)/%.sh
	$(RUNSCRIPT)

$(BYSAMPLE_DIR)/%.sh: $(MAF_DIR)/%.maf
	$(eval SAMPLE := $(shell basename $@ | cut -d '.' -f 1))
	echo ./augment_maf.py \
	$(addprefix --normal-bam , $(wildcard $(BAM_DIR)/$(SAMPLE)/normal/*.bam)) \
	$(addprefix --tumour-bam , $(wildcard $(BAM_DIR)/$(SAMPLE)/tumour/*.bam)) \
	$^ $(REF) $(patsubst %.sh, %.maf, $@) > $@
