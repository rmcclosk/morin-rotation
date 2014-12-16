TITAN_DIR=titan
INTERM_DIR=intermediate
BAM_DIR=bams
LOG_DIR=logs
INTERVALS=/extscratch/morinlab/shared/common/SureSelect_regions.list
NCPUS=4

ALL_SAMPLES = $(wildcard $(BAM_DIR)/*)
ALL_PARAMS = $(patsubst $(BAM_DIR)/%, \
					    $(TITAN_DIR)/%_params.txt.bz2, \
                	  	$(ALL_SAMPLES))

ALL_SEGS = $(patsubst $(BAM_DIR)/%, \
					  $(TITAN_DIR)/%_segs.txt.bz2, \
                	  $(ALL_SAMPLES))

ALL_BEDS = $(patsubst $(BAM_DIR)/%, \
					  $(TITAN_DIR)/%.bed.bz2, \
                	  $(ALL_SAMPLES))

MAKESCRIPT=mkdir -p $(dir $@); \
	echo "\#!/bin/bash" > $@; \
	echo $(SCRIPT) >> $@

RUNSCRIPT=echo "`date`: starting $(notdir $^)"; \
	mkdir -p $(dir $@); \
	mkdir -p $(patsubst $(INTERM_DIR)/%/, $(LOG_DIR)/%/, $(dir $@)); \
	qsub -V -sync yes -cwd -S /bin/sh -N $(notdir $^) \
		 -o $(patsubst $(INTERM_DIR)/%.sh, $(LOG_DIR)/%.stdout, $^) \
		 -e $(patsubst $(INTERM_DIR)/%.sh, $(LOG_DIR)/%.stderr, $^) \
		 -pe ncpus $(NCPUS) \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^ && \
	echo "`date`: finished $(notdir $^)"

all: titan.seg $(ALL_BEDS)

# make BED files
%.bz2: %
	bzip2 $^

%.bed: %_segs.txt.bz2
	bzcat $^ | tail -n +2 | awk '{print $$2"\t"$$3"\t"$$4"\tcopy.number\t"$$10}' > $@

# compile into one file
titan.seg: purity.tsv
	rm -f $@
	echo -e "sample\tchr\tstart\tend\tcopy.number\tprevalence" > $@
	for SEGS in $(ALL_SEGS); do \
		SAMPLE=`basename $$SEGS | cut -d '_' -f 1`; \
		bzcat $$SEGS | tail -n +2 | sed s/tumour/$$SAMPLE/ | cut -f 1-4,10,14 >> $@; \
	done

purity.tsv: $(ALL_PARAMS)
	rm -f $@
	for PARAMS in $^; do \
		SAMPLE=`basename $$PARAMS | cut -d '_' -f 1`; \
		CONTAM=`bzgrep contamination $$PARAMS | cut -f 2`; \
		PURITY=`echo 1-$$CONTAM | bc`; \
		echo -e "$$SAMPLE\t$$PURITY" >> purity.tsv; \
	done

# pull out best cluster
$(TITAN_DIR)/%_params.txt.bz2: $(INTERM_DIR)/%/titan
	mkdir -p $(dir $@)
	./parse_titan.py $^ $(patsubst %_params.txt.bz2, %, $@)

# run titan
$(INTERM_DIR)/%/titan: $(INTERM_DIR)/%/input.txt $(INTERM_DIR)/%/tumour.bam $(INTERM_DIR)/%/normal.bam $(INTERM_DIR)/%/tumour.bam.bai $(INTERM_DIR)/%/normal.bam.bai
	echo "`date`: started TITAN"
	python TITANRunner/TitanRunner.py -i $(word 1, $^) \
                                      --project-name titan \
                                      --project-path $(dir $@) \
                                      -c config.cfg
	echo "`date`: finished TITAN"
	mv $@_R* $@
	rm $(word 2, $^)
	rm $(word 3, $^)

# make titan input
$(INTERM_DIR)/%/input.txt: $(INTERM_DIR)/%/tumour.bam $(INTERM_DIR)/%/normal.bam $(INTERM_DIR)/%/tumour.bam.bai $(INTERM_DIR)/%/normal.bam.bai
	echo -e "tumour\ttumour\t$(word 1, $^)\tnormal\tnormal\t$(word 2, $^)" > $@

# index bams
$(INTERM_DIR)/%.bam.bai: MEM := 1G
$(INTERM_DIR)/%.bam.bai: NCPUS := 1
$(INTERM_DIR)/%.bam.bai: $(INTERM_DIR)/%_index.sh
	$(RUNSCRIPT)

$(INTERM_DIR)/%_index.sh: $(INTERM_DIR)/%.bam
	echo "#!/bin/bash" > $@
	echo samtools index $^ >> $@

# pull out exome reads and merge bams
$(INTERM_DIR)/%_extract.sh: NCPUS := 4
$(INTERM_DIR)/%_extract.sh: $(BAM_DIR)/%
	mkdir -p $(dir $@)
	echo "#!/bin/bash" > $@
	for BAM in $(wildcard $^/*.bam); do \
		echo "./extract_exome.py $$BAM $(INTERVALS) `echo $$BAM | sed s/$(BAM_DIR)/$(INTERM_DIR)/` $(NCPUS)" >> $@; \
	done
	if [[ $(words $(wildcard $^/*.bam)) -eq 1 ]]; then \
		echo mv $(patsubst $(BAM_DIR)/%.bam, $(INTERM_DIR)/%.bam,  $(wildcard $^/*.bam)) $(patsubst $(INTERM_DIR)/%_extract.sh, $(INTERM_DIR)/%.bam, $@) >> $@; \
	else \
		echo samtools merge -@ $(NCPUS) $(patsubst $(INTERM_DIR)/%_extract.sh, $(INTERM_DIR)/%.bam, $@) $(patsubst $(BAM_DIR)/%.bam, $(INTERM_DIR)/%.bam,  $(wildcard $^/*.bam)) >> $@; \
	fi

$(INTERM_DIR)/%.bam: NCPUS := 4
$(INTERM_DIR)/%.bam: MEM := 1G
$(INTERM_DIR)/%.bam: $(INTERM_DIR)/%_extract.sh
	mkdir -p $(basename $@)
	mkdir -p $(dir $@)
	$(RUNSCRIPT)
	rm -rf $(patsubst $(INTERM_DIR)/%.bam, $(INTERM_DIR)/%/*.bam, $@)
