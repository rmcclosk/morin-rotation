TITAN_DIR=titan
INTERM_DIR=intermediate
NORM_DIR=normal-bam
TUM_DIR=tumour-bam
LOG_DIR=logs

ALL_PARAMS = $(patsubst $(TUM_DIR)/%, \
					    $(TITAN_DIR)/%_params.txt.bz2, \
                	  	$(wildcard $(TUM_DIR)/*))

ALL_SEGS = $(patsubst $(TUM_DIR)/%, \
					  $(TITAN_DIR)/%_segs.txt.bz2, \
                	  $(wildcard $(TUM_DIR)/*))

#.SECONDARY: $(wildcard $(INTERM_DIR)/*) $(ALL_SEGS) $(ALL_PARAMS)

all: purity.tsv segments.tsv

# merge bams
$(INTERM_DIR)/%_normal_merge.sh: $(NORM_DIR)/%
	echo "#!/bin/bash" > $@
	echo samtools merge $(patsubst $(INTERM_DIR)/%_normal_merge.sh, $(INTERM_DIR)/%_normal.bam, $@) $(wildcard $^/*) >> $@

$(INTERM_DIR)/%_normal.bam: MEM := 1G
$(INTERM_DIR)/%_normal.bam: $(INTERM_DIR)/%_normal_merge.sh
	echo "`date`: started merging normal bams"
	qsub -V -sync yes -cwd -S /bin/sh -N MergeNormal \
	     -o $(LOG_DIR)/$(notdir $(basename $^)).stdout \
		 -e $(LOG_DIR)/$(notdir $(basename $^)).stderr \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^
	echo "`date`: finished merging normal bams"

$(INTERM_DIR)/%_tumour_merge.sh: $(NORM_DIR)/%
	echo "#!/bin/bash" > $@
	echo samtools merge $(patsubst $(INTERM_DIR)/%_tumour_merge.sh, $(INTERM_DIR)/%_tumour.bam, $@) $(wildcard $^/*) >> $@

$(INTERM_DIR)/%_tumour.bam: MEM := 1G
$(INTERM_DIR)/%_tumour.bam: $(INTERM_DIR)/%_tumour_merge.sh
	echo "`date`: started merging tumour bams"
	qsub -V -sync yes -cwd -S /bin/sh -N MergeTumour \
	     -o $(LOG_DIR)/$(notdir $(basename $^)).stdout \
		 -e $(LOG_DIR)/$(notdir $(basename $^)).stderr \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^
	echo "`date`: finished merging tumour bams"

# index bams
$(INTERM_DIR)/%_normal_index.sh: $(INTERM_DIR)/%_normal.bam
	echo "#!/bin/bash" > $@
	echo samtools index $^ >> $@

$(INTERM_DIR)/%_normal.bam.bai: MEM := 1G
$(INTERM_DIR)/%_normal.bam.bai: $(INTERM_DIR)/%_normal_index.sh
	echo "`date`: started indexing normal bams"
	qsub -V -sync yes -cwd -S /bin/sh -N IndexNormal \
	     -o $(LOG_DIR)/$(notdir $(basename $^)).stdout \
		 -e $(LOG_DIR)/$(notdir $(basename $^)).stderr \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^
	echo "`date`: finished indexing normal bams"

$(INTERM_DIR)/%_tumour_index.sh: $(INTERM_DIR)/%_tumour.bam
	echo "#!/bin/bash" > $@
	echo samtools index $^ >> $@

$(INTERM_DIR)/%_tumour.bam.bai: MEM := 1G
$(INTERM_DIR)/%_tumour.bam.bai: $(INTERM_DIR)/%_tumour_index.sh
	echo "`date`: started indexing tumour bams"
	qsub -V -sync yes -cwd -S /bin/sh -N IndexTumour \
	     -o $(LOG_DIR)/$(notdir $(basename $^)).stdout \
		 -e $(LOG_DIR)/$(notdir $(basename $^)).stderr \
		 -l h_vmem=$(MEM) -l mem_token=$(MEM) -l mem_free=$(MEM) $^
	echo "`date`: finished indexing tumour bams"

# make titan input
$(INTERM_DIR)/%_input.txt: $(INTERM_DIR)/%_tumour.bam $(INTERM_DIR)/%_normal.bam $(INTERM_DIR)/%_tumour.bam.bai $(INTERM_DIR)/%_normal.bam.bai
	echo -e "tumour\ttumour\t$(word 1, $^)\tnormal\tnormal\t$(word 2, $^)" > $@

# run titan
$(INTERM_DIR)/%: $(INTERM_DIR)/%_input.txt $(INTERM_DIR)/%_tumour.bam $(INTERM_DIR)/%_normal.bam
	echo "`date`: started TITAN"
	python TITANRunner/TitanRunner.py -i $(word 1, $^) \
                                      --project-name $(notdir $@) \
                                      --project-path $(INTERM_DIR) \
                                      -c config.cfg
	echo "`date`: finished TITAN"
	mv $@_R* $@
	rm $(word 2, $^)
	rm $(word 3, $^)

# pull out segments and parameters
$(TITAN_DIR)/%_segs.txt.bz2: $(INTERM_DIR)/%
	./parse_titan.py $^ $(TITAN_DIR)/$(notdir $^)

$(TITAN_DIR)/%_params.txt.bz2: $(INTERM_DIR)/%
	./parse_titan.py $^ $(TITAN_DIR)/$(notdir $^)

# compile into one file
segments.tsv: $(ALL_SEGS)
	echo -e "sample\tchrom\tstart\tend\tcopy.number" > segments.tsv
	for SEGS in $^; do \
		bzcat $$SEGS | tail -n +2 | cut -f 1-4,10 >> segments.tsv; \
	done

purity.tsv: $(ALL_PARAMS)
	for PARAMS in $^; do \
		SAMPLE=`basename $$PARAMS | cut -d '_' -f 1`; \
		CONTAM=`bzgrep contamination $$PARAMS | cut -f 2`; \
		PURITY=`echo 1-$$CONTAM | bc`; \
		echo -e "$$SAMPLE\t$$PURITY" >> purity.tsv; \
	done
