include modules/Makefile.inc

LOGDIR ?= log/check_bam.$(NOW)
PHONY += check_bam

check_bam : $(foreach sample,$(SAMPLES),check_bam/$(sample).txt) check_bam/summary.txt

CHECK_BAM ?= $(wildcard $(foreach set,$(SAMPLES),check_bam/$(set).txt))

define check-bam
check_bam/%.txt : bam/%.bam
	$$(call RUN,-c -n 1 -s 2G -m 4G,"echo $$(*) > check_bam/$$(*).txt && \
									 echo ' ' >> check_bam/$$(*).txt && \
									 if [ -f $$(<) ]; then echo '1' >> check_bam/$$(*).txt; else echo '0' >> check_bam/$$(*).txt; fi")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call check-bam,$(sample))))
		
check_bam/summary.txt: $(wildcard check_bam/$(SAMPLES).txt)
	$(call RUN,-n 1 -s 2G -m 4G,"cat $(CHECK_BAM) > check_bam/summary.txt")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)