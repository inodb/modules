include modules/Makefile.inc

LOGDIR ?= log/hotspot_summary.$(NOW)
PHONY += hotspot summary summary/tsv

HOTSPOT ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),summary/tsv/$(set).tsv))

hotspot_summary : summary/hotspot_summary.xlsx

summary/tsv/hotspot_summary.tsv : $(wildcard hotspot/$(SAMPLE_PAIRS).txt)
	$(call RUN,-n 1 -s 4G -m 4G,"$(RSCRIPT) modules/summary/hotspotsummary.R --in_file '$(HOTSPOT)' --out_file summary/tsv/hotspot_summary.tsv")
		
summary/hotspot_summary.xlsx : summary/tsv/hotspot_summary.tsv
	$(INIT) python modules/summary/hotspot_summary_excel.py

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
