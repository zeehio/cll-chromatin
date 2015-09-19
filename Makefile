.DEFAULT_GOAL := all

requirements:
	pip install https://github.com/epigen/pipelines.git

	# python stuff:
	pip install numpy pandas scipy pysam pybedtools matplotlib seaborn parmap sklearn statsmodels

	# R stuff:
	# lola seq2pathway

	# other tools:
	# bwa macs samtools sambamba bedtools meme

preprocess: requirements
	pipelines cll-project metadata/sequencing_sample_annotation.csv

external_files:
	python src/prepare_external_files.py

analysis: preprocess external_files
	python src/cohort_description.py
	python src/seq_stats.py
	python src/analysis.py
	python src/grn_analysis.py
	python src/call_variants.py

manuscript:
	cd manuscript && $(MAKE) pdf

all: requirements preprocess external_files analysis manuscript

clean:
	cd manuscript && $(MAKE) clean

.PHONY: requirements preprocess external_files analysis manuscript all
