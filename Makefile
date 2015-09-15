.DEFAULT_GOAL := all

requirements:
	pip install https://github.com/epigen/pipelines.git

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
	cd manuscript
	make pdf

all: requirements preprocess external_files analysis manuscript

.PHONY: requirements preprocess external_files analysis manuscript all
