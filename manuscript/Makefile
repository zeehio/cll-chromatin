latex:
	scholdoc settings.yaml main.md --natbib --output=main.tex

# Need the `latexmk` program to automate LaTeX builds (included with texlive)
pdf: latex
	latexmk -pdf main.tex

# Scholdoc native PDF output doesn't call bibtex
autopdf:
	scholdoc settings.yaml main.md --citeproc --to=latex --output=main_generatedByScholdoc.pdf

doc: latex
	latex2rtf main.tex  # outputs main.rst
	pandoc -f latex -t docx -o main.docx main.tex  # does not see references and other stuff

clean:
	rm -f main.aux
	rm -f main.blg
	rm -f main.out
	rm -f main.tex
	rm -f main.log
	rm -f main.fls
	rm -f main.fdb_latexmk
	rm -f main.bbl
	rm -f main.pdf
	rm -f main.html
	rm -f main_generatedByScholdoc.pdf
	rm -f main.rst
	rm -f main.doc

.PHONY: html latex pdf autopdf clean
