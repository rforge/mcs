all:
	pdflatex --shell-escape mcsSubset.tex
	bibtex mcsSubset.aux
	pdflatex --shell-escape mcsSubset.tex
	pdflatex --shell-escape mcsSubset.tex

clean:
	rm -f mcsSubset.tex
	rm -f mcsSubset.aux
	rm -f mcsSubset.bbl
	rm -f mcsSubset.blg
	rm -f mcsSubset.log
	rm -f mcsSubset.out
	rm -f mcsSubset-*.pdf
