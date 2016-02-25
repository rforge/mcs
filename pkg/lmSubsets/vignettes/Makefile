all:
	pdflatex --shell-escape lmSubsets.tex
	bibtex lmSubsets.aux
	pdflatex --shell-escape lmSubsets.tex
	pdflatex --shell-escape lmSubsets.tex

clean:
	rm -f lmSubsets.tex
	rm -f lmSubsets.aux
	rm -f lmSubsets.bbl
	rm -f lmSubsets.blg
	rm -f lmSubsets.log
	rm -f lmSubsets.out
	rm -f lmSubsets-*.pdf
