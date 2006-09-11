all:
	pdflatex nem.tex
	pdflatex nem.tex
	cp nem.pdf nem.tmp
	cp markowetz-thesis-2006.pdf markowetz-thesis-2006.tmp
	-rm *.pdf *.eps *.bak *.log *.dvi *.aux 
	cp nem.tmp nem.pdf
	cp markowetz-thesis-2006.tmp markowetz-thesis-2006.pdf
