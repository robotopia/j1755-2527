TARGETS = main.pdf
#IMAGES = 

all: $(TARGETS)

%.pdf: %.tex
	pdflatex $<
	bibtex $(basename $<)
	pdflatex $<
	pdflatex $<
