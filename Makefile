#
# Makefile for LARCC framework
#

NAME = larcc
LANGUAGE = py
BIBFILE = $(NAME).bib

IDIR = src/tex/
ODIR = lib/
DOCTEX = doc/tex/
DOCPDF = doc/pdf/
TESTDIR = test/py/$(NAME)/
TESTFILES = TESTDIR/*.py

all: 
	echo building $(NAME)
	make pdf
	make clean
	open $(DOCPDF)$(NAME).pdf

exec:
	cp $(IDIR)macros.tex macros.tex
	cp $(IDIR)bib.bib $(BIBFILE)
	cp -R $(IDIR)images . 
	cp $(IDIR)$(NAME).tex $(NAME).w
	
	nuweb $(NAME).w

pdf: $(IDIR)$(NAME).tex
	make exec
	
	pdflatex $(NAME).tex
	nuweb $(NAME)
	bibtex $(NAME)
	
	pdflatex $(NAME).tex
	pdflatex $(NAME).tex

html:
	make pdf
	
	rm -dfR $(NAME)/*
	rm -dfR $(NAME)
	mkdir $(NAME)
	cp src/html/css.cfg $(NAME).cfg
	makeindex $(NAME).tex
	htlatex $(NAME).tex "$(NAME).cfg,TocLink,html,index=2,3"
	mv -fv images $(NAME).html $(NAME).css $(NAME)
	rm -fv $(NAME).* macros.tex
	mv -fv $(NAME)*.* $(NAME)
	if [ -d doc/html/$(NAME) ] ; then rm -R doc/html/$(NAME) ; fi
	mv $(NAME) doc/html/
	open doc/html/$(NAME)/$(NAME).html

tests: $(TESTFILES)
	@echo $^
clean:
	mv -fv $(NAME).tex $(NAME).bbl macros.tex $(DOCTEX)
	mv -fv $(NAME).pdf $(DOCPDF)
	mv -fv $(NAME).w $(ODIR)w
	if [ -d $(DOCTEX)images ] ; then rm -R $(DOCTEX)images ; fi
	mv -fv images $(DOCTEX)
	rm $(NAME).*
