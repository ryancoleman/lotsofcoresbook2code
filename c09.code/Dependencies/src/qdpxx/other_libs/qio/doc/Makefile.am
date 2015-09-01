# For making Postscript and html documentation

docdir=@prefix@/doc

DOCSRCDIR=@top_srcdir@/doc

LATEX = latex

VER     = 2p3
HTMLDIR = html
PS      = qio_${VER}.ps
SRC     = qio_${VER}.tex
DVI     = qio_${VER}.dvi

.PHONY: ps html
ps: ${PS}
html: ${HTMLDIR}


${PS}: ${DOCSRCDIR}/${SRC}
	${LATEX} ${DOCSRCDIR}/${SRC}
	${LATEX} ${DOCSRCDIR}/${SRC}
	dvips ${DVI} -o

${HTMLDIR}: ${DOCSRCDIR}/${SRC}
	mkdir -p $@
	ln -sf ${DOCSRCDIR}/${SRC} ${SRC}
	latex2html -math -dir ${HTMLDIR} ${SRC}

doc: ps html

doc-install: ${PS} ${HTMLDIR}
	mkdir -p ${docdir}
	cp -r ${PS} ${HTMLDIR} ${docdir}

doc-uninstall:
	-/bin/rm -rf ${docdir}

# Removes LaTeX debris
doc-clean:
	-/bin/rm -f *.aux *.log *.err *.dvi

# Removes all but the .tex file and Makefile
doc-distclean: clean
	-/bin/rm -f *.ps
	-/bin/rm -rf ${HTMLDIR}

# also clean on regular make clean
clean-local: doc-clean

EXTRA_DIST = ${SRC}

.PHONY: doc doc-install doc-uninstall doc-clean doc-distclean


