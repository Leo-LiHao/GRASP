# Makefile for lgrind, a LaTeX prettyprinter

# $Id: makefile,v 1.4 1998/01/24 22:59:03 ballen Exp $

# You might want to change these
DESTDIR=`pwd`
DEFSFILE=`pwd`/lgrindefs
TEXINPUTS=`pwd`

SOURCES=lgrind.c regexp.c retest.c vgrindefs.c v2lg.c
TEXFILES=lgrind.doc lgrind.sty fancyhead.doc fancyhead.sty \
	 docstrip.tex docstrip.cmd
MANPAGES=lgrind.1 lgrindefs.5
EXAMPLES=doc-lgrind.lg
CMDS=lgrind
OBJS=lgrind.o regexp.o retest.o vgrindefs.o v2lg.o

CFLAGS=-g
CC=cc

all:	$(CMDS)

lgrind: lgrind.o vgrindefs.o regexp.o
	$(CC) ${CFLAGS} -o lgrind lgrind.o vgrindefs.o regexp.o

lgrind.o: lgrind.c
	$(CC) ${CFLAGS} -DDEFSFILE=\"${DEFSFILE}\" -c lgrind.c

v2lg:	v2lg.o
	$(CC) ${CFLAGS} -o v2lg v2lg.o

tar:
	tar cvf /tmp/lgrind.tar \
	   Makefile ${MANPAGES} lgrindefs ${SOURCES} ${TEXFILES} ${EXAMPLES}
	compress /tmp/lgrind.tar

shar:
	shar Makefile ${MANPAGES} lgrindefs ${SOURCES} ${TEXFILES} ${EXAMPLES}
		> lgrind.shar

retest: retest.o regexp.o
	$(CC) ${CFLAGS} -o retest retest.o regexp.o

install: all
	install -s lgrind ${DESTDIR}/lgrind
	install -s lgrindefs ${DEFSFILE}
	install -c lgrind.doc ${TEXINPUTS}/lgrind.doc
	install -c lgrind.sty ${TEXINPUTS}/lgrind.sty
	install -c fancyhead.doc ${TEXINPUTS}/fancyhead.doc
	install -c fancyhead.sty ${TEXINPUTS}/fancyhead.sty
	install -c lgrind.1 /cs/man/man1/lgrind.1
	install -c lgrindefs.5 /cs/man/man5/lgrindefs.5

clean:
	/bin/rm -f ${CMDS} ${OBJS} retest v2lg

