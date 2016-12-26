bin_PROGRAMS=east katool trecall
east_SOURCES=recursion.c traceback.c seq.c statistics.c east.c opts.c output.c matrix.c
katool_SOURCES=katool.c statistics.c
trecall_SOURCES=trecall.c recursion.c traceback.c seq.c statistics.c opts.c output.c matrix.c ambig.c

dist_pkgdata_DATA=BLOSUM

include_HEADERS=east_types.h

AM_CFLAGS=-DDATADIR=\"$(pkgdatadir)\"
