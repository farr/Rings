AM_CFLAGS=$(GSL_CFLAGS)
AM_LDFLAGS=$(GSL_LIBS)

include_HEADERS = rings.h

bin_PROGRAMS = rings

lib_LTLIBRARIES = librings.la

librings_la_SOURCES = vectors.c body.c raw_average.c		\
analytic_average.c advancer.c read_write.c tides.c

rings_SOURCES = rings.c
rings_LDADD = librings.la