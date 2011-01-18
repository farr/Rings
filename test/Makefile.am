check_PROGRAMS = vector-check body-test projection-test rotate-test \
init-body-from-elements-test \
raw-average-zero-adot-test raw-average-orthogonal-test raw-average-Ltot-test \
raw-average-solar-system raw-average-A-L-mags-test \
raw-average-jupiter-saturn \
raw-average-yoram-elements \
analytic-average-test \
analytic-average-angular-momentum \
average-rhs-test \
quad-test \
extreme-eccentricity-test \
analytic-average-A-L \
advancer-test

check_HEADERS = checks.h
check_LTLIBRARIES = libcheck.la

libcheck_la_SOURCES = checks.c

TESTS = $(check_PROGRAMS)

LDADD = ../src/librings.la libcheck.la