check_PROGRAMS = vector-check body-test projection-test rotate-test \
init-body-from-elements-test \
raw-average-zero-adot-test raw-average-orthogonal-test raw-average-Ltot-test \
raw-average-rhs-earth-jupiter

check_HEADERS = checks.h
check_LTLIBRARIES = libcheck.la

libcheck_la_SOURCES = checks.c

TESTS = $(check_PROGRAMS)

vector_check_SOURCES = vector-check.c
vector_check_LDADD = ../src/librings.la libcheck.la

body_test_SOURCES = body-test.c
body_test_LDADD = ../src/librings.la libcheck.la

projection_test_SOURCES = projection-test.c
projection_test_LDADD = ../src/librings.la libcheck.la

rotate_test_SOURCES = rotate-test.c
rotate_test_LDADD = ../src/librings.la libcheck.la

init_body_from_elements_test_SOURCES = init-body-from-elements-test.c
init_body_from_elements_test_LDADD = ../src/librings.la libcheck.la

raw_average_zero_adot_test_SOURCES = raw-average-zero-adot-test.c
raw_average_zero_adot_test_LDADD = ../src/librings.la libcheck.la

raw_average_orthogonal_test_SOURCES = raw-average-orthogonal-test.c
raw_average_orthogonal_test_LDADD = ../src/librings.la libcheck.la

raw_average_Ltot_test_SOURCES = raw-average-Ltot-test.c
raw_average_Ltot_test_LDADD = ../src/librings.la libcheck.la

raw_average_rhs_earth_jupiter_SOURCES = raw-average-rhs-earth-jupiter.c
raw_average_rhs_earth_jupiter_LDADD = ../src/librings.la libcheck.la