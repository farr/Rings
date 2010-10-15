check_PROGRAMS = vector-check body-test projection-test rotate-test \
init-body-from-elements-test \
raw-average-zero-adot-test raw-average-orthogonal-test raw-average-Ltot-test \
raw-average-solar-system raw-average-A-L-mags-test \
raw-average-jupiter-saturn \
raw-average-yoram-elements \
analytic-average-test \
average-rhs-test \
quad-test

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

raw_average_solar_system_SOURCES = raw-average-solar-system.c
raw_average_solar_system_LDADD = ../src/librings.la libcheck.la

raw_average_A_L_mags_test_SOURCES = raw-average-A-L-mags-test.c
raw_average_A_L_mags_test_LDADD = ../src/librings.la libcheck.la

raw_average_jupiter_saturn_SOURCES = raw-average-jupiter-saturn.c
raw_average_jupiter_saturn_LDADD = ../src/librings.la libcheck.la

raw_average_yoram_elements_SOURCES = raw-average-yoram-elements.c
raw_average_yoram_elements_LDADD = ../src/librings.la libcheck.la

analytic_average_test_SOURCES = analytic-average-test.c
analytic_average_test_LDADD = ../src/librings.la libcheck.la

average_rhs_test_SOURCES = average-rhs-test.c
average_rhs_test_LDADD = ../src/librings.la libcheck.la

quad_test_SOURCES = quad-test.c
quad_test_LDADD = ../src/librings.la libcheck.la