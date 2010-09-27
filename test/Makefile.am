check_PROGRAMS = vector-check body-test projection-test rotate-test \
alloc-body-from-elements-test

check_HEADERS = checks.h

TESTS = $(check_PROGRAMS)

vector_check_SOURCES = vector-check.c checks.c
vector_check_LDADD = ../src/librings.la

body_test_SOURCES = body-test.c checks.c
body_test_LDADD = ../src/librings.la

projection_test_SOURCES = checks.c projection-test.c
projection_test_LDADD = ../src/librings.la

rotate_test_SOURCES = checks.c rotate-test.c
rotate_test_LDADD = ../src/librings.la

alloc_body_from_elements_test_SOURCES = checks.c alloc-body-from-elements-test.c
alloc_body_from_elements_test_LDADD = ../src/librings.la