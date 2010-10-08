include_HEADERS = rings.h

bin_PROGRAMS = rings rings_elements_to_bodies rings_output_to_elements

lib_LTLIBRARIES = librings.la

librings_la_SOURCES = vectors.c body.c raw_average.c analytic_average.c advancer.c read_write.c

rings_SOURCES = rings.c
rings_LDADD = librings.la

rings_elements_to_bodies_SOURCES = rings_elements_to_bodies.c
rings_elements_to_bodies_LDADD = librings.la

rings_output_to_elements_SOURCES = rings_output_to_elements.c
rings_output_to_elements_LDADD = librings.la