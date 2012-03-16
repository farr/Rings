#!/bin/sh

set -e

autoreconf || (libtoolize; autoreconf; automake --add-missing; autoreconf)