#!/bin/sh

set -e

mkdir -p m4
libtoolize
automake --add-missing
autoreconf