#!/bin/sh

SRC="./src/*.c ./src/sim/*.c"
out_name="a.out"

CFLAGS="-Wall"
LFLAGS="-flto"
LIB="-lm -lSDL2 -lx264"

set -xe
gcc  $CFLAGS $LFLAGS -o $out_name $SRC $LIB


