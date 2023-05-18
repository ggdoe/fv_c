#!/bin/sh

SRC=./src/*.c
out_name="a.out"

LIB="-lSDL2"

set -xe
gcc  $SRC -o $out_name $LIB


