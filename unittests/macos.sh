#!/bin/bash
set -e
for ccfile in unittests/src/*; do
   exefile="${ccfile%.*}".exe
   clang++ "$ccfile" -Iinclude -O3 -Wall -Wpedantic -Wextra -std=c++17 -o "$exefile"
   "$exefile"
done
