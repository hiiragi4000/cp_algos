#!/bin/bash
set -e
for ccfile in unittests/src/*; do
   exefile="${ccfile%.*}".exe
   g++ "$ccfile" -Iinclude -O3 -Wall -Wpedantic -Wextra -std=c++17 -o "$exefile"
   "$exefile"
done
