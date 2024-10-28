#!/bin/bash

DIR=build_debug

# [ -d "$DIR" ] && rm -rf "$DIR"
[ ! -d "$DIR" ] && mkdir "$DIR"
cd "$DIR"
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j
