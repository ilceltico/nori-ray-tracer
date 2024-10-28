#!/bin/bash

DIR=build

# [ -d "$DIR" ] && rm -rf "$DIR"
[ ! -d "$DIR" ] && mkdir "$DIR"
cd "$DIR"
cmake ..
make -j
