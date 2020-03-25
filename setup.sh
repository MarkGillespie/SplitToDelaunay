#!/bin/bash

mkdir deps
cd deps
git submodule add https://github.com/MarkGillespie/polyscope.git
git submodule add https://github.com/MarkGillespie/geometry-central.git
git submodule add https://github.com/google/googletest
git submodule update --init --recursive
cd ..
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release

cd ../src
mkdir build
cd build
ln -s ../../build/compile_commands.json compile_commands.json
