#!/bin/bash

rm -rf ct-build
mkdir ct-build
cd ct-build
cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_CLANG_TIDY:STRING="clang-tidy;-checks=*" ..
make
