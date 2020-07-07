#!/bin/bash

#export TAPIR_PREFIX=/efs/tools/opencilk/llvm/

#export PATH=$TAPIR_PREFIX/bin:$PATH

#export C_INCLUDE_PATH=$TAPIR_PREFIX/lib/clang/5.0.0/include/
#export CPLUS_INCLUDE_PATH=$TAPIR_PREFIX/lib/clang/5.0.0/include/
export CXX=clang++
#export LD_LIBRARY_PATH=x86_64:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$TAPIR_PREFIX/lib:$LD_LIBRARY_PATH
#export EXTRA_CFLAGS="-fcilkplus -Wall -Werror"

export N_TEMPORARY_BYTES=500000000
$@
