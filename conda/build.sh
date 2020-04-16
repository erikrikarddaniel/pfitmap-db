#!/bin/sh

[ ! -d $PREFIX/bin ] && mkdir -p $PREFIX/bin
cp $SRC_DIR/src/R/*.r $PREFIX/bin
