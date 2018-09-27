#!/bin/sh
pushd /home/ester/git/synthlisa
python setup.py install --prefix=/home/ester/anaconda2 $*
popd
