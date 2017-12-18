#!/usr/bin/env bash
#
# Runs all unit tests included in PINTS.
#
# To run all tests:
#   $ ./run-test.sh
#
# To run a specific test:
#   $ test/test_logistic_model.py
#
# This file is part of PINTS.
#  Copyright (c) 2017, University of Oxford.
#  For licensing information, see the LICENSE file distributed with the PINTS
#  software package.
#
pyv="$(python -c 'import sys; print(sys.version_info[0])' 2>&1)"
if [ $pyv -eq 2 ]
then
    electrochemistry_dir=`pwd`
    export set PYTHONPATH=$PYTHONPATH:${electrochemistry_dir}
    ls
    mkdir -p build
    cd build
    cmake -DCMAKE_BUILD_TYPE=$PINTS_BUILD_TYPE ..
    make
    python2 -m unittest discover -v ../test
    exit_code=$?
    exit $exit_code
else
    echo "Electrochemistry tests require Python 2"
fi
