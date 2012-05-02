#!/bin/bash
set -o errexit

git clone git://github.com/scikit-learn/scikit-learn.git
patch -p1 < ./setup/nonnegative_lars_scikit.patch
cd scikit-learn/
python setup.py build_ext --inplace

# this tests that the installation worked
cd "../$(dirname "$0")"
python least_angle_test.py
