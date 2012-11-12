#!/usr/bin/python

import sys
sys.path.insert( 0, "../scikit-learn" )
from sklearn import linear_model
from sklearn import datasets

if __name__ == '__main__':
    diabetes=datasets.load_diabetes()
    X=diabetes.data
    print X.shape
    y=diabetes.target
    print y.shape
    alphas, added_predictor, coefs = linear_model.lars_path( X, y, method='lasso', \
                                                   verbose=True, non_negative = True)

    print zip( alphas, added_predictor  )
    print coefs.shape

    print "\n\ncoefs:"
    print coefs.T
    print "\n\nalphas"
    print alphas
