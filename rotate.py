from numpy import array, dot
from math import sin, cos

def euler_matrix(a,b,g):
    sa,sb,sg = sin(a), sin(b), sin(g)
    ca,cb,cg = cos(a), cos(b), cos(g)
    return array([ [-sa*cb*sg+ca*cg,    ca*cb*sg+sa*cg, sb*sg],
                   [-sa*cb*cg-ca*sg,    ca*cb*cg-sa*sg, sb*cg],
                   [          sa*sb,            -ca*sb,    cb]])

def rotate(M,x):
    return dot(M,x.T).T

def transform(M,x):
    return dot(M,x.T).T

def transform_list(M,x):
    return dot(M,array(x).T).T

def euler_rotate(a,b,g,x):
    return rotate(euler_matrix(a,b,g),x)
