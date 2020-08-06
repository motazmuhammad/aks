import sympy
from sympy.abc import x
from sympy import rem, poly, trunc
from sympy import rootof, I
import math
import numpy as np

def powmod(p,n,r):
    ret=np.zeros(r,dtype=int);
    ret[0]=ret[1]=1;
    if n==0 :
        return ret;
    ret=np.zeros(r,dtype=int);
    for i in range(len(p)):
        ret[i%r]+=p[i];
    if n==1 :
        return ret;
    if(n%2==0):
        p=powmod(np.polymul(p,p),n//2,r);
        ret=np.zeros(r,dtype=int);
        for i in range(len(p)):
            ret[i%r]+=p[i%r];
        return ret;
    p=np.polymul(p,powmod(np.polymul(p,p),n//2,r));
    ret=np.zeros(r,dtype=int);
    for i in range(len(p)):
        ret[i%r]+=p[i%r];
    return ret;
     

p1=np.zeros(5,dtype=int);
p1[0]=p1[1]=1;
p=powmod(p1,5,5);
print(p);