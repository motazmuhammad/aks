# the aim of this project is to create a python code for primality testing based on this paper
#  https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf

from sympy import perfect_power
from sympy.ntheory import isprime
import math
import mpmath
import sympy
from mpmath import mp 

import time
import copy
#from sympy import fft, ifft
from sympy import convolution
from sympy import re
from sympy.discrete.convolutions import convolution_fft
from sympy.discrete.convolutions import convolution_ntt
import copy
import cProfile
import re
from mpmath import sin, cos
# f
mp.dps = 100
def fft( b,  invert=False) :
    a=[mp.mpc(0)]*len(b);
    a[0:len(b)]=b;
    n=len(a);
    if (n == 1):
        return a;
    a0=[mp.mpc(0)]*(n//2);
    a1=[mp.mpc(0)]*(n//2);
    for i in range(n//2) :
        a0[i] = mp.mpc(a[2*i]);
        a1[i] = mp.mpc(a[2*i+1]);
    a0=fft(a0, invert);
    a1=fft(a1, invert);
    factor=[-1,1];
    ang = (mp.mpf(2) * mp.pi / n) * factor[int(invert)];
    w=mp.mpc(1);
    wn=mp.exp(mp.mpc(real=0,imag=ang));
    
    for i in range(n//2) :
        a[i] = a0[i] + w * a1[i];
        a[i + n//2] = a0[i] - w * a1[i];
        if (invert) :
            a[i] /= mp.mpc(2);
            a[i + n//2] /= mp.mpc(2);
        w *= wn;

    
    return a;
    

def checkRCandidate(r,n):
    for i in range(1,4*n.bit_length()+2):
        if pow(n,i,r)==1 :
             return False;
    return True;
def findR(n):# step 2
    r=2;
    while r<n :
        if n%r==0 :
            return False;
        if isprime(r):
            if checkRCandidate(r,n) :
                return r;
        r=r+1;
    return r;

def isPerfectPower(n):
    if perfect_power(n)==False :
        return False;
    return True;

def polyMul(p1,p2,n):
    r=len(p1);
    result = [0]*r;
    for i in range(r):
        for j in range(r):
            result[(i+j)%r]+=p1[i]*p2[j];
            result[(i+j)%r]%=n;
    return result;

# def polyMulFast(p1,p2,n):
#     p1F=sympy.discrete.transforms.fft( p1)
#     p2F=sympy.discrete.transforms.fft( p2,15);
#     result=[0]*len(p1F);
#     for i in range(len(p1F)):
#         result[i]=p1F[i]*p2F[i];
#     result=sympy.discrete.transforms.ifft( result,15);
#     res=[0]*len(p1);
#     for i in range(len(result)):
#         res[i%len(p1)]+=int(round(sympy.re(result[i])));
#         res[i%len(p1)]%=n;
#     return res;

# def polyMulFast(p1,p2,n):
#     res=convolution_ntt(p1,p2,prime=5*2**75 + 1);
#     result=[0]*len(p1);
#     for i in range(len(res)):
#         result[i%len(p1)]+=int(res[i]);
#         result[i%len(p1)]%=n;
#     return result;

def polyMulFast(p1,p2,n):
    r=len(p1);
    m=2**(r.bit_length());
    p1F=[0]*m;
    p2F=[0]*m;
    p1F[0:len(p1)]=[mp.mpc(x) for x in p1];
    p2F[0:len(p1)]=[mp.mpc(x) for x in p2];
    p1F=fft(p1F);
    p2F=fft(p2F);
    result=[0]*m;
    for i in range(m):
        result[i]=p1F[i]*p2F[i];
    result=fft(result,True);
    res=[0]*r;
    for i in range(len(result)):
        res[i%r]+=int(round(result[i].real));
       # res[i%r]%=n;
    #res[i%len(p1)]%=n;
    return res;

# def fft(vector<cd> & a, bool invert) {
#     int n = a.size();
#     int lg_n = 0;
#     while ((1 << lg_n) < n)
#         lg_n++;

#     for (int i = 0; i < n; i++) {
#         if (i < reverse(i, lg_n))
#             swap(a[i], a[reverse(i, lg_n)]);
#     }

#     for (int len = 2; len <= n; len <<= 1) {
#         double ang = 2 * PI / len * (invert ? -1 : 1);
#         cd wlen(cos(ang), sin(ang));
#         for (int i = 0; i < n; i += len) {
#             cd w(1);
#             for (int j = 0; j < len / 2; j++) {
#                 cd u = a[i+j], v = a[i+j+len/2] * w;
#                 a[i+j] = u + v;
#                 a[i+j+len/2] = u - v;
#                 w *= wlen;
#             }
#         }
#     }

#     if (invert) {
#         for (cd & x : a)
#             x /= n;
#     }
# }



def polypow(a,n,r,m,fast=False): # calculates (x+a)**n %(n,x**r-1)
    x = [0]*r;
    result=[0]*r;
    x[0]=a;    
    x[1]=1;
    result[0]=1;
    while n > 0 :
        if n%2== 1 :
            if fast:
                result=polyMulFast(result,x,m);
            else :
                result=polyMul(result,x,m);
        n = n//2;
        if fast:  
            x = polyMulFast(x,x,m);
        else:
            x = polyMul(x,x,m);
    return result; 

def checkWitness(a,n,r,fast=False): 
    LHS=polypow(a,n,r,n,fast);
    RHS = [0]*r;
    RHS[0]=a;
    RHS[n%r]=1;
    for i in range(r):
        if LHS[i]!=RHS[i] :
             return True;
    return False;


def aks(n,fast=False):
    if n==1 :
         return False;# step 0

    if( isPerfectPower(n)) :# step 1
        return False;

    r=findR(n);

    if r==False :
        return False
    if r==n :# step 2
        return True;
    
    end=min(math.ceil(math.sqrt(r))*(n.bit_length()+1),n);
    for a in range(1,end) :
        if math.gcd(a,n)>1 :
            return False
        if checkWitness(a,n,r,fast) :
            return False
    return True;


n=1;
# while n>=1:
#     n =int(input());
#     a=1;
#     r=int(input());
#     print(polypow(a,n,r,100,True));
end=100
start_time = time.time()
for i in range(1,end):
    withAks=aks(i);
    withSympy=isprime(i);
    if i%100==0 :
        print(i)
    if withAks!=withSympy:
        print('come on something is wrong slow')
print("--- %s seconds ---" % (time.time() - start_time))

# start_time = time.time()
# for i in range(37,end):
#     withAks=aks(i,True);
#     withSympy=isprime(i);
#     if i%100==0 :
#         print(i)
#     if withAks!=withSympy:
#         print('come on something is wrong Fast')
# print("Fast--- %s seconds ---" % (time.time() - start_time))