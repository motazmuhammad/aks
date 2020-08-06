# the aim of this project is to create a python code for primality testing based on this paper
#  https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf

from sympy import perfect_power
from sympy.ntheory import isprime
import math
def checkRCandidate(r,n):
    for i in range(1,4*n.bit_length()+2):
        if math.pow(n,i,r)==1 :
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
    return 0;

def isPerfectPower(n):
    if perfect_power(n)==False :
        return False;
    return True;

def aks(n):
    if n==1 :
         return False;# step 0

    if( isPerfectPower(n)) :# step 1
        return False;

    r=findR(n);
    if r==n :# step 2
        return True;
        
    end=2*ceil(sqrt(r))*(n.bit_length()+1);
    for a in range(1,end) :
        if math.gcd(a,n)>1 :
            return False
        if checkWitness(a) :
            return False
    return True;


n = int(input());
print(aks(n));