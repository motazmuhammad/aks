# the aim of this project is to create a python code for primality testing based on this paper
#  https://www.cse.iitk.ac.in/users/manindra/algebra/primality_v6.pdf

from sympy import perfect_power
from sympy.ntheory import isprime
import math
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
    if len(p1) != len(p2):
        return False;
    result = [0]*r;
    for i in range(r):
        for j in range(r):
            result[(i+j)%r]+=p1[i]*p2[j];
            result[(i+j)%r]%=n;
    return result;

def polypow(a,n,r,m): # calculates (x+a)**n %(n,x**r-1)
    x = [0]*r;
    result=[0]*r;
    x[0]=a;    
    x[1]=1;
    result[0]=1;

    

    while n > 0 :
        if n%2== 1 :
            result=polyMul(result,x,m);
        n = n//2;  
        x = polyMul(x,x,m);
    return result; 

def checkWitness(a,n,r): 
    LHS=polypow(a,n,r,n);
    RHS = [0]*r;
    RHS[0]=a;
    RHS[n%r]=1;
    for i in range(r):
        if LHS[i]!=RHS[i] :
             return True;
    return False;


def aks(n):
    if n==1 :
         return False;# step 0

    if( isPerfectPower(n)) :# step 1
        return False;

    r=findR(n);

    if r==False :
        return False
    if r==n :# step 2
        return True;
    
    end=min(2*math.ceil(math.sqrt(r))*(n.bit_length()+1),n);
    for a in range(1,end) :
        if math.gcd(a,n)>1 :
            return False
        if checkWitness(a,n,r) :
            return False
    return True;



for i in range(31,10000):
    withAks=aks(i);
    withSympy=isprime(i);
    if i%100==0 :
        print(i)
    if withAks!=withSympy:
        print('come on something is wrong')