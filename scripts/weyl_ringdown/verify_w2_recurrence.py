#!/usr/bin/env python3
"""Exact check of the massive three-term recurrence; requires sympy==1.14.0."""
import json
import sympy as s

x,w,k,a,b,n=s.symbols('x omega k a b n',nonzero=True)
r=2/(1-x);xd=(1-x)**2/2
chi=s.I*(w*w+k*k)/k
lp=s.I*k+chi/r-2*s.I*w*xd/x
u=s.Function('u')(x)
ur=xd*s.diff(u,x)+lp*u
urr=xd*s.diff(ur,x)+lp*ur
operator=x*x*urr+x*(2/r**2)*ur+(w*w-x*(w*w-k*k+a/r**2+b/r**3))*u
P=x*(1-x)**2
b0=1-4*s.I*w
b1=6*s.I*k+8*s.I*w-4+2*s.I*w*w/k
b2=-2*s.I*k-4*s.I*w+3-2*s.I*w*w/k
c0=4*k*k+12*k*w+3*s.I*k+12*w*w+4*s.I*w+4*w**3/k+s.I*w*w/k
c1=-k*k-4*k*w-2*s.I*k-6*w*w-4*s.I*w-4*w**3/k-2*s.I*w*w/k-w**4/k**2
expected=P*s.diff(u,x,2)+(b0+b1*x+b2*x*x)*s.diff(u,x)+(c0-a-b/2+(c1+b/2)*x)*u
assert s.simplify(operator/(x*(1-x)**2/4)-expected)==0
# The scalar differential pieces multiply identity matrices; a,b may therefore
# be replaced by the constant, noncommuting W1 potential matrices A,B without
# multiplying them together. No assumption of simultaneous diagonalization.
alpha=(n+1)*(n+b0)
beta=-2*n*(n-1)+b1*n+c0-a-b/2
gamma=(n-1)*(n-2)+b2*(n-1)+c1+b/2
for index in range(1,5):
    terms=[]
    for power in (index+1,index,index-1):
        v=x**power
        transformed=s.expand(P*s.diff(v,x,2)+(b0+b1*x+b2*x*x)*s.diff(v,x)+(c0-a-b/2+(c1+b/2)*x)*v)
        terms.append(transformed.coeff(x,index))
    assert all(s.simplify(v-formula.subs(n,index))==0 for v,formula in zip(terms,(alpha,beta,gamma)))
print(json.dumps({'sympy':s.__version__,'checks_passed':True,
                  'operator_transformation':'exact symbolic identity for arbitrary a,b,omega,k',
                  'coefficient_extraction_indices':[1,2,3,4],
                  'mu_squared_substitution':'omega^2-k^2',
                  'matrix_extension':'scalar differential terms times identity; potential linear in A and B'},indent=2))
