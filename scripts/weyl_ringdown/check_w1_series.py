#!/usr/bin/env python3
"""Test local series against the differential operator, not against stored coefficients."""
import json
from axial_system import horizon_series, infinity_series, massive_W, mv, add, scale


def norm(v):
    return max(map(abs,v))


def polynomial(coeff,x,inverse=False):
    values=[]
    for derivative in range(3):
        value=(0j,0j)
        for n,c in enumerate(coeff):
            power=-n if inverse else n
            factor=1
            for j in range(derivative): factor*=power-j
            if factor: value=add(value,scale(factor*x**(power-derivative),c))
        values.append(value)
    return values


def main():
    w,mu,ell=.5-.08j,.2,2
    rows=[]
    for order in (2,4,6):
        h=horizon_series(w,mu,ell,(1,.3),order)
        k,chi,b=infinity_series(w,mu,ell,(1,.3),order)
        horizons=[]; infinities=[]
        for u in (.1,.05):
            r=2+u; f=1-2/r; fp=2/r**2
            v,d,dd=polynomial(h,u)
            residual=add(add(scale(f,dd),scale(fp-2j*w,d)),scale(-1,mv(massive_W(r,ell,mu),v)))
            horizons.append({'r_minus_2':u,'factored_operator_residual':norm(residual)})
        for r in (20.,40.):
            f=1-2/r; fp=2/r**2
            v,d,dd=polynomial(b,r,True)
            lp=1j*k+chi/r
            lpp=-chi/r**2
            residual=add(scale(f*f,dd),scale(2*f*f*lp+f*fp,d))
            residual=add(residual,scale(f*f*(lp*lp+lpp)+f*fp*lp+w*w,v))
            residual=add(residual,scale(-f,mv(massive_W(r,ell,mu),v)))
            infinities.append({'r':r,'factored_operator_residual':norm(residual)})
        rows.append({'order':order,'horizon':horizons,'infinity':infinities})
    # These are local-series convergence checks at a generic frequency, not modes.
    assert rows[-1]['horizon'][0]['factored_operator_residual'] < rows[0]['horizon'][0]['factored_operator_residual']/1000
    assert rows[-1]['infinity'][0]['factored_operator_residual'] < rows[0]['infinity'][0]['factored_operator_residual']/100
    assert all(row['horizon'][1]['factored_operator_residual'] < row['horizon'][0]['factored_operator_residual'] for row in rows)
    assert all(row['infinity'][1]['factored_operator_residual'] < row['infinity'][0]['factored_operator_residual'] for row in rows)
    print(json.dumps({'test':'local massive axial boundary series, not a QNM',
                      'omega':[w.real,w.imag],'mu':mu,'ell':ell,'runs':rows,
                      'checks_passed':True},indent=2))


if __name__=='__main__': main()
