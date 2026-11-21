#!/usr/bin/env python3
"""W1: vacuum GR axial QNM, not a massive-mode or sourced-response solver.

Python standard library only. Units r_g=GM=1; exp(-i omega t).
Run from any directory; stdout is JSON unless --output is supplied.
See notes/WEYL_RINGDOWN_AXIAL_DERIVATION.md for the derived recurrence.
"""
import argparse
import cmath
import json
import platform
from pathlib import Path


def recurrence(n, w, ell=2):
    return ((n + 1) * (n + 1 - 4j * w),
            -2*n*n + (-2 + 16j*w)*n + 32*w*w + 8j*w + 3-ell*(ell+1),
            n*n - 8j*w*n - 16*w*w - 4)


def ratios(w, depth, tail=0j):
    """R[n]=a[n+1]/a[n], starting with a[depth+1]/a[depth]=tail."""
    result = [0j] * (depth + 1)
    result[depth] = tail
    for n in range(depth-1, -1, -1):
        a, b, c = recurrence(n+1, w)
        result[n] = -c / (b + a*result[n+1])
    return result


def continued_fraction(w, depth, tail=0j):
    a, b, _ = recurrence(0, w)
    return b + a*ratios(w, depth, tail)[0]


def solve(depth, tail=0j, guesses=(0.37-0.09j, 0.38-0.08j)):
    x, y = guesses
    for iteration in range(80):
        fx, fy = continued_fraction(x, depth, tail), continued_fraction(y, depth, tail)
        if abs(fy-fx) < 1e-30:
            raise ArithmeticError('Secant denominator vanished')
        z = y - fy*(y-x)/(fy-fx)
        if abs(z-y) < 2e-14 and abs(continued_fraction(z, depth, tail)) < 1e-11:
            if not (0.3 < z.real < 0.45 and -0.15 < z.imag < -0.03):
                raise ArithmeticError('Root left the requested fundamental-mode neighborhood')
            return z, iteration+1
        x, y = y, z
    raise ArithmeticError('Root finder did not converge')


def coefficients(w, depth):
    values = [1+0j]
    for ratio in ratios(w, depth)[:-1]:
        values.append(values[-1]*ratio)
    return values


def radial(w, coeff, r):
    """Return psi, d_r psi, d_r^2 psi by analytic series differentiation."""
    x = 1-2/r
    xp, xpp = 2/r**2, -4/r**3
    a = sum(c*x**n for n,c in enumerate(coeff))
    ax = sum(n*coeff[n]*x**(n-1) for n in range(1,len(coeff)))
    axx = sum(n*(n-1)*coeff[n]*x**(n-2) for n in range(2,len(coeff)))
    p = cmath.exp(1j*w*r + 2j*w*cmath.log(r) - 2j*w*cmath.log(x))
    lp = 1j*w + 2j*w/r - 2j*w*xp/x
    lpp = -2j*w/r**2 - 2j*w*(xpp/x-(xp/x)**2)
    return (p*a, p*(lp*a+xp*ax),
            p*((lp*lp+lpp)*a+(2*lp*xp+xpp)*ax+xp*xp*axx))


def metric(w, coeff, r):
    psi, dpsi, _ = radial(w, coeff, r)
    f = 1-2/r
    # H_tA=h0 S_A, H_rA=h1 S_A in RW gauge, S_A fixed in note.
    return 1j*f*(psi+r*dpsi)/w, r*psi/f


def diagnostics(w, depth):
    coeff = coefficients(w, depth)
    ode = []
    for r in (2.1, 3., 6., 10., 20.):
        psi, dp, ddp = radial(w, coeff, r)
        f, fp = 1-2/r, 2/r**2
        terms = (f*f*ddp, f*fp*dp, (w*w-f*(6/r**2-6/r**3))*psi)
        ode.append({'r_over_rg':r, 'normalized_RW_residual':abs(sum(terms))/sum(map(abs,terms))})
    constraints=[]
    for step in (0.02,0.01,0.005,0.0025):
        worst = 0.
        for r in (3.,6.,10.):
            h0,_ = metric(w,coeff,r)
            def fh1(x):
                return (1-2/x)*metric(w,coeff,x)[1]
            derivative = (fh1(r+step)-fh1(r-step))/(2*step)
            first,second=1j*w*h0,(1-2/r)*derivative
            worst=max(worst,abs(first+second)/(abs(first)+abs(second)))
        constraints.append({'centered_difference_step_over_rg':step,
                            'max_normalized_angular_Einstein_constraint':worst})
    return {'radial_equation':ode,'metric_constraint':constraints}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    # Berti s2l2.dat first row: 2M Re(omega), positive 2M damping magnitude.
    # Used only for comparison, never inside root finder or boundary data.
    reference=complex(0.7473433688360838/2,-0.1779246313778714/2)
    runs=[]
    for depth in (12,25,50,100,200,400):
        w,iterations=solve(depth)
        runs.append({'depth':depth,'omega_rg':[w.real,w.imag],
                     'abs_error_vs_reference':abs(w-reference),
                     'continued_fraction_residual':abs(continued_fraction(w,depth)),
                     'iterations':iterations,'diagnostics':diagnostics(w,depth)})
    tails=[]
    for tail in (0j,0.5+0j,1+0j):
        w,_=solve(200,tail)
        tails.append({'tail':[tail.real,tail.imag],'omega_rg':[w.real,w.imag],
                      'abs_error_vs_reference':abs(w-reference)})
    seeds=[]
    for pair in ((.36-.1j,.39-.07j),(.35-.08j,.4-.1j)):
        w,_=solve(200,guesses=pair)
        seeds.append({'guesses':[[z.real,z.imag] for z in pair],
                      'omega_rg':[w.real,w.imag]})
    assert max(r['abs_error_vs_reference'] for r in runs[-3:]) < 1e-11
    assert max(r['abs_error_vs_reference'] for r in tails) < 1e-11
    assert max(p['normalized_RW_residual'] for p in runs[-1]['diagnostics']['radial_equation']) < 1e-10
    c=runs[-1]['diagnostics']['metric_constraint']
    assert c[-1]['max_normalized_angular_Einstein_constraint'] < c[0]['max_normalized_angular_Einstein_constraint']/30
    result={'gate':'W1','model':'vacuum massless GR axial l=2 fundamental',
            'python':platform.python_version(),'units':'r_g=1; exp(-i omega t)',
            'method':'derived scalar continued fraction; analytic horizon/infinity factors',
            'reference_url':'https://pages.jh.edu/eberti2/ringdown/s2l2.dat',
            'reference_omega_rg':[reference.real,reference.imag],
            'precision_note':'double-precision numerical agreement, not certified error bounds',
            'runs':runs,'tail_variations':tails,'seed_variations':seeds,
            'checks_passed':True,'massive_roots_computed':False}
    output=json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(output)
    else:
        print(output,end='')


if __name__=='__main__':
    main()
