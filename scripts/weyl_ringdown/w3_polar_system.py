"""Coordinate-derived constrained polar operator; GM=1, exp(-i omega t).

Requires SymPy. The fields are (a,b,c,d,e,K,G), with q_tt=F a Y,
q_tr=b Y, q_rr=c/F Y, q_tA=d D_A Y, q_rA=e D_A Y, and
q_AB=r^2(K Omega_AB Y + 2G Y_AB). No gauge is imposed on q.
"""
import argparse
import json
from functools import lru_cache
from pathlib import Path
import sympy as s
from sympy.core.function import AppliedUndef
from verify_w1_algebra import r,theta,w,mu,coords,f,g,gi,gamma,riemann,harmonic_diff


def trig_coefficient(expr):
    # Exact rational reduction on sin^2(theta)+cos^2(theta)=1 avoids
    # heuristic trigsimp leaving spurious angular dependence in projections.
    u,v=s.symbols('_sin _cos')
    expr=s.expand_trig(expr).subs({s.tan(theta):s.sin(theta)/s.cos(theta),
                                  s.cot(theta):s.cos(theta)/s.sin(theta)})
    expr=expr.subs({s.sin(theta):u,s.cos(theta):v})
    num,den=s.fraction(s.cancel(expr))
    num=s.rem(num,u*u+v*v-1,u);den=s.rem(den,u*u+v*v-1,u)
    return s.factor(s.cancel(num/den)).subs({u:s.sin(theta),v:s.cos(theta)})


def clean(expr):
    # Factor coefficient by coefficient, never a giant polynomial in field jets.
    jets=sorted(expr.atoms(s.Derivative)|expr.atoms(AppliedUndef),key=str)
    if not jets:return trig_coefficient(expr)
    matrix,rhs=s.linear_eq_to_matrix([s.expand(expr)],jets)
    return sum(trig_coefficient(c)*j for c,j in zip(matrix,jets))-rhs[0]


def derive(ell):
    names=['a','b','c','K'] if ell==0 else ['a','b','c','d','e','K'] if ell==1 else ['a','b','c','d','e','K','G']
    fields={n:s.Function(n)(r) for n in names}
    a,b,c,K=[fields[n] for n in ('a','b','c','K')]
    d=fields.get('d',s.S.Zero);e=fields.get('e',s.S.Zero);G=fields.get('G',s.S.Zero)
    Y=s.legendre(ell,s.cos(theta));Yp=s.diff(Y,theta)
    W=s.diff(Y,theta,2)-s.cot(theta)*Yp
    q=s.zeros(4)
    q[0,0]=f*a*Y;q[0,1]=q[1,0]=b*Y;q[1,1]=c*Y/f
    q[0,2]=q[2,0]=d*Yp;q[1,2]=q[2,1]=e*Yp
    q[2,2]=r*r*(K*Y+G*W);q[3,3]=r*r*s.sin(theta)**2*(K*Y-G*W)

    @lru_cache(None)
    def nabla(i,j,k):
        return clean(harmonic_diff(q[j,k],i)-sum(gamma(z,i,j)*q[z,k]+gamma(z,i,k)*q[j,z] for z in range(4)))

    def op(j,k):
        box=0
        for i in range(4):
            value=harmonic_diff(nabla(i,j,k),i)-sum(gamma(z,i,i)*nabla(z,j,k)+gamma(z,i,j)*nabla(i,z,k)+gamma(z,i,k)*nabla(i,j,z) for z in range(4))
            box+=gi[i,i]*value
        curvature=sum(2*g[j,j]*riemann(j,i,k,z)*gi[i,i]*gi[z,z]*q[i,z] for i in range(4) for z in range(4) if q[i,z]!=0)
        return clean(box+curvature-mu*mu*q[j,k])

    trace=clean(sum(gi[i,i]*q[i,i] for i in range(4))/Y)
    divergence={}
    for j in ([0,1] if ell==0 else [0,1,2]):
        divergence[str(j)]=clean(sum(gi[i,i]*nabla(i,i,j) for i in range(4))/(Yp if j==2 else Y))
    eq={'tt':clean(op(0,0)/Y),'tr':clean(op(0,1)/Y),'rr':clean(op(1,1)/Y)}
    if ell>=1: eq.update({'tA':clean(op(0,2)/Yp),'rA':clean(op(1,2)/Yp)})
    angular1,angular2=op(2,2),op(3,3)/s.sin(theta)**2
    eq['K']=clean((angular1+angular2)/(2*r*r*Y))
    if ell>=2:eq['G']=clean((angular1-angular2)/(2*r*r*W))
    all_expr=list(eq.values())+[trace]+list(divergence.values())
    assert all(not expression.has(theta) for expression in all_expr)
    # The general spherical identities predict these constraints independently.
    L=ell*(ell+1)
    expected_t=s.I*w*a+f*s.diff(b,r)+(s.diff(f,r)+2*f/r)*b-L*d/r**2
    expected_r=s.I*w*b/f+s.diff(c,r)+s.diff(f,r)*(a+c)/(2*f)+2*(c-K)/r-L*e/r**2
    assert s.simplify(divergence['0']-expected_t)==0
    assert s.simplify(divergence['1']-expected_r)==0
    if ell>=1:
        expected_A=s.I*w*d/f+f*s.diff(e,r)+(s.diff(f,r)+2*f/r)*e+K-(L-2)*G
        assert s.simplify(divergence['2']-expected_A)==0
    assert s.simplify(trace-(-a+c+2*K))==0
    return names,fields,eq,trace,divergence


def serialize(ell):
    names,fields,eq,trace,divergence=derive(ell)
    jets=[s.diff(fields[n],r,j) if j else fields[n] for j in range(3) for n in names]
    rows=list(eq.values())+[trace,s.diff(trace,r),s.diff(trace,r,2)]
    row_names=list(eq)+['trace','trace_prime','trace_second']
    for name,C in divergence.items():
        rows.extend([C,s.diff(C,r)]);row_names.extend(['div_'+name,'div_'+name+'_prime'])
    matrix,rhs=s.linear_eq_to_matrix(rows,jets)
    assert rhs==s.zeros(len(rows),1)
    return {'ell':ell,'fields':names,'rows':row_names,
            'columns':[f'{n}_{j}' for j in range(3) for n in names],
            'matrix':[[s.sstr(s.factor(v)) for v in row] for row in matrix.tolist()],
            'constraints_verified':True,'angular_dependence_cancelled':True}


if __name__=='__main__':
    import sys
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result={'conventions':'GM=1, exp(-i omega t), mu=m2_pole, q_AB=r^2(K Omega Y+2G Y_AB)',
            'sympy':s.__version__,'systems':[]}
    for ell in (0,1,2,3,4):
        print(f'Deriving polar ell={ell}',file=sys.stderr,flush=True)
        result['systems'].append(serialize(ell))
    # The sphere's covariant Laplacian supplies L=ell(ell+1). Extract its
    # affine coefficient and verify against a third independent projection.
    L=s.Symbol('L')
    M2=s.Matrix(result['systems'][2]['matrix']);M3=s.Matrix(result['systems'][3]['matrix'])
    general=M2+(L-6)*(M3-M2)/6
    M4=s.Matrix(result['systems'][4]['matrix'])
    assert all(s.cancel(v)==0 for v in general.subs(L,20)-M4)
    result['general_l_ge_2']={**result['systems'][2],'ell':'general',
        'matrix':[[s.sstr(s.factor(v)) for v in row] for row in general.tolist()],
        'angular_eigenvalue':'L=ell*(ell+1); affine dependence verified at ell=4'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
