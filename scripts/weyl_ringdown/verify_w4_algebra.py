"""Exact null-source elimination and radial/tortoise convention audit."""
import json
import sympy as s
from w3_polar import DEFAULT_OPERATORS,monopole_coefficients


def verify():
    record=json.loads(DEFAULT_OPERATORS.read_text())['systems'][0]
    matrix=s.Matrix(record['matrix']);jets=s.Matrix(s.symbols(' '.join(record['columns'])))
    r,w,mu=s.symbols('r omega mu');F=1-2/r
    a,b,K,ap,bp,Kp,app,bpp,Kpp=s.symbols('a_0 b_0 K_0 a_1 b_1 K_1 a_2 b_2 K_2')
    c,cp,cpp=s.symbols('c_0 c_1 c_2')
    derivative=lambda z:s.diff(z,r)+s.diff(z,a)*ap+s.diff(z,b)*bp+s.diff(z,K)*Kp+s.diff(z,Kp)*Kpp
    # The 2 K' term must be retained before differentiating this constraint.
    ap_expr=2*Kp-s.I*w*b/F-s.diff(F,r)*(a-K)/F-2*(a-3*K)/r
    bp_expr=-s.I*w*a/F-(s.diff(F,r)/F+2/r)*b
    eq=[]
    for e in (matrix*jets)[:3]:
        e=e.subs({c:a-2*K,cp:ap-2*Kp,cpp:app-2*Kpp})
        e=e.subs({app:derivative(ap_expr),bpp:derivative(bp_expr)})
        e=e.subs({ap:ap_expr,bp:bp_expr})
        # Isolate the forcing: K=K'=0; the angular wave equation gives K''.
        eq.append(s.cancel(e.subs({K:0,Kp:0,Kpp:-2*a/r**2})))
    M,y=s.linear_eq_to_matrix([eq[1]-1/(r*r*F),eq[2]-1/(r*r*F*F)],(a,b))
    det=s.cancel(M.det())
    aa=s.cancel((y[0]*M[1,1]-y[1]*M[0,1])/det)
    bb=s.cancel((M[0,0]*y[1]-M[1,0]*y[0])/det)
    source=s.factor(-2*aa/r**2)
    D=r*(2*mu**2*F+mu**4*r*r+4*w*w)
    expected=2*(mu*mu*r+2*s.I*w)/(r*r*F*D)
    assert s.cancel(source-expected)==0
    assert s.cancel(eq[0].subs({a:aa,b:bb})-1/r**2)==0
    P,V=monopole_coefficients(r,w,mu)
    assert s.cancel(P-(s.diff(F,r)/F+5/r-s.diff(D,r)/D))==0
    Q0=4*mu*mu*r*r-mu**4*r**5-2*r**3*(2*w*w+mu*mu)
    assert s.cancel(Q0+r*r*D)==0
    # Strip the paper's 4 sqrt(pi) E factor; Eq.(35) is an r_* source.
    paper=-2*F*(mu*mu*r+2*s.I*w)/Q0
    assert s.cancel(paper-F*F*source)==0
    horizon=s.cancel(F*F*source/(r-2)).subs(r,2)
    assert s.cancel(horizon-1/(16*(mu*mu-s.I*w)))==0
    return {'null_source_from_W3':'exact identity; all three remaining wave rows consistent',
            'radial_source_per_kappa_E_Y00_phase':s.sstr(source),
            'Wronskian_r':'D/(r^5 F), up to a constant',
            'Wronskian_tortoise':'D/r^5, up to the same constant',
            'AGP_equation_35_mapping':'equals F^2 times radial source after stripping 4 sqrt(pi) E; symmetric Fourier prefactor must be handled separately',
            'horizon_forcing_coefficient':'1/[16 (mu^2-i omega)]',
            'passed':True}


if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
