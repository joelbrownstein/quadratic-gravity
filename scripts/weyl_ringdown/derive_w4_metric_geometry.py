"""Independent coordinate derivation of Ricci, mass and observer tidal maps."""
from functools import lru_cache
from pathlib import Path
import json
import sympy as s
from verify_w1_algebra import r,theta,w,coords,f,g,gi,gamma,riemann,harmonic_diff
from w3_polar_system import clean


def derive():
    names=('T','B','C','D')
    T,B,C,D=[s.Function(n)(r) for n in names]
    h=s.zeros(4);h[0,0]=T;h[0,1]=h[1,0]=B;h[1,1]=C
    h[2,2]=r*r*D;h[3,3]=r*r*s.sin(theta)**2*D
    @lru_cache(None)
    def nab(a,b,c):
        return clean(harmonic_diff(h[b,c],a)-sum(gamma(z,a,b)*h[z,c]+gamma(z,a,c)*h[b,z] for z in range(4)))
    @lru_cache(None)
    def dg(a,b,c):
        return clean(gi[a,a]*(nab(b,a,c)+nab(c,a,b)-nab(a,b,c))/2)
    def dric(b,c):
        ans=sum(harmonic_diff(dg(a,b,c),a)-harmonic_diff(dg(a,a,b),c) for a in range(4))
        ans+=sum(dg(a,a,z)*gamma(z,b,c)+gamma(a,a,z)*dg(z,b,c)
                 -dg(a,c,z)*gamma(z,a,b)-gamma(a,c,z)*dg(z,a,b)
                 for a in range(4) for z in range(4))
        return clean(ans)
    def dR(a,b,c,d):
        raised=harmonic_diff(dg(a,d,b),c)-harmonic_diff(dg(a,c,b),d)
        raised+=sum(dg(a,c,z)*gamma(z,d,b)+gamma(a,c,z)*dg(z,d,b)
                    -dg(a,d,z)*gamma(z,c,b)-gamma(a,d,z)*dg(z,c,b) for z in range(4))
        return clean(g[a,a]*raised+sum(h[a,z]*riemann(z,b,c,d) for z in range(4)))
    Rtt,Rtr,Rrr,Rang=dric(0,0),dric(0,1),dric(1,1),clean(dric(2,2)/r**2)
    print('Derived independent Ricci components',flush=True)
    scalar=clean(-Rtt/f+f*Rrr+2*Rang)
    mass=clean(D/2+r*f*f*C/2-r*f*(D+r*s.diff(D,r))/2)
    # Orthnormalization, then displacement to the same physical areal radius.
    Er=clean(dR(0,1,0,1)-2/r**3*(T/f-f*C)-3*D/r**3)
    Et=clean(dR(0,2,0,2)/(f*r*r)+T/(f*r**3)+D/(2*r**3))
    assert clean(Er-(2*Rang-2*mass/r**3-scalar/2))==0
    assert clean(Et-(Rtt/(2*f)-Rang+mass/r**3+scalar/4))==0
    print('Checked relational tidal identities',flush=True)
    fields=[T,B,C,D]
    jets=[s.diff(q,r,j) for j in range(3) for q in fields]
    columns=[n+'_'+str(j) for j in range(3) for n in names]
    rows={'Ricci_tt':Rtt,'Ricci_tr':Rtr,'Ricci_rr':Rrr,'Ricci_angular':Rang,
          'Ricci_scalar':scalar,'mass':mass,'tidal_radial':Er,'tidal_transverse':Et}
    matrix,rhs=s.linear_eq_to_matrix(list(rows.values()),jets)
    assert rhs==s.zeros(len(rows),1)
    # Pure radial/time coordinate changes must leave the relational outputs zero.
    xt,xr=s.Function('xt')(r),s.Function('xr')(r)
    gauge=[-s.diff(f,r)*xr+2*s.I*w*f*xt,
           -f*s.diff(xt,r)-s.I*w*xr/f,
           -s.diff(f,r)*xr/f**2+2*s.diff(xr,r)/f,2*xr/r]
    gj=s.Matrix([s.diff(q,r,j) for j in range(3) for q in gauge])
    assert all(clean(e)==0 for e in matrix*gj)
    result={'conventions':'signature -+++; W0 Riemann; fixed areal radius; instantaneous Kodama observer',
            'fields':list(names),'columns':columns,'rows':list(rows),
            'matrix':[[s.sstr(x) for x in row] for row in matrix.tolist()],
            'checks':{'general_geometric_tidal_identities':True,'pure_gauge_Ricci_mass_and_relational_tides':True}}
    path=Path(__file__).resolve().parents[2]/'notes/validation/w4_metric_geometry.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('Saved',path,flush=True)


if __name__=='__main__':derive()
