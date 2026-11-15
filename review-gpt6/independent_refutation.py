#!/usr/bin/env python3
"""GPT-6 adversarial exact verifier, independently authored 2026-11-15.

Do not import linear_finite_size_law.py, product_bound.py, their
adjacency-cap or polytope code. Rebuild S3 from permutations; solve
each rational polytope by independent 3-plane intersection; evaluate
all 21 published L=4..6 points and one genuine torus Hamiltonian
configuration (all 96 raw electric moves at L=4).
"""
from fractions import Fraction as F
from itertools import permutations, combinations
from math import isqrt

G=tuple(permutations((0,1,2)))
index={p:i for i,p in enumerate(G)}
unit=index[(0,1,2)]


def require(ok,reason):
    if not ok:
        raise ValueError("INDEPENDENT_LINEAR_REFUSE_"+reason)


def compose(a,b):
    aa,bb=G[a],G[b]
    return index[tuple(aa[bb[i]] for i in range(3))]


def inverse(a):
    p=G[a]
    return index[tuple(p.index(i) for i in range(3))]


def kind(a):
    fixed=sum(G[a][i]==i for i in range(3))
    return 0 if fixed==3 else (1 if fixed==1 else 2)


trans=tuple(i for i in range(6) if kind(i)==1)
cycles=tuple(i for i in range(6) if kind(i)==2)
require(len(trans)==3 and len(cycles)==2,"S3_classes")


def electric_pair(a,b,e,h,flat_special=False):
    """One independently generated 3-transposition electric row."""
    oldT=(kind(a)==1)+(kind(b)==1)
    oldC=(kind(a)==2)+(kind(b)==2)
    result=F(0)
    for t in trans:
        aa,bb=compose(t,a),compose(t,b)
        newT=(kind(aa)==1)+(kind(bb)==1)
        newC=(kind(aa)==2)+(kind(bb)==2)
        require((newT-oldT)%2==0,"frame_transposition_parity")
        factor=e**((newT-oldT)//2)*h**(newC-oldC)
        if flat_special and oldT==2 and a==b==t:
            factor=F(0)
        result+=1-factor
    return result


def weights_from_group(e,h,r):
    w00=electric_pair(unit,unit,e,h)
    w0t=electric_pair(unit,trans[0],e,h)
    w0c=electric_pair(unit,cycles[0],e,h)
    wtc=electric_pair(trans[0],cycles[0],e,h)
    wcc=electric_pair(cycles[0],cycles[1],e,h)
    same=electric_pair(trans[0],trans[0],e,h)
    different=electric_pair(trans[0],trans[1],e,h)
    special=electric_pair(trans[0],trans[0],e,h,True)
    require(same<=different,"same_TT_not_smaller")
    return (2*r+4*(w0t-w00),
            3*r+4*(w0c-w00),
            min(same,different)+w00-2*w0t,
            min(special,different)+w00-2*w0t,
            wtc+w00-w0t-w0c,
            wcc+w00-2*w0c)


def independent_simplification(e,h,r):
    return (2*r-4-8*h+12*e,
            3*r-12*e/h+12*e,
            2-3*e+4*h-(1+2*h*h)/e,
            2-3*e+4*h-(2*h+h*h)/e,
            -1-1/h+2*h+3*e*(1/h-1),
            -3*e*(1/h-1)**2)


def ceil2root(n):
    z=4*n
    v=isqrt(z)
    return v+(v*v<z)


def adjacent(s,L):
    P=L*L
    if s==0: return 0
    if s==P: return 2*P
    return 2*s-min(ceil2root(s),L,ceil2root(P-s))


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def dot(a,b):
    return sum(x*y for x,y in zip(a,b))


def vertices(a,b,L):
    # Independent cross-product solution, not the copied Cramer
    # numerator code in the claimed law and product_bound.py.
    planes=(((2,1,0),4*a),((0,1,2),4*b),
            ((1,1,1),adjacent(a+b,L)),((1,0,0),adjacent(a,L)),
            ((0,0,1),adjacent(b,L)),
            ((-1,0,0),0),((0,-1,0),0),((0,0,-1),0))
    result=set()
    for ii,jj,kk in combinations(range(8),3):
        n1,d1=planes[ii];n2,d2=planes[jj];n3,d3=planes[kk]
        c23,c31,c12=cross(n2,n3),cross(n3,n1),cross(n1,n2)
        determinant=dot(n1,c23)
        if determinant==0: continue
        xyz=tuple(F(d1*c23[i]+d2*c31[i]+d3*c12[i],determinant)
                  for i in range(3))
        if all(dot(normal,xyz)<=limit for normal,limit in planes):
            result.add(xyz)
    require(bool(result),"empty_legal_lump")
    return tuple(result)


def all_lumps(L):
    P=L*L
    return [(a,b,vertices(a,b,L)) for a in range(0,P+1,2)
            for b in range(P-a+1) if a+b]


def prove_cases_and_minimum(L,r,lumps):
    P=L*L
    e=F(9*(L+1),8*L)/r
    h=F(1,4)
    A,B,D,Dprime,Y,Z=weights_from_group(e,h,r)
    require((A,B,D,Dprime,Y,Z)==independent_simplification(e,h,r),
            "S3_weights_do_not_equal_claim")
    require(0<e<F(1,3) and D<=0 and Y<=0 and Z<=0,
            "weight_sign_conditions")
    KT=A+2*D
    KC=B+4*min(Y,Z/2)
    require(KT>0 and KC>=B+4*Y+2*Z==3*r-18-54*e,
            "C_case_step1")
    require(2*A+min(0,Dprime)>=F(7,2)*r-9,
            "special_two_T_not_bounded")
    require(2*A+D>=3*r-9,
            "two_T_with_C_not_bounded")
    require(4*(-D)>=F(16,5)*r-12,
            "middle_a_defects_not_bounded")
    require((P-3)*KT>=3*r,"near_all_T_not_bounded")
    simple=3*r-18-F(27*(P+9)*(L+1),4*L)/r
    closed=min(KC,2*A+min(0,Dprime),4*(-D),(P-3)*KT)-6*P*e
    require(closed>=simple>0,"claimed_strict_gap_fails")
    best=None;where=None
    for a,b,vs in lumps:
        mu=2*a-adjacent(a,L)
        if 4<=a<=P-4:
            require(mu>=4,"middle_defect_adjacency")
        if a==2:
            require(adjacent(2,L)==1,"two_defect_adjacency")
        for x,y,z in vs:
            value=a*A+b*B+x*(Dprime if (a,b)==(2,0) else D)+y*Y+z*Z
            if (a,b)==(2,0):
                lower=2*A+min(0,Dprime)
            elif a==0:
                lower=KC
            elif a==2:
                lower=2*A+D+KC
            elif a<=P-4:
                lower=4*(-D)
            else:
                lower=(P-3)*KT
            require(value>=lower>=3*r-18-54*e,
                    "lump_case_"+str((L,a,b)))
            if best is None or value<best:
                best,where=value,(a,b)
    require(where==(0,1),"unexpected_exact_minimizer")
    return best-6*P*e,closed,simple,where


def plaquettes(links,L):
    P=L*L
    def i(x,y):return (x%L)*L+(y%L)
    result=[]
    for x in range(L):
        for y in range(L):
            q=compose(links[i(x,y)],links[P+i(x+1,y)])
            q=compose(q,inverse(links[i(x,y+1)]))
            q=compose(q,inverse(links[P+i(x,y)]))
            result.append(q)
    return result


def genuine_torus_row():
    # Noncommuting x/y Wilson seams localize a single 3-cycle
    # at the corner of a periodic 4x4 torus.
    L=4;P=L*L;r=F(12);h=F(1,4);e=F(15,128)
    links=[unit]*(2*P)
    def i(x,y):return (x%L)*L+(y%L)
    for y in range(L):
        links[i(L-1,y)]=trans[0]
    for x in range(L):
        links[P+i(x,L-1)]=trans[1]
    initial=plaquettes(links,L)
    a=sum(kind(v)==1 for v in initial)
    b=sum(kind(v)==2 for v in initial)
    require((a,b)==(0,1),"seam_is_not_single_C")
    row=r*(2*a+3*b)
    checked=0
    for link in range(2*P):
        for t in trans:
            moved=links[:]
            moved[link]=compose(t,moved[link])
            after=plaquettes(moved,L)
            a1=sum(kind(v)==1 for v in after)
            b1=sum(kind(v)==2 for v in after)
            require(a1%2==0,"real_config_parity")
            ratio=F(0) if (a1,b1)==(0,0) else (
                e**((a1-a)//2)*h**(b1-b))
            row+=1-ratio
            checked+=1
    require(checked==3*2*P,"electric_moves_not_complete")
    require(row==F(3729,32),"single_C_direct_ratio_mismatch")
    print("REVIEW_REAL|L=4|r=12|configuration=NONCOMMUTING_BOUNDARY_SEAMS"
          "|plaquette_lump=0,1|raw_moves=96|ratio=3729/32"
          "|ratio_minus_flat_bound=657/32|pass=1|json=0",flush=True)


EXPECTED={
 4:(20531250,21065721,24721153,28741071,
    43687500,102843750,478839843),
 5:(21900000,22445642,26182758,32306250,
    51042857,129420000,598744500),
 6:(22950000,23507231,27328125,35625000,
    58218750,155875000,718621875),
}


def run():
    for e,h,r in ((F(15,128),F(1,4),F(12)),
                  (F(1,5),F(2,3),F(19)),
                  (F(1,17),F(7,10),F(13))):
        require(weights_from_group(e,h,r)==
                independent_simplification(e,h,r),"weight_reconstruction")
    print("REVIEW_FRAME|independent_group_S3=PASS"
          "|same_TT_smaller=PASS|special_flat_weight=PASS|json=0",flush=True)
    for L in (4,5,6):
        lumps=all_lumps(L)
        r0=F(3*L,2)+6
        rs=(r0,r0+F(1,8),r0+1,F(2*L+6),
            F(3*L+6),F(9*L),F(40*L))
        for j,r in enumerate(rs):
            val,closed,simple,where=prove_cases_and_minimum(L,r,lumps)
            micro=(val*1000000).__floor__()
            require(micro==EXPECTED[L][j],"published_micro_mismatch")
            print("REVIEW_EXACT|L=%s|r=%s|lumps=%s|minimum_micro=%s"
                  "|at=0,1|closed_micro=%s|simple_micro=%s|pass=1|json=0"
                  %(L,r,len(lumps),micro,
                    (closed*1000000).__floor__(),
                    (simple*1000000).__floor__()),flush=True)
    for L in (4,5,6,7,8,16,64,200,1000):
        r=F(3*L,2)+6
        require(L*L-3*L-3>0,"L4_threshold_algebra")
        e=F(9*(L+1),8*L)/r
        gap=3*r-18-F(27*(L*L+9)*(L+1),4*L)/r
        require(0<e<F(1,3) and gap>0,"formula_sign")
    genuine_torus_row()
    print("REVIEW_DONE|S3_weights=PASS|cases=PASS"
          "|exact_L4_L5_L6_points=21|real_4x4=PASS"
          "|threshold_algebra=PASS|counterexample=NONE_FOUND"
          "|prior_art=UNVERIFIED|whole_Q_from_prior=INHERITED|json=0",flush=True)


if __name__=="__main__":
    run()
