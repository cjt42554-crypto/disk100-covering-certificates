"""Exact standard-library SOSC acceptance for the frozen disk100 KKT root.

The proposal provides only pivots, a dyadic inverse Y and a dyadic Cholesky
transform W. Root, geometry, constraints, Hessian, residuals and all norm bounds
are freshly checked. No floating point enters any acceptance inequality.
"""
import json,gzip,hashlib,time,sys,io,contextlib
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,localcontext
ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'source'

def require(condition,message):
    if not condition:raise ValueError(message)

def qlist(q):return [q.numerator,q.denominator]

def dec(q):
    with localcontext() as c:
        c.prec=40
        return str(Decimal(q.numerator)/Decimal(q.denominator))

def norm(a,den=1):return Q(max(sum(abs(v) for v in row) for row in a),den)

def tnorm(a,den=1):return Q(max(sum(abs(row[j]) for row in a) for j in range(len(a[0]))),den)

def dense_product(a,b):
    bt=list(zip(*b))
    return [[sum(x*y for x,y in zip(row,col) if x and y) for col in bt] for row in a]

def exact_derivatives(system,mid):
    n=system['n_base'];m=system['n_constraints'];D=mid['denominator'];x=[int(v) for v in mid['numerators']]
    require(type(D) is int and D>0,'Bad midpoint denominator')
    require(len(x)==system['n_variables']==n+m,'Bad midpoint dimension')
    B=[{} for _ in range(m)];Blip=[{} for _ in range(m)]
    H=[{} for _ in range(n)];Hlip=[{} for _ in range(n)]
    for row,terms in enumerate(system['polynomials'][:m]):
        for a,vs in terms:
            if not vs:continue
            require(type(a) is int,'Derivative coefficient must be integer')
            require(all(v<n for v in vs),'Constraint depends on a multiplier')
            if len(vs)==1:
                j=vs[0];B[row][j]=B[row].get(j,0)+a*D
            elif len(vs)==2:
                i,j=vs
                for u,v in ((i,j),(j,i)):
                    B[row][u]=B[row].get(u,0)+a*x[v]
                    Blip[row][u]=Blip[row].get(u,0)+abs(a)
                    H[u][v]=H[u].get(v,0)+a*x[n+row]
                    Hlip[u][v]=Hlip[u].get(v,0)+abs(a)
            else:raise ValueError('Constraint degree exceeds two')
    for i,row in enumerate(H):
        for j,v in row.items():require(H[j].get(i,0)==v,'Hessian is not symmetric')
    return D,x,B,Blip,H,Hlip

def verify_second_order(system,mid,witness,rho):
    start=time.time();n=system['n_base'];m=system['n_constraints'];k=n-m
    require(witness['format']=='disk100-sosc-dyadic-v1','Unknown witness format')
    require(witness['n_base']==n and witness['n_constraints']==m,'Witness dimension mismatch')
    p=witness['pivot_columns'];f=witness['free_columns']
    require(len(p)==m and len(f)==k and all(type(v) is int for v in p+f),'Bad pivot dimensions')
    require(sorted(p+f)==list(range(n)),'Pivot and free columns must form a partition')
    P=witness['inverse_denominator'];Y=witness['inverse_matrix']
    V=witness['cholesky_transform_denominator'];W=witness['cholesky_transform']
    require(type(P) is int and P>0 and type(V) is int and V>0,'Bad dyadic denominator')
    require(P&(P-1)==0 and V&(V-1)==0,'Denominators must be powers of two')
    require(len(Y)==m and all(len(row)==m and all(type(v) is int for v in row) for row in Y),'Bad inverse matrix')
    require(len(W)==k and all(len(row)==k and all(type(v) is int for v in row) for row in W),'Bad Cholesky transform')
    D,x,B,Blip,H,Hlip=exact_derivatives(system,mid)
    PD=P*D
    cols=[[] for _ in range(n)]
    for i,row in enumerate(B):
        for j,v in row.items():
            if v:cols[j].append((i,v))
    residual=[0]*m
    for j,c in enumerate(p):
        for i,row in enumerate(Y):
            value=sum(row[l]*v for l,v in cols[c])-(PD if i==j else 0)
            residual[i]+=abs(value)
    epsilon=Q(max(residual),PD)
    ynorm=norm(Y,P)
    lp=max(sum(row.get(j,0) for j in p) for row in Blip)
    lf=max(sum(row.get(j,0) for j in f) for row in Blip)
    epstar=epsilon+ynorm*lp*rho
    require(epstar<1,'Pivot inverse Neumann condition failed')
    inv_bound=ynorm/(1-epstar)
    Z=[[0]*k for _ in range(n)]
    for i,ci in enumerate(p):
        row=Y[i]
        for j,cj in enumerate(f):Z[ci][j]=-sum(row[l]*v for l,v in cols[cj])
    for j,cj in enumerate(f):Z[cj][j]=PD
    pnorm=norm([Z[ci] for ci in p],PD)
    znorm=norm(Z,PD);ztnorm=tnorm(Z,PD)
    dz=(epstar*pnorm+ynorm*lf*rho)/(1-epstar)
    dzt=m*dz
    hnorm=Q(max(sum(abs(v) for v in row.values()) for row in H),D)
    hlip=max(sum(row.values()) for row in Hlip)
    dh=hlip*rho
    dq=dzt*hnorm*znorm+ztnorm*hnorm*dz+dzt*hnorm*dz+(ztnorm+dzt)*dh*(znorm+dz)
    HZ=[[sum(v*Z[j][c] for j,v in row.items()) for c in range(k)] for row in H]
    Qnum=dense_product(list(map(list,zip(*Z))),HZ)
    Qden=PD*PD*D
    require(all(Qnum[i][j]==Qnum[j][i] for i in range(k) for j in range(k)),'Reduced Hessian is not symmetric')
    Snum=dense_product(dense_product(W,Qnum),list(map(list,zip(*W))))
    Sden=V*V*Qden
    eps_s=Q(max(sum(abs(v-(Sden if i==j else 0)) for j,v in enumerate(row)) for i,row in enumerate(Snum)),Sden)
    wn=norm(W,V);wtn=tnorm(W,V)
    pert=wn*wtn*dq
    total=eps_s+pert
    require(total<1,'Verified transformed reduced Hessian is not strictly positive definite')
    ambient_sigma=(1-total)/(wn*wtn*(znorm+dz)*(ztnorm+dzt))
    quantities={'pivot_inverse_residual':epsilon,'pivot_neumann_uniform_bound':epstar,
        'pivot_inverse_infinity_bound':inv_bound,'null_basis_midpoint_infinity_norm':znorm,
        'null_basis_midpoint_transpose_infinity_norm':ztnorm,'null_basis_error_infinity_bound':dz,
        'null_basis_error_transpose_infinity_bound':dzt,'hessian_midpoint_infinity_norm':hnorm,
        'hessian_root_perturbation_infinity_bound':dh,'reduced_hessian_perturbation_infinity_bound':dq,
        'cholesky_transform_infinity_norm':wn,'cholesky_transform_transpose_infinity_norm':wtn,
        'transformed_midpoint_identity_residual':eps_s,'transformed_root_perturbation_bound':pert,
        'transformed_total_identity_residual':total,'transformed_smallest_eigenvalue_lower_bound':1-total,
        'ambient_tangent_quadratic_form_lower_bound':ambient_sigma}
    return {'verified':True,'n_base':n,'n_constraints':m,'tangent_dimension':k,
        'exact_bounds':{name:qlist(v) for name,v in quantities.items()},
        'decimal_diagnostics':{name:dec(v) for name,v in quantities.items()},'elapsed_seconds':time.time()-start}

def main():
    started=time.time()
    binding=json.loads((ROOT/'SOURCE_BINDING.json').read_text())
    manifest_bytes=(ROOT/binding['upstream_manifest_file']).read_bytes()
    require(hashlib.sha256(manifest_bytes).hexdigest()==binding['upstream_manifest_sha256'],'Upstream frozen manifest binding failed')
    upstream_manifest=json.loads(manifest_bytes)
    for name,sha in binding['files'].items():
        require(sha==upstream_manifest['files'][name]['sha256'],'Input differs from frozen upstream manifest: '+name)
        require(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==sha,'Source binding failed: '+name)
    system=json.loads((SOURCE/'polynomial_system.json').read_text())
    witness=json.loads(gzip.decompress((ROOT/'sosc_witness.json.gz').read_bytes()))
    require(witness['source_system_sha256']==binding['files']['polynomial_system.json'],'Witness source mismatch')
    sys.path.insert(0,str(SOURCE))
    import exact_krawczyk,exact_algebraic_cover
    require(exact_krawczyk.ROOT==SOURCE and exact_algebraic_cover.ROOT==SOURCE,'Source verifier path mismatch')
    with contextlib.redirect_stdout(io.StringIO()):
        root_report=exact_krawczyk.verify(generate_requested=False,emit=False)
        cover_report=exact_algebraic_cover.main()
    mid=json.loads((SOURCE/'rational_midpoint.json').read_text())
    rho=Q(*root_report['box_radius'])
    report=verify_second_order(system,mid,witness,rho)
    # Every hole must define a smooth local function of all center coordinates.
    I=exact_algebraic_cover.I;sub=exact_algebraic_cover.sub;cross=exact_algebraic_cover.cross
    D=mid['denominator'];interval=[I(Q(int(v),D)-rho,Q(int(v),D)+rho) for v in mid['numerators']]
    nt=len(system['internal_active_triangles']);nc=system['n_centers'];n=system['n_base']
    C=[interval[2*i:2*i+2] for i in range(nc)]
    trans=[]
    for k,(a,b) in enumerate(system['boundary_active_pairs']):
        q=interval[2*nc+2*(nt+k):2*nc+2*(nt+k)+2]
        det=4*cross(sub(C[b],C[a]),q)
        require(det.lo>0 or det.hi<0,'Boundary hole IFT determinant contains zero')
        trans.append(min(abs(det.lo),abs(det.hi)))
    require(not system['unused_center_anchors'],'Unanchored full center participation required')
    gx=C[system['gauge_center']][0]
    require(gx.lo>0 or gx.hi<0,'Rotation gauge transversality failed')
    distance=[];collapsed={}
    for i,meta in enumerate(system['constraint_metadata']):
        if meta['kind']=='distance':
            value=interval[n+i]
            require(value.lo>0,'Distance multiplier is not strictly positive')
            distance.append(value.lo);collapsed[meta['hole']]=collapsed.get(meta['hole'],I(0))+value
    require(len(collapsed)==nt+len(system['boundary_active_pairs']),'Contact grouping incomplete')
    require(all(v.lo>0 for v in collapsed.values()),'Collapsed contact multiplier is not positive')
    report.update({'statement':'The certified 100-center covering is a strict local minimizer of radius modulo simultaneous rotation. No global optimality conclusion.',
        'source_binding':binding,'certified_radius_interval':cover_report['certified_radius_interval'],
        'radius_interval_decimal':cover_report['radius_interval_decimal'],
        'boundary_hole_ift_determinant_absolute_lower_bound':qlist(min(trans)),
        'minimum_distance_multiplier_lower_bound':qlist(min(distance)),
        'minimum_collapsed_contact_multiplier_lower_bound':qlist(min(v.lo for v in collapsed.values())),
        'root_freshly_verified':root_report['verified'],'whole_closed_disk_freshly_verified':cover_report['verified'],
        'total_elapsed_seconds':time.time()-started})
    (ROOT/'replay_verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    return report

if __name__=='__main__':
    try:
        require(len(sys.argv)==1,'This acceptance entry point takes no arguments')
        main()
    except Exception as exc:
        print(json.dumps({'verified':False,'error':type(exc).__name__+': '+str(exc)}))
        raise SystemExit(1)
