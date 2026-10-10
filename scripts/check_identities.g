# Exact checks of three-copy identities for n = 1 and n = 2.
# Run from the repository root: gap -q -b --quitonbreak scripts/check_identities.g
# Only GAP core is needed. All scalars are rational numbers or cyclotomics.
# Sparse matrices use lists of columns; a column is a sorted list [row, value].
# Thus the n = 2 calculation never performs dense 729 x 729 multiplication.
# This finite check supports the proofs; it does not prove statements for all n.

CheckThreeCopyIdentities := function()
    local check, count, addColumn, plus, scale, product, adjoint, identity,
          tensor, cube, trace, equal, index, dot, weyl, permutation, css,
          n, d, dim, labels, omega, id, zero, s, t, p, ps, piS, piA, piM,
          q, rp, rm, kraw, ranks, expected, a, b, c, e, w, x,
          ws, points, parity, averageW, averageA, j, h, pairing, label,
          sumW, pointFourier, rho, rho3, t2, t3, um, vm, weights, chi,
          wigner, val, l1, l2, h1, h2, xx, yy, xp, yp, aa, bb, det1,
          det2, jac, z;

    count := 0;
    check := function(condition, message)
        count := count + 1;
        if not condition then
            Error(Concatenation("Identity check failed: ", message));
        fi;
    end;

    addColumn := function(a, b)
        local out, i, j, value;
        out := []; i := 1; j := 1;
        while i <= Length(a) or j <= Length(b) do
            if j > Length(b) or (i <= Length(a) and a[i][1] < b[j][1]) then
                Add(out, a[i]); i := i + 1;
            elif i > Length(a) or b[j][1] < a[i][1] then
                Add(out, b[j]); j := j + 1;
            else
                value := a[i][2] + b[j][2];
                if value <> 0 then Add(out, [a[i][1], value]); fi;
                i := i + 1; j := j + 1;
            fi;
        od;
        return out;
    end;
    plus := function(a, b)
        return List([1..Length(a)], i -> addColumn(a[i], b[i]));
    end;
    scale := function(a, c)
        if c = 0 then return List(a, col -> []); fi;
        return List(a, col -> List(col, entry -> [entry[1], c * entry[2]]));
    end;
    product := function(a, b)
        local out, col, entry, value;
        out := [];
        for col in b do
            value := [];
            for entry in col do
                value := addColumn(value,
                    List(a[entry[1]], x -> [x[1], entry[2] * x[2]]));
            od;
            Add(out, value);
        od;
        return out;
    end;
    adjoint := function(a)
        local out, j, entry;
        out := List(a, col -> []);
        for j in [1..Length(a)] do
            for entry in a[j] do
                Add(out[entry[1]], [j, ComplexConjugate(entry[2])]);
            od;
        od;
        return out;
    end;
    identity := dim -> List([1..dim], j -> [[j, 1]]);
    tensor := function(a, b)
        local out, ac, bc, av, bv, col, db;
        out := []; db := Length(b);
        for ac in a do
            for bc in b do
                col := [];
                for av in ac do
                    for bv in bc do
                        Add(col, [(av[1]-1)*db+bv[1], av[2]*bv[2]]);
                    od;
                od;
                Add(out, col);
            od;
        od;
        return out;
    end;
    cube := a -> tensor(tensor(a, a), a);
    trace := function(a)
        local answer, j, entry;
        answer := 0;
        for j in [1..Length(a)] do
            for entry in a[j] do
                if entry[1] = j then answer := answer + entry[2]; fi;
            od;
        od;
        return answer;
    end;
    equal := function(a, b, message) check(a = b, message); end;
    index := function(x)
        local answer, entry;
        answer := 0;
        for entry in x do answer := 3*answer + (entry mod 3); od;
        return answer + 1;
    end;
    dot := function(a, b)
        return Sum([1..Length(a)], j -> a[j]*b[j]);
    end;
    omega := E(3);
    weyl := function(a, b, labels)
        # W_(a,b)|x> = omega^(a.x + 2 a.b)|x+b>.
        return List(labels, x -> [[index(x+b), omega^(dot(a,x)+2*dot(a,b))]]);
    end;
    permutation := function(d, order)
        local out, x, y;
        out := [];
        for x in Tuples([0..d-1], 3) do
            y := x{order};
            Add(out, [[d*d*y[1]+d*y[2]+y[3]+1, 1]]);
        od;
        return out;
    end;
    css := function(labels)
        local d, out, x, y, support, row, col;
        d := Length(labels);
        out := List([1..d^3], j -> []);
        for x in labels do
            support := [];
            for y in labels do
                Add(support, (index(x+y)-1)*d*d+(index(y-x)-1)*d+index(y));
            od;
            Sort(support);
            for col in support do
                out[col] := List(support, row -> [row, 1/d]);
            od;
        od;
        return out;
    end;

    Print("Exact three-copy checks (GAP ", GAPInfo.Version, ")\n");
    for n in [1,2] do
        d := 3^n; dim := d^3; labels := Tuples([0,1,2], n);
        id := identity(dim); zero := scale(id, 0);
        s := permutation(d, [2,1,3]); t := permutation(d, [3,1,2]);
        p := css(labels); ps := product(p, s);

        equal(product(s,s), id, "S^2 = I");
        equal(product(product(t,t),t), id, "T^3 = I");
        equal(product(product(s,t),s), product(t,t), "STS = T^2");
        equal(product(p,p), p, "P^2 = P");
        equal(adjoint(p), p, "P is Hermitian");
        equal(product(p,t), p, "PT = P");
        equal(product(s,p), ps, "SP = PS");
        equal(trace(p), d, "rank P = 3^n");
        equal(trace(ps), 1, "Tr(PS) = 1");

        piS := scale(plus(plus(plus(id,t),product(t,t)),
                    plus(plus(s,product(t,s)),product(product(t,t),s))), 1/6);
        piA := scale(plus(plus(plus(id,t),product(t,t)),
                    scale(plus(plus(s,product(t,s)),product(product(t,t),s)),-1)), 1/6);
        piM := plus(id, scale(plus(piS,piA),-1));
        q := [scale(plus(p,ps),1/2)];
        Add(q, plus(piS,scale(q[1],-1)));
        Add(q, scale(plus(p,scale(ps,-1)),1/2));
        Add(q, plus(piA,scale(q[3],-1)));
        Add(q, piM);
        z := zero;
        for j in [1..5] do
            equal(adjoint(q[j]),q[j],"Q_j is Hermitian");
            equal(product(q[j],q[j]),q[j],"Q_j^2 = Q_j");
            for h in [1..j-1] do
                equal(product(q[j],q[h]),zero,"Q_j Q_k = 0 for j <> k");
            od;
            z := plus(z,q[j]);
        od;
        equal(z,id,"sum of five sector projectors = I");
        equal(product(p,piM),zero,"P Pi_m = 0");
        ranks := List(q,trace);
        expected := [(d+1)/2,(d^2-1)*(d+3)/6,(d-1)/2,
                     (d^2-1)*(d-3)/6,2*d*(d^2-1)/3];
        equal(ranks,expected,"five sector ranks");
        rp := scale(product(piM,plus(id,s)),1/2);
        rm := scale(product(piM,plus(id,scale(s,-1))),1/2);
        for x in [rp,rm] do
            equal(product(x,x),x,"R_+ and R_- are projections");
            equal(adjoint(x),x,"R_+ and R_- are Hermitian");
            equal(trace(x),d*(d^2-1)/3,"primitive mixed-sector ranks");
        od;
        equal(product(rp,rm),zero,"R_+ R_- = 0");
        equal(plus(rp,rm),piM,"R_+ + R_- = Pi_m");
        kraw := plus(t,scale(product(t,t),-1));
        # These are the identities for K=(T-T^2)/sqrt(3), with sqrt(3) cleared.
        equal(adjoint(kraw),scale(kraw,-1),"(T-T^2)^* = -(T-T^2)");
        equal(product(kraw,kraw),scale(piM,-3),"(T-T^2)^2 = -3 Pi_m");
        equal(product(product(s,kraw),s),scale(kraw,-1),"S(T-T^2)S = -(T-T^2)");
        equal(product(kraw,rp),product(rm,kraw),"intertwiner between mixed copies");

        ws := [];
        for a in labels do for b in labels do Add(ws,weyl(a,b,labels)); od; od;
        for j in [1..Length(ws)] do
            a := labels[QuoInt(j-1,d)+1]; b := labels[((j-1) mod d)+1];
            equal(product(adjoint(ws[j]),ws[j]),identity(d),"Weyl unitarity");
            for h in [1..Length(ws)] do
                c := labels[QuoInt(h-1,d)+1]; e := labels[((h-1) mod d)+1];
                pairing := dot(a,e)-dot(b,c);
                label := (index(a+c)-1)*d+index(b+e);
                equal(product(ws[j],ws[h]),scale(ws[label],omega^(2*pairing)),
                      "W_u W_v = omega^([u,v]/2) W_(u+v)");
            od;
        od;
        parity := List(labels,x -> [[index(-x),1]]);
        sumW := scale(identity(d),0);
        for w in ws do sumW := plus(sumW,w); od;
        equal(scale(sumW,1/d),parity,"J = D^(-1) sum_u W_u");
        points := List(ws,w -> product(product(w,parity),adjoint(w)));
        averageW := zero; averageA := zero;
        for j in [1..Length(ws)] do
            averageW := plus(averageW,scale(cube(ws[j]),1/d^2));
            averageA := plus(averageA,scale(cube(points[j]),1/d^2));
            a := labels[QuoInt(j-1,d)+1]; b := labels[((j-1) mod d)+1];
            pointFourier := scale(identity(d),0);
            for h in [1..Length(ws)] do
                c := labels[QuoInt(h-1,d)+1]; e := labels[((h-1) mod d)+1];
                pointFourier := plus(pointFourier,
                    scale(ws[h],omega^(dot(a,e)-dot(b,c))/d));
            od;
            equal(points[j],pointFourier,"phase-point Fourier expansion");
        od;
        equal(averageW,p,"P = D^(-2) sum_u W_u^(tensor 3)");
        equal(averageA,ps,"PS = D^(-2) sum_u A_u^(tensor 3)");

        # Independent exact density example: the normalized negative-parity projector.
        rho := scale(plus(identity(d),scale(parity,-1)),1/(d-1));
        rho3 := cube(rho);
        t2 := trace(product(rho,rho)); t3 := trace(product(product(rho,rho),rho));
        um := trace(product(p,rho3)); vm := trace(product(ps,rho3));
        val := (d-3)/(d*(d-1)^2);
        equal([t2,t3,um,vm],[2/(d-1),4/(d-1)^2,val,-val],
              "negative-parity density invariants");
        chi := List(ws,w -> trace(product(w,rho)));
        wigner := List(points,a -> trace(product(a,rho))/d);
        equal(Sum(chi,x -> x^3)/d^2,um,"Weyl cubic moment");
        equal(d*Sum(wigner,x -> x^3),vm,"Wigner cubic moment");
        equal(Sum(wigner),1,"Wigner normalization");
        equal(wigner,Concatenation([-1/d],List([2..d^2],j -> 1/(d*(d-1)))),
              "negative-parity Wigner values");
        weights := [(um+vm)/2,(1+3*t2+2*t3)/6-(um+vm)/2,
                    (um-vm)/2,(1-3*t2+2*t3)/6-(um-vm)/2,2*(1-t3)/3];
        equal(List(q,x -> trace(product(x,rho3))),weights,"five density sector weights");
        if n = 2 then equal(weights,[0,5/16,1/96,5/96,5/8],"n=2 density example"); fi;
        Print("n=",n,": sector ranks ",ranks,"; Weyl, CSS, mixed-sector and density checks passed.\n");
    od;

    # Appendix A: differentiate its displayed polynomials and evaluate the two
    # diagonal Jacobian blocks using exact rational arithmetic.
    l1 := 1/2; l2 := 1/3; h1 := 1/4; h2 := 1/4;
    xx := l -> (1+2*l^2)/3; yy := l -> (1+6*l^2+2*l^3)/9;
    xp := l -> 4*l/3; yp := l -> (4*l+2*l^2)/3;
    aa := function(l,h) return (1-l^3)/9+l^3*h; end;
    bb := function(l,h) return aa(l,h)-2*(1-l)^2*(1+2*l)/27; end;
    det1 := xp(l1)*xx(l2)*yy(l1)*yp(l2)-xx(l1)*xp(l2)*yp(l1)*yy(l2);
    det2 := l1^3*l2^3*(aa(l2,h2)*bb(l1,h1)-aa(l1,h1)*bb(l2,h2));
    jac := det1*det2;
    equal([det1,det2,jac],[142/19683,1/78732,71/774840978],"Appendix A Jacobian");
    Print("Appendix A: det J_11 = ",det1,"; det J_22 = ",det2,"; det J = ",jac,".\n");
    Print("PASS: ",count," exact equalities. No floating-point tolerance was used.\n");
    return true;
end;

CheckThreeCopyIdentities();
QUIT_GAP(0);
