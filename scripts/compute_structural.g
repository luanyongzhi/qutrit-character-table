# Independent little-group computation for n=2. No checkpoints are reused.
# Usage: gap -q -A -b --quitonbreak scripts/compute_structural.g
if not IsBound(QutritTablesRoot) then QutritTablesRoot := "."; fi;
Read(Filename(Directory(QutritTablesRoot), "init.g"));
Read(QutritPath("scripts/affine.g"));

QutritComputeStructuralN2 := function()
    local n, start, f, m, vecs, powers, index, matrixGroup, qgens, q, p, g,
          classData, qd, pd, classes, ys, queue, qi, u, a, v, tdual,
          qkeys, pkeys, qcache, pcache, classPosition, columns, ci,
          element, translation, matrix, dual, qp, weights, h, j, rows,
          d, t, frozen, output;
    n := 2; start := Runtime(); f := GF(3); m := 2*n;
    Print("STRUCTURAL_START n=2 GAP=", GAPInfo.Version, "\n");
    vecs := Tuples(List([0..2], i -> i*One(f)), m);
    powers := List([1..m], i -> 3^(m-i));
    index := v -> 1 + List(v, IntFFE)*powers;
    matrixGroup := Sp(m, 3);
    # Contragredient action on additive characters indexed by row vectors.
    qgens := List(GeneratorsOfGroup(matrixGroup),
        a -> PermList(List(vecs, v -> index(v*TransposedMat(a^-1)))));
    q := Group(qgens); p := Stabilizer(q, 2); g := QutritAffineGroup(n);
    if Size(g) <> QutritExpectedOrder(n) or Size(q) <> 51840
       or Size(p) <> 648 then Error("computed group or stabilizer order failed"); fi;
    classData := function(group)
        local cc, ct, rr;
        cc := ConjugacyClasses(group); ct := CharacterTable(group);
        rr := List(Irr(ct), ValuesOfClassFunction);
        return rec(generators := GeneratorsOfGroup(group),
            representatives := List(cc, Representative),
            classSizes := List(cc, Size),
            classOrders := List(cc, c -> Order(Representative(c))),
            irreducibles := rr, order := Size(group));
    end;
    qd := classData(q); pd := classData(p);
    if Length(qd.irreducibles) <> 34 or Length(pd.irreducibles) <> 24 then
        Error("unexpected auxiliary character count");
    fi;
    Print("AUXILIARIES quotient=34 stabilizer=24\n");
    classes := ConjugacyClasses(g);
    # ys[u] maps the chosen nonzero dual vector (point 2) to point u.
    ys := []; ys[2] := (); queue := [2]; qi := 1;
    while qi <= Length(queue) do
        u := queue[qi];
        for a in GeneratorsOfGroup(q) do
            v := u^a;
            if not IsBound(ys[v]) then
                ys[v] := ys[u]*a; Add(queue, v);
            fi;
        od;
        qi := qi+1;
    od;
    if Length(queue) <> 3^m-1 then Error("nonzero dual orbit is not transitive"); fi;
    qkeys := List(qd.representatives, CycleStructurePerm);
    pkeys := List(pd.representatives, CycleStructurePerm);
    qcache := NewDictionary((), true); pcache := NewDictionary((), true);
    classPosition := function(element, group, data, keys, cache)
        local pos;
        pos := LookupDictionary(cache, element);
        if pos <> fail then return pos; fi;
        for pos in Positions(keys, CycleStructurePerm(element)) do
            if IsConjugate(group, data.representatives[pos], element) then
                AddDictionary(cache, element, pos); return pos;
            fi;
        od;
        Error("auxiliary class matching failed");
    end;
    columns := [];
    for ci in [1..Length(classes)] do
        element := Representative(classes[ci]);
        translation := vecs[1^element];
        matrix := List([1..m], i -> vecs[(1+powers[i])^element]-translation);
        dual := TransposedMat(matrix^-1);
        tdual := PermList(List(vecs, a -> index(a*dual)));
        qp := classPosition(tdual, q, qd, qkeys, qcache);
        weights := List(pd.representatives, x -> 0);
        for u in [2..Length(vecs)] do
            if u^tdual = u then
                h := ys[u]*tdual/ys[u];
                j := classPosition(h, p, pd, pkeys, pcache);
                weights[j] := weights[j]+E(3)^IntFFE(vecs[u]*translation);
            fi;
        od;
        Add(columns, Concatenation(List(qd.irreducibles, r -> r[qp]),
                                   List(pd.irreducibles, r -> r*weights)));
    od;
    rows := TransposedMat(columns);
    d := rec(n := n, schema := 1, order := Size(g),
        convention := "projective qutrit Clifford = F3^(2n):Sp(2n,3)",
        generators := GeneratorsOfGroup(g),
        representatives := List(classes, Representative),
        classSizes := List(classes, Size),
        classOrders := List(classes, c -> Order(Representative(c))),
        irreducibles := rows, degrees := List(rows, r -> r[1]),
        quotientRows := Length(qd.irreducibles), auxiliary := rec(Q := qd, P := pd));
    # Only this structural table receives the constructed rows. The generic
    # comparison in verify.g uses a separate, fresh permutation group.
    t := CharacterTable(g); SetIrr(t, List(rows, r -> Character(t, r)));
    d.powerPrimes := Set(FactorsInt(d.order));
    d.powerMaps := List(d.powerPrimes, p -> PowerMap(t, p));
    QutritCheckData(d);
    frozen := ReadAsFunction(QutritPath("data/qutrit_n2.g"))();
    QutritCompareCharacterRows(d, frozen);
    output := QutritPath("generated/qutrit_structural_n2.g");
    PrintTo(output, "return ", d, ";\n");
    Print("PASS structural computation agrees with Appendix B after class matching\n");
    Print("SAVED ", output, "; GAP CPU milliseconds=", Runtime()-start, "\n");
end;
QutritComputeStructuralN2();
QUIT_GAP(0);
