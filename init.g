# Exact character data for the two-qutrit projective Clifford group.
# Run from the repository root, or set QutritTablesRoot before Read().
if not IsBound(QutritTablesRoot) then
    QutritTablesRoot := ".";
fi;

QutritPath := function(relative)
    return Filename(Directory(QutritTablesRoot), relative);
end;

QutritExpectedOrder := function(n)
    if not IsInt(n) or n < 1 then Error("rank must be a positive integer"); fi;
    return 3^(n*n + 2*n) * Product([1..n], i -> 3^(2*i)-1);
end;

# This function checks data algebraically, without installing cached attributes.
QutritCheckData := function(d)
    local k, rows, sizes, j, p, weighted, adjoint, columns, expected;
    if d.n <> 2 or d.schema <> 1 or d.order <> QutritExpectedOrder(2) then
        Error("unexpected rank, schema, or group order");
    fi;
    rows := d.irreducibles; sizes := d.classSizes; k := Length(rows);
    if k <> 58 or Length(d.representatives) <> k or Length(sizes) <> k
       or Length(d.classOrders) <> k or Length(d.degrees) <> k then
        Error("expected 58 rows, representatives, and class metadata entries");
    fi;
    if not ForAll(rows, r -> IsList(r) and Length(r) = k
                  and ForAll(r, IsCyclotomic)) then
        Error("character matrix must be square and cyclotomic");
    fi;
    if not ForAll(sizes, s -> IsInt(s) and s > 0 and d.order mod s = 0)
       or not ForAll(d.classOrders, o -> IsInt(o) and o > 0)
       or not ForAll(d.degrees, x -> IsInt(x) and x > 0) then
        Error("invalid class metadata or degrees");
    fi;
    if d.representatives[1] <> () or sizes[1] <> 1
       or d.classOrders[1] <> 1 or List(rows, r -> r[1]) <> d.degrees
       or not List([1..k], j -> 1) in rows then
        Error("identity column, degrees, or trivial character is invalid");
    fi;
    if Sum(sizes) <> d.order or Sum(d.degrees, x -> x^2) <> d.order then
        Error("class-size sum or degree-square sum failed");
    fi;
    weighted := List(rows, r -> List([1..k], j -> r[j]*sizes[j]));
    adjoint := TransposedMat(List(rows, r -> List(r, ComplexConjugate)));
    if weighted*adjoint <> d.order*IdentityMat(k) then
        Error("exact row orthogonality failed");
    fi;
    columns := TransposedMat(rows)*List(rows, r -> List(r, ComplexConjugate));
    expected := List([1..k], i -> List([1..k], j -> 0));
    for j in [1..k] do expected[j][j] := d.order/sizes[j]; od;
    if columns <> expected then Error("exact column orthogonality failed"); fi;
    if d.powerPrimes <> Set(FactorsInt(d.order))
       or Length(d.powerMaps) <> Length(d.powerPrimes) then
        Error("expected a power map for every prime dividing the group order");
    fi;
    for p in d.powerMaps do
        if Length(p) <> k or not ForAll(p, j -> IsInt(j) and j >= 1 and j <= k) then
            Error("invalid power-map indices");
        fi;
    od;
    return true;
end;

# Locate elements in a concrete ordered class list; cycle structure is only
# a filter. IsConjugate performs the actual group-theoretic comparison.
QutritClassPositions := function(g, representatives, elements)
    local keys, result, x, candidates, pos, j;
    keys := List(representatives, CycleStructurePerm); result := [];
    for x in elements do
        candidates := Positions(keys, CycleStructurePerm(x)); pos := fail;
        for j in candidates do
            if IsConjugate(g, representatives[j], x) then pos := j; break; fi;
        od;
        if pos = fail then Error("no matching conjugacy class"); fi;
        Add(result, pos);
    od;
    return result;
end;

QutritTable := function(n)
    local d, g, classes, j, t;
    if n <> 2 then Error("only rank 2 data are included"); fi;
    d := ReadAsFunction(QutritPath("data/qutrit_n2.g"))();
    QutritCheckData(d);
    g := Group(d.generators);
    classes := List(d.representatives, x -> ConjugacyClass(g, x));
    # Loader only: the deep validator below independently recalculates these.
    for j in [1..Length(classes)] do SetSize(classes[j], d.classSizes[j]); od;
    SetConjugacyClasses(g, classes); t := CharacterTable(g);
    SetOrdersClassRepresentatives(t, d.classOrders);
    SetIrr(t, List(d.irreducibles, r -> Character(t, r)));
    for j in [1..Length(d.powerPrimes)] do
        ComputedPowerMaps(t)[d.powerPrimes[j]] := d.powerMaps[j];
    od;
    d.group := g; d.table := t;
    return d;
end;

# deep=true recomputes on a FRESH group: stored class sizes and Irr are never
# installed on that group. In particular, the order/centralizer checks below
# are computations, not comparisons against SetSize caches.
QutritValidate := function(d, deep)
    local g, k, j, ell, keys, actual, p, ci, classes, rows, map;
    QutritCheckData(d);
    if deep then
        g := Group(d.generators); k := Length(d.representatives);
        if Size(g) <> d.order then Error("independent group order failed"); fi;
        if not ForAll(d.representatives, x -> x in g) then
            Error("a stored representative is outside the group");
        fi;
        keys := List(d.representatives, CycleStructurePerm);
        for j in [1..k] do
            if Order(d.representatives[j]) <> d.classOrders[j] then
                Error("actual representative order failed");
            fi;
            if Size(g)/Size(Centralizer(g, d.representatives[j])) <> d.classSizes[j] then
                Error("actual centralizer/class size failed");
            fi;
            for ell in [1..j-1] do
                if keys[j] = keys[ell]
                   and IsConjugate(g, d.representatives[j], d.representatives[ell]) then
                    Error("stored classes are not distinct");
                fi;
            od;
        od;
        Print("PASS actual group order, class sizes, distinctness, and element orders\n");
        for ci in [1..Length(d.powerPrimes)] do
            p := d.powerPrimes[ci];
            actual := QutritClassPositions(g, d.representatives,
                                           List(d.representatives, x -> x^p));
            if actual <> d.powerMaps[ci] then Error("actual prime power map failed"); fi;
        od;
        Print("PASS actual power maps for primes ", d.powerPrimes, "\n");
        classes := ConjugacyClasses(g);
        rows := List(Irr(CharacterTable(g)), ValuesOfClassFunction);
        map := QutritClassPositions(g, List(classes, Representative), d.representatives);
        rows := List(rows, r -> List(map, j -> r[j]));
        if Set(rows) <> Set(d.irreducibles) then
            Error("fresh generic irreducible characters do not match stored rows");
        fi;
        Print("PASS fresh generic GAP character calculation after class matching\n");
    fi;
    Print("PASS rank=2, 58 classes, exact row/column orthogonality, deep=", deep, "\n");
    return true;
end;

QutritCompareCharacterRows := function(a, b)
    local g, h, map, rows;
    QutritCheckData(a); QutritCheckData(b);
    g := Group(a.generators); h := Group(b.generators);
    if Size(g) <> Size(h) or not ForAll(b.generators, x -> x in g) then
        Error("the concrete permutation groups differ");
    fi;
    map := QutritClassPositions(g, a.representatives, b.representatives);
    if Length(Set(map)) <> Length(map) then Error("class matching is not bijective"); fi;
    rows := List(a.irreducibles, r -> List(map, j -> r[j]));
    if Set(rows) <> Set(b.irreducibles) then Error("complete character row sets differ"); fi;
    return true;
end;
