# Concrete affine symplectic permutation action; no asserted group-size cache.
QutritAffineGroup := function(n)
    local f, m, vecs, powers, index, linear, gens, b;
    if not IsInt(n) or n < 1 then Error("rank must be a positive integer"); fi;
    f := GF(3); m := 2*n;
    vecs := Tuples(List([0..2], i -> i*One(f)), m);
    powers := List([1..m], i -> 3^(m-i));
    index := v -> 1 + List(v, IntFFE)*powers;
    linear := Sp(m, 3);
    gens := List(GeneratorsOfGroup(linear),
                 a -> PermList(List(vecs, v -> index(v*a))));
    for b in Basis(f^m) do
        Add(gens, PermList(List(vecs, v -> index(v+b))));
    od;
    return Group(gens);
end;
