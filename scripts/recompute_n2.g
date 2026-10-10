# Usage from repository root:
# gap -q -A -b --quitonbreak scripts/recompute_n2.g
if not IsBound(QutritTablesRoot) then QutritTablesRoot := "."; fi;
Read(Filename(Directory(QutritTablesRoot), "init.g"));
Read(QutritPath("scripts/affine.g"));
QutritRecomputeN2 := function()
    local start, g, classes, t, rows, d, frozen, output;
    start := Runtime();
    Print("DIRECT_START n=2 GAP=", GAPInfo.Version, "\n");
    g := QutritAffineGroup(2);
    if Size(g) <> QutritExpectedOrder(2) then Error("computed group order failed"); fi;
    classes := ConjugacyClasses(g); t := CharacterTable(g);
    rows := List(Irr(t), ValuesOfClassFunction);
    d := rec(n := 2, schema := 1, order := Size(g),
        convention := "projective qutrit Clifford = F3^(2n):Sp(2n,3)",
        generators := GeneratorsOfGroup(g),
        representatives := List(classes, Representative),
        classSizes := List(classes, Size),
        classOrders := List(classes, c -> Order(Representative(c))),
        irreducibles := rows, degrees := List(rows, r -> r[1]));
    d.powerPrimes := Set(FactorsInt(d.order));
    d.powerMaps := List(d.powerPrimes, p -> PowerMap(t, p));
    QutritCheckData(d);
    frozen := ReadAsFunction(QutritPath("data/qutrit_n2.g"))();
    QutritCompareCharacterRows(d, frozen);
    output := QutritPath("generated/qutrit_n2_recomputed.g");
    PrintTo(output, "return ", d, ";\n");
    Print("PASS direct recomputation agrees with Appendix B after class matching\n");
    Print("SAVED ", output, "; GAP CPU milliseconds=", Runtime()-start, "\n");
end;
QutritRecomputeN2();
QUIT_GAP(0);
