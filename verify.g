# Usage: gap -q -A -b --quitonbreak verify.g
if not IsBound(QutritTablesRoot) then QutritTablesRoot := "."; fi;
Read(Filename(Directory(QutritTablesRoot), "init.g"));
QutritRunVerification := function()
    local d, structural;
    Print("GAP ", GAPInfo.Version, "; exact character-table verification\n");
    d := QutritTable(2);
    QutritValidate(d, true);
    structural := ReadAsFunction(QutritPath("data/qutrit_structural_n2.g"))();
    QutritCompareCharacterRows(structural, d);
    Print("PASS stored little-group computation agrees after actual class matching\n");
    Print("ALL_CHECKS_PASS\n");
end;
QutritRunVerification();
QUIT_GAP(0);
