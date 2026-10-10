# Usage and data

The data describe the projective two-qutrit Clifford group
`G = F_3^4 : Sp(4,3)`, with scalar phases quotiented out. Its order is
`4199040`; its action has 81 points and its character table has 58 rows
and 58 columns.

## Load the table

Use GAP 4.15.1. Start GAP from this directory, then enter:

```gap
Read("init.g");;
D := QutritTable(2);;
Display(D.table);
D.degrees;
D.classOrders;
D.classSizes;
PowerMap(D.table, 3);
```

If GAP is started elsewhere, set the directory before reading `init.g`:

```gap
QutritTablesRoot := "/path/to/qutrit_character_table";;
Read(Concatenation(QutritTablesRoot, "/init.g"));;
D := QutritTable(2);;
```

## Meaning of the data

`data/qutrit_n2.g` preserves the row and column order of the printed table.
Rows are ordinary irreducible characters; columns are conjugacy classes.
The first column is the identity class and gives the character degrees.
Indices start at 1.

| Record field | Meaning |
|---|---|
| `n` | Number of qutrits, 2 |
| `order` | Group order |
| `generators` | Permutation generators |
| `representatives` | Conjugacy-class representatives, in column order |
| `classSizes` | Conjugacy-class sizes |
| `classOrders` | Orders of the representatives |
| `irreducibles` | Exact character-value matrix |
| `degrees` | Character degrees |
| `powerPrimes` | Primes `[2,3,5]` |
| `powerMaps` | For each prime `p`, the class index containing each representative's `p`th power |
| `group`, `table` | GAP objects added by `QutritTable(2)` |

The centralizer order for column `j` is `D.order / D.classSizes[j]`.

`data/qutrit_n2.json` contains the same character values and class metadata.
Each pair `[a,b]` means `a + b*zeta`, where `zeta = exp(2*pi*i/3)` and
`zeta^2 + zeta + 1 = 0`. In GAP, `zeta` is `E(3)`. These values are exact.

`data/qutrit_structural_n2.g` contains an independent table with its own
ordering. Its first 34 rows are inflated from `Sp(4,3)`; the remaining
24 are induced from the stabilizer of a nonzero dual vector. The fields
`auxiliary.Q` and `auxiliary.P` contain those two auxiliary character tables.
Use `QutritCompareCharacterRows(a,b)` to compare records by actual class
matching, rather than matching indices or degrees. These files contain
characters, not matrices for all irreducible representations.

## Verify and recompute

Run these commands from this directory. Replace `gap` by the path to your
executable if necessary, such as `~/gap-4.15.1/gap`.

```bash
gap -q -A -b --quitonbreak verify.g
gap -q -A -b --quitonbreak scripts/recompute_n2.g
gap -q -A -b --quitonbreak scripts/compute_structural.g
```

- `verify.g` checks the concrete group, classes, centralizers, prime power
  maps, orthogonality and agreement with freshly computed characters.
- `recompute_n2.g` uses the generic character-table calculation.
- `compute_structural.g` uses the quotient and stabilizer character tables.

Recomputed tables are written to `generated/`; the files in `data/` are
unchanged. New computations may return different row and column orders.
The scripts compare complete rows after matching conjugacy classes.

For exact checks of the three-copy identities at `n=1,2` and the Jacobian:

```bash
gap -q -A -b --quitonbreak scripts/check_identities.g
```

For independent integer-arithmetic checks of the GAP and JSON tables,
use Python 3.9 or later; no additional packages are required:

```bash
python3 scripts/verify_paper_table.py --data-only
```

To also compare the printed table in a LaTeX source:

```bash
python3 scripts/verify_paper_table.py /path/to/paper.tex
```

The latter command expects the six character-table panels in Appendix B.
All checks report failure with a nonzero exit code. The GAP command-line
scripts terminate GAP when finished; `init.g` only loads the interface.
