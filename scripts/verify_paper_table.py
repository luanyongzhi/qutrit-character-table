#!/usr/bin/env python3
"""Check the frozen n=2 data and optionally compare the paper's Appendix B.

Use --data-only without a manuscript, or supply the path to its LaTeX source.
Python standard library only. The manuscript is not bundled with the repository.
This is a narrow parser for the paper's six LaTeX longtable panels and frozen
GAP data expressions, not a general TeX or GAP interpreter.
"""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def multiply(x, y):
    """(a,b) represents a+b*zeta with zeta^2+zeta+1=0."""
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c-b*d


def inner_product(x, y, weights):
    result_a = result_b = 0
    for (a, b), (c, d), weight in zip(x, y, weights):
        real, coeff = multiply((a, b), (c-d, -d))
        result_a += weight*real
        result_b += weight*coeff
    return result_a, result_b


def parse_entry(text):
    text = text.strip().strip('$').strip()
    if re.fullmatch(r'-?\d+', text):
        return int(text), 0
    match = re.fullmatch(r'\\langle\s*(-?\d+)\s*,\s*(-?\d+)\s*\\rangle', text)
    require(match is not None, 'Unrecognized LaTeX entry: '+repr(text))
    return tuple(map(int, match.groups()))


def parse_paper(text, count):
    text = re.sub(r'(?<!\\)%[^\n]*', '', text)
    marker = r'\label{app:character-table-n2}'
    require(text.count(marker) == 1, 'Expected one active Appendix B label: '+marker)
    section = text.split(marker, 1)[1]
    marker = r'\subsection*{Class metadata}'
    require(marker in section, 'Class metadata subsection not found')
    metadata = section.split(marker, 1)[1].split(r'\end{longtable}', 1)[0]
    classes = {}
    for line in metadata.splitlines():
        if not re.match(r'^\s*\d+\s*&', line):
            continue
        values = [int(value.strip().rstrip('\\')) for value in line.split('&')]
        require(len(values) == 8, 'Expected two four-field class records: '+line)
        for offset in (0, 4):
            index, order, size, centralizer = values[offset:offset+4]
            require(index not in classes, 'Duplicate class index '+str(index))
            classes[index] = (order, size, centralizer)
    require(sorted(classes) == list(range(1, count+1)), 'Incomplete class metadata')

    table = [[None]*count for _ in range(count)]
    panels = 0
    for block in section.split(r'\begin{longtable}')[1:]:
        block = block.split(r'\end{longtable}', 1)[0]
        headers = [line for line in block.splitlines() if '$K_' in line]
        if not headers:
            continue
        columns = [int(value) for value in re.findall(r'\$K_\{(\d+)\}\$', headers[0])]
        require(columns and all(1 <= j <= count for j in columns), 'Invalid panel columns')
        rows = []
        for line in block.splitlines():
            match = re.match(r'\s*\$\\chi_\{(\d+)\}\$\s*&\s*(.*?)\\\\\s*$', line)
            if not match:
                continue
            row = int(match[1])
            require(1 <= row <= count, 'Row index out of range')
            rows.append(row)
            values = [parse_entry(value) for value in match[2].split('&')]
            require(len(values) == len(columns), 'Wrong entry count in row '+str(row))
            for column, value in zip(columns, values):
                require(table[row-1][column-1] is None, 'Duplicate table cell')
                table[row-1][column-1] = value
        require(sorted(rows) == list(range(1, count+1)), 'Incomplete or duplicate panel rows')
        panels += 1
    require(panels == 6, 'Expected six character-value panels, found '+str(panels))
    require(all(value is not None for row in table for value in row), 'Incomplete value matrix')
    return classes, table, panels


def gap_list(source, field):
    """Extract one bracketed data field; never execute the GAP source."""
    match = re.search(r'\b'+re.escape(field)+r'\s*:=\s*\[', source)
    require(match is not None, 'Missing GAP list field '+field)
    start = match.end()-1
    depth = 0
    for end in range(start, len(source)):
        if source[end] == '[':
            depth += 1
        elif source[end] == ']':
            depth -= 1
        if depth == 0:
            return source[start:end+1]
    raise ValueError('Unclosed GAP list field '+field)


def cyclotomic_literal(node):
    """Whitelist only lists, integers, zeta, +, -, *, and squares."""
    if isinstance(node, ast.List):
        return [cyclotomic_literal(child) for child in node.elts]
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return node.value, 0
    if isinstance(node, ast.Name) and node.id == 'zeta':
        return 0, 1
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        a, b = cyclotomic_literal(node.operand)
        return -a, -b
    if isinstance(node, ast.BinOp):
        a, b = cyclotomic_literal(node.left)
        c, d = cyclotomic_literal(node.right)
        if isinstance(node.op, ast.Add):
            return a+c, b+d
        if isinstance(node.op, ast.Sub):
            return a-c, b-d
        if isinstance(node.op, ast.Mult):
            return multiply((a, b), (c, d))
        if isinstance(node.op, ast.Pow) and (c, d) == (2, 0):
            return multiply((a, b), (a, b))
    raise ValueError('Unsupported GAP cyclotomic expression: '+ast.dump(node))


def verify(args):
    data = json.loads(args.data.read_text(encoding='utf-8'))
    require(isinstance(data, dict), 'Expected a JSON object')
    require(type(data['classes']) is int and type(data['order']) is int
            and data['classes'] == 58 and data['order'] == 4199040,
            'Expected the frozen 58-class two-qutrit table')
    count, order = data['classes'], data['order']
    for field in ('class_orders', 'class_sizes', 'centralizer_orders', 'degrees'):
        values = data[field]
        require(isinstance(values, list) and len(values) == count,
                'Expected 58 entries in JSON '+field)
        require(all(type(value) is int and value > 0 for value in values),
                'Expected positive integers in JSON '+field)
    require(isinstance(data['values'], list) and len(data['values']) == count
            and all(isinstance(row, list) and len(row) == count for row in data['values']),
            'Invalid JSON value-matrix dimensions')
    require(all(isinstance(value, list) and len(value) == 2
                and all(type(coefficient) is int for coefficient in value)
                for row in data['values'] for value in row),
            'Expected integer pairs for all JSON character values')
    table = [[tuple(value) for value in row] for row in data['values']]
    if not args.data_only:
        classes, paper_table, panels = parse_paper(args.manuscript.read_text(encoding='utf-8'), count)
        require(paper_table == table,
                'Paper values differ from JSON in their saved row/column order')
        for offset, field in enumerate(('class_orders', 'class_sizes', 'centralizer_orders')):
            require([classes[j][offset] for j in range(1, count+1)] == data[field],
                    'Paper metadata differs from JSON: '+field)

    gap_source = args.gap_data.read_text(encoding='utf-8')
    for field, key in (('classOrders', 'class_orders'), ('classSizes', 'class_sizes'),
                       ('degrees', 'degrees')):
        require(ast.literal_eval(gap_list(gap_source, field)) == data[key],
                'GAP data differs from JSON: '+field)
    expression = gap_list(gap_source, 'irreducibles').replace('E(3)', 'zeta').replace('^', '**')
    require(cyclotomic_literal(ast.parse(expression, mode='eval').body) == table,
            'GAP character values differ from JSON')
    order_match = re.search(r'\border\s*:=\s*(\d+)\s*,', gap_source)
    require(order_match is not None and int(order_match[1]) == order, 'GAP group order differs')

    sizes, centralizers = data['class_sizes'], data['centralizer_orders']
    require(data['class_orders'][0] == sizes[0] == 1 and centralizers[0] == order,
            'First class is not the identity class')
    require(data['class_orders'].count(1) == 1, 'More than one class has element order 1')
    require(sum(sizes) == order, 'Class sizes do not sum to the group order')
    require(all(s*c == order for s, c in zip(sizes, centralizers)), 'Wrong centralizer orders')
    require(all(c % o == 0 for o, c in zip(data['class_orders'], centralizers)),
            'Invalid element-order divisibility')
    require(all(row[0] == (degree, 0) for row, degree in zip(table, data['degrees'])),
            'First column is not the positive degree list')
    require(sum(degree**2 for degree in data['degrees']) == order, 'Wrong degree-square sum')
    require(all(value == (1, 0) for value in table[0]), 'First row is not the trivial character')
    columns = list(zip(*table))
    for i in range(count):
        for j in range(count):
            require(inner_product(table[i], table[j], sizes) == (order if i == j else 0, 0),
                    'Row orthogonality failed at '+str((i+1, j+1)))
            require(inner_product(columns[i], columns[j], [1]*count) == (centralizers[i] if i == j else 0, 0),
                    'Column orthogonality failed at '+str((i+1, j+1)))
    report = {
        'status': 'PASS', 'mode': 'data-only' if args.data_only else 'paper',
        'matrix_shape': [count, count], 'group_order': order,
        'json_gap_same_order': True, 'metadata_match': True,
        'row_pairs_checked': count**2, 'column_pairs_checked': count**2,
        'arithmetic': 'exact integer pairs in Z[zeta_3]',
        'scope': 'Saved value matrices and metadata, plus the printed table if supplied. No group reconstruction '
                 'or check of actual class representatives/power maps in this Python script.'
    }
    if not args.data_only:
        report.update({
            'manuscript': args.manuscript.name,
            'manuscript_sha256': hashlib.sha256(args.manuscript.read_bytes()).hexdigest(),
            'panels': panels, 'paper_json_gap_same_order': True
        })
    return report


def main():
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manuscript', nargs='?', type=Path, help='Path to the paper LaTeX source')
    parser.add_argument('--data-only', action='store_true', help='Check frozen GAP/JSON data without a manuscript')
    parser.add_argument('--data', type=Path, default=root/'data'/'qutrit_n2.json')
    parser.add_argument('--gap-data', type=Path, default=root/'data'/'qutrit_n2.g')
    parser.add_argument('--report', type=Path, help='Optional JSON report to write')
    args = parser.parse_args()
    if args.manuscript is None and not args.data_only:
        parser.error('supply a manuscript path or use --data-only')
    if args.manuscript is not None and args.data_only:
        parser.error('supply either a manuscript path or --data-only, not both')
    try:
        report = verify(args)
        result = json.dumps(report, indent=2)+'\n'
        if args.report:
            args.report.write_text(result, encoding='utf-8')
        print(result, end='')
    except (OSError, ValueError, KeyError, TypeError, SyntaxError) as error:
        print('FAIL: '+str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
