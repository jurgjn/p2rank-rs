#!/usr/bin/env python3
"""Extract the chem-comp fields P2Rank's BioJava parsing depends on from mmCIF chem comp files.

Output: one line per component:
  ID <TAB> type <TAB> one_letter_code <TAB> mon_nstd_parent_comp_id <TAB> three_letter_code <TAB> atom_id=type_symbol;...
Usage: gen_chemcomp.py OUT FILE.cif[.gz]...   (a file may hold many data_ blocks, e.g. components.cif.gz)
"""
import gzip, sys


def tokens(line):
    i, n = 0, len(line)
    while i < n:
        c = line[i]
        if c.isspace():
            i += 1
            continue
        if c in "'\"" :
            j = i + 1
            while j < n:
                if line[j] == c and (j + 1 == n or line[j + 1].isspace()):
                    break
                j += 1
            yield line[i + 1:j]
            i = j + 1
        else:
            j = i
            while j < n and not line[j].isspace():
                j += 1
            yield line[i:j]
            i = j


def blocks(path):
    op = gzip.open if path.endswith(".gz") else open
    cur = None
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("data_"):
                if cur is not None:
                    yield cur
                cur = []
            if cur is not None:
                cur.append(line)
    if cur is not None:
        yield cur


def parse_block(lines):
    items = {}
    loops = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()
        if s == "loop_":
            cols = []
            i += 1
            while i < n and lines[i].strip().startswith("_"):
                cols.append(lines[i].strip())
                i += 1
            vals = []
            while i < n:
                s2 = lines[i].strip()
                if s2 == "#" or s2 == "loop_" or (s2.startswith("_") ) or lines[i].startswith("data_"):
                    break
                if lines[i].startswith(";"):
                    txt = [lines[i][1:]]
                    i += 1
                    while i < n and not lines[i].startswith(";"):
                        txt.append(lines[i])
                        i += 1
                    vals.append("\n".join(txt))
                    i += 1
                    continue
                vals.extend(tokens(lines[i]))
                i += 1
            rows = [vals[k:k + len(cols)] for k in range(0, len(vals) - len(cols) + 1, len(cols))] if cols else []
            loops.append((cols, rows))
            continue
        if s.startswith("_"):
            toks = list(tokens(s))
            key = toks[0]
            if len(toks) >= 2:
                items[key] = toks[1]
                i += 1
            else:
                i += 1
                if i < n and lines[i].startswith(";"):
                    txt = [lines[i][1:]]
                    i += 1
                    while i < n and not lines[i].startswith(";"):
                        txt.append(lines[i])
                        i += 1
                    items[key] = "\n".join(txt)
                    i += 1
                elif i < n:
                    t = list(tokens(lines[i]))
                    items[key] = t[0] if t else ""
                    i += 1
            continue
        i += 1
    atoms = []
    for cols, rows in loops:
        if "_chem_comp_atom.atom_id" in cols:
            ia = cols.index("_chem_comp_atom.atom_id")
            it = cols.index("_chem_comp_atom.type_symbol")
            for r in rows:
                atoms.append((r[ia], r[it]))
    return items, atoms


def clean(v):
    return (v or "").replace("\t", " ").replace("\n", " ")


def main():
    out = open(sys.argv[1], "w")
    for path in sys.argv[2:]:
        for b in blocks(path):
            items, atoms = parse_block(b)
            cid = items.get("_chem_comp.id")
            if cid is None:
                continue
            for a, t in atoms:
                if "=" in a or ";" in a:
                    raise SystemExit("bad atom id %r in %s" % (a, cid))
            out.write("\t".join([clean(cid), clean(items.get("_chem_comp.type", "")),
                                 clean(items.get("_chem_comp.one_letter_code", "")),
                                 clean(items.get("_chem_comp.mon_nstd_parent_comp_id", "")),
                                 clean(items.get("_chem_comp.three_letter_code", "")),
                                 ";".join("%s=%s" % (a, t) for a, t in atoms)]) + "\n")


main()
