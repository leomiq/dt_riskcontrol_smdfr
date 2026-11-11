"""Convert docs/monolithic_formulation_plan.md to LaTeX (tailored converter, no pandoc)."""
import re, sys

src, dst = sys.argv[1], sys.argv[2]
lines = open(src, encoding="utf8").read().splitlines()

SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6",
       "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+"}
TXTMATH = {"≈": r"\approx", "≤": r"\le", "≥": r"\ge", "≪": r"\ll", "→": r"\to",
           "↔": r"\leftrightarrow", "×": r"\times", "±": r"\pm", "−": "-", "∩": r"\cap",
           "✓": r"\checkmark", "✗": r"\times", "λ": r"\lambda", "Δ": r"\Delta",
           "α": r"\alpha", "Λ": r"\Lambda", "β": r"\beta", "η": r"\eta",
           "γ": r"\gamma", "ρ": r"\rho", "°": r"^\circ"}


def esc_text(s):
    s = s.replace("\\", "\x00")
    for a, b in [("&", r"\&"), ("%", r"\%"), ("#", r"\#"), ("_", r"\_"),
                 ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\^{}")]:
        s = s.replace(a, b)
    s = s.replace("\x00", r"\textbackslash{}")
    s = re.sub("[" + "".join(SUP) + "]+",
               lambda m: "$^{" + "".join(SUP[c] for c in m.group(0)) + "}$", s)
    s = "".join("$" + TXTMATH[c] + "$" if c in TXTMATH else c for c in s)
    s = s.replace("…", r"\ldots{}").replace("§", r"\S{}")
    return s


def esc_code(s):
    out = []
    for c in s:
        out.append({"\\": r"\textbackslash{}", "_": r"\_", "#": r"\#", "%": r"\%",
                    "&": r"\&", "~": r"\textasciitilde{}", "^": r"\^{}", "{": r"\{",
                    "}": r"\}", "$": r"\$"}.get(c, c))
    return r"\code{" + "".join(out) + "}"


def fix_math(m):
    # unicode inside math
    for c, r in TXTMATH.items():
        m = m.replace(c, r + " ")
    return m.replace("·", r"\cdot ").replace("…", r"\ldots ")


def inline(s):
    """Inline markup: code, math, links, bold, italic."""
    tokens = []

    def keep(t):
        tokens.append(t)
        return f"\x01{len(tokens)-1}\x02"

    s = re.sub(r"`([^`]+)`", lambda m: keep(esc_code(m.group(1))), s)
    s = re.sub(r"(?<!\\)\$([^$]+?)\$", lambda m: keep("$" + fix_math(m.group(1)) + "$"), s)
    # Unicode Greek letter with a plain-text subscript, e.g. "β_i" -> $\beta_i$
    s = re.sub(r"([λΔαΛβηγρ])_(\w+)", lambda m: keep("$" + fix_math(m.group(1)) + "_{" + m.group(2) + "}$"), s)
    links = []

    def link(m):
        links.append((m.group(1), m.group(2)))
        return f"\x03{len(links)-1}\x04"

    s = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link, s)
    s = esc_text(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", s)

    def unlink(m):
        t, u = links[int(m.group(1))]
        u = u.replace("%", r"\%").replace("#", r"\#")
        return r"\href{" + u + "}{" + re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", esc_text(t)) + "}"

    s = re.sub("\x03(\\d+)\x04", unlink, s)
    while "\x01" in s:
        s = re.sub("\x01(\\d+)\x02", lambda m: tokens[int(m.group(1))], s)
    return s


def split_row(r):
    r = r.strip()
    if r.startswith("|"):
        r = r[1:]
    if r.endswith("|"):
        r = r[:-1]
    cells, cur, inm, inc = [], "", False, False
    i = 0
    while i < len(r):
        c = r[i]
        if c == "`" and not inm:
            inc = not inc
        elif c == "$" and not inc and (i == 0 or r[i-1] != "\\"):
            inm = not inm
        if c == "|" and not inm and not inc:
            cells.append(cur.strip()); cur = ""
        else:
            cur += c
        i += 1
    cells.append(cur.strip())
    return cells


def table(rows):
    head = split_row(rows[0])
    align = split_row(rows[1])
    body = [split_row(r) for r in rows[2:]]
    n = len(head)
    spec = []
    for a in align:
        a = a.strip()
        spec.append("R" if a.endswith(":") and not a.startswith(":") else "Y")
    size = r"\footnotesize" if n >= 4 else r"\small"
    out = [r"{" + size, r"\begin{xltabular}{\linewidth}{@{}" + "".join(spec) + "@{}}", r"\toprule",
           " & ".join(r"\textbf{" + inline(h) + "}" for h in head) + r" \\", r"\midrule",
           r"\endhead"]
    for b in body:
        b = (b + [""] * n)[:n]
        out.append(" & ".join(inline(c) for c in b) + r" \\")
    out += [r"\bottomrule", r"\end{xltabular}}", ""]
    return out


def heading(level, text):
    if level == 1:
        return None
    cmd = {2: "section", 3: "subsection", 4: "subsubsection"}[level]
    t = inline(text)
    return [rf"\{cmd}*{{{t}}}", rf"\addcontentsline{{toc}}{{{cmd}}}{{{t}}}"]


out, title = [], ""
i = 0
para = []


def flush():
    global para
    if para:
        out.append(inline(" ".join(para)))
        out.append("")
        para = []


list_re = re.compile(r"^(\s*)([-*]|\d+\.)\s+(.*)$")
while i < len(lines):
    L = lines[i]
    if L.startswith("# ") and not title:
        title = L[2:].strip(); i += 1; continue
    m = re.match(r"^(#{2,4})\s+(.*)$", L)
    if m:
        flush(); out += heading(len(m.group(1)), m.group(2)); i += 1; continue
    if L.strip() == "---":
        flush(); out.append(r"\medskip"); i += 1; continue
    if L.strip().startswith("```"):
        flush(); i += 1; buf = []
        while not lines[i].strip().startswith("```"):
            buf.append(lines[i]); i += 1
        i += 1
        out += [r"{\footnotesize", r"\begin{verbatim}"] + buf + [r"\end{verbatim}}", ""]
        continue
    if L.strip() == "$$":
        flush(); i += 1; buf = []
        while lines[i].strip() != "$$":
            buf.append(lines[i]); i += 1
        i += 1
        out += [r"\[", fix_math("\n".join(buf)), r"\]", ""]
        continue
    if L.startswith("|"):
        flush(); rows = []
        while i < len(lines) and lines[i].startswith("|"):
            rows.append(lines[i]); i += 1
        out += table(rows); continue
    m = list_re.match(L)
    if m:
        flush()
        # collect list block (items + continuation lines + nested)
        stack = []  # (indent, env)
        while i < len(lines):
            L = lines[i]
            m = list_re.match(L)
            if m:
                ind = len(m.group(1))
                env = "enumerate" if m.group(2)[0].isdigit() else "itemize"
                while stack and stack[-1][0] > ind:
                    out.append(rf"\end{{{stack.pop()[1]}}}")
                if not stack or stack[-1][0] < ind:
                    stack.append((ind, env)); out.append(rf"\begin{{{env}}}")
                item = [m.group(3)]
                i += 1
                while i < len(lines) and lines[i].strip() and not list_re.match(lines[i]) \
                        and lines[i].startswith(" ") and not lines[i].strip() == "$$":
                    item.append(lines[i].strip()); i += 1
                # display math inside item
                while i < len(lines) and lines[i].strip() == "$$":
                    i += 1; buf = []
                    while lines[i].strip() != "$$":
                        buf.append(lines[i]); i += 1
                    i += 1
                    item.append("\x05" + fix_math("\n".join(buf)) + "\x06")
                txt = " ".join(item)
                parts = re.split("\x05|\x06", txt)
                res = ""
                for k, p in enumerate(parts):
                    res += (r"\[" + p + r"\]") if k % 2 else inline(p)
                out.append(r"\item " + res)
                continue
            if not L.strip():
                # blank: list continues only if next non-blank is a list item
                j = i
                while j < len(lines) and not lines[j].strip():
                    j += 1
                if j < len(lines) and list_re.match(lines[j]) and (lines[j].startswith(" ") or stack[0][0] == 0):
                    i = j; continue
            break
        while stack:
            out.append(rf"\end{{{stack.pop()[1]}}}")
        out.append("")
        continue
    if not L.strip():
        flush(); i += 1; continue
    para.append(L.strip()); i += 1
flush()

PRE = r"""% Generado desde docs/monolithic_formulation_plan.md con tools/md2tex_plan.py (no editar a mano).
% Regenerar: python tools/md2tex_plan.py docs/monolithic_formulation_plan.md docs/monolithic_formulation_plan.tex
% Compilar desde docs/ con Tectonic:  tectonic -X compile monolithic_formulation_plan.tex
\documentclass[11pt,a4paper]{article}
\usepackage[spanish,es-noquoting,es-nodecimaldot]{babel}
\usepackage{fontspec}
\setmainfont[Path=../assets/fonts/,BoldFont=InterVariable.ttf,
  BoldFeatures={FakeBold=2.2},ItalicFont=InterVariable.ttf,ItalicFeatures={FakeSlant=0.15},
  BoldItalicFont=InterVariable.ttf,BoldItalicFeatures={FakeBold=2.2,FakeSlant=0.15}]{InterVariable.ttf}
\usepackage[margin=22mm]{geometry}
\usepackage{amsmath,amssymb,booktabs,tabularx,xltabular,array}
\usepackage{xcolor,enumitem,titlesec}
\usepackage[skip=4pt plus 1pt]{parskip}
\usepackage{hyperref}
\definecolor{ink}{HTML}{12233B}
\definecolor{accent}{HTML}{003876}
\hypersetup{colorlinks=true,linkcolor=accent,urlcolor=accent}
\setlist{itemsep=2pt,topsep=3pt}
\renewcommand{\arraystretch}{1.2}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\newcolumntype{R}{>{\raggedleft\arraybackslash}X}
\newcommand{\codeus}{\textunderscore\allowbreak}
\newcommand{\code}[1]{{\small\ttfamily\let\_\codeus #1}}
\setlength{\emergencystretch}{3em}
\titleformat{\section}{\large\bfseries\color{accent}}{\thesection}{0.6em}{}
\titleformat{\subsection}{\normalsize\bfseries\color{ink}}{\thesubsection}{0.6em}{}
\titleformat{\subsubsection}{\normalsize\itshape\color{ink}}{\thesubsubsection}{0.6em}{}
\title{\color{accent}\bfseries TITLE}
\author{Proyecto DT--SMDFR Risk Control}
\date{Versión con la reorientación del 11 de noviembre de 2026}
\begin{document}
\maketitle
\tableofcontents
\bigskip
"""
doc = PRE.replace("TITLE", inline(title)) + "\n".join(out) + "\n\\end{document}\n"
open(dst, "w", encoding="utf8").write(doc)
print("ok", len(out))
