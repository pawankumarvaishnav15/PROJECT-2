"""

Modes (MODE key): COMP, CMPLX, STAT, BASE-N, EQN, MATRIX, TABLE, VECTOR
SHIFT+MODE = SETUP,  SHIFT+7 CONST,  SHIFT+8 CONV,  SHIFT+log▫ = Σ
SHIFT+1 STAT, SHIFT+2 CMPLX menu, SHIFT+3 BASE menu, SHIFT+4 MATRIX, SHIFT+5 VECTOR
"""
import cmath
import math
import random
import re
import tkinter as tk
from fractions import Fraction
from tkinter import messagebox, simpledialog

VARS = "ABCDEFXYM"

# ---------------- display token -> python ----------------
REPL = [
    ("sinh⁻¹(", "asinh("), ("cosh⁻¹(", "acosh("), ("tanh⁻¹(", "atanh("),
    ("sin⁻¹(", "asin("), ("cos⁻¹(", "acos("), ("tan⁻¹(", "atan("),
    ("⁻¹", "**(-1)"), ("²", "**2"), ("³", "**3"), ("^", "**"),
    ("×", "*"), ("÷", "/"), ("−", "-"),
    ("ˣ√(", "rtx("), ("√(", "sqrt("), ("∛(", "cbrt("),
    ("π", "(pi)"), ("ℯ", "(e)"), ("Ran#", "ran()"), ("Abs(", "abs("),
    ("Pol(", "pol("), ("Rec(", "rec("), ("Rnd(", "rnd("),
    ("Arg(", "arg("), ("Conjg(", "conjg("), ("Re(", "cre("), ("Im(", "cim("),
    ("det(", "det("), ("Trn(", "trn("), ("Ident(", "ident("),
    ("Dot(", "dot("), ("Angle(", "angle("), ("UnitV(", "unitv("),
    ("∫(", "integ(lambda X:"), ("d/dx(", "deriv(lambda X:"), ("Σ(", "sigma(lambda X:"),
]

CONSTS = [
    ("mp", "Proton mass", 1.672621777e-27), ("mn", "Neutron mass", 1.674927351e-27),
    ("me", "Electron mass", 9.10938291e-31), ("mμ", "Muon mass", 1.883531475e-28),
    ("a0", "Bohr radius", 5.2917721092e-11), ("h", "Planck constant", 6.62606957e-34),
    ("μN", "Nuclear magneton", 5.05078353e-27), ("μB", "Bohr magneton", 9.27400968e-24),
    ("ħ", "Reduced Planck const", 1.054571726e-34), ("α", "Fine-structure const", 7.2973525698e-3),
    ("re", "Classical electron radius", 2.8179403267e-15), ("λc", "Compton wavelength", 2.4263102389e-12),
    ("γp", "Proton gyromagnetic ratio", 2.675222005e8), ("λcp", "Proton Compton wavelength", 1.32140985623e-15),
    ("λcn", "Neutron Compton wavelength", 1.319590906e-15), ("R∞", "Rydberg constant", 10973731.568539),
    ("u", "Atomic mass unit", 1.660538921e-27), ("μp", "Proton magnetic moment", 1.410606743e-26),
    ("μe", "Electron magnetic moment", -9.28476430e-24), ("μn", "Neutron magnetic moment", -9.6623647e-27),
    ("μμ", "Muon magnetic moment", -4.49044807e-26), ("F", "Faraday constant", 96485.3365),
    ("e", "Elementary charge", 1.602176565e-19), ("NA", "Avogadro constant", 6.02214129e23),
    ("k", "Boltzmann constant", 1.3806488e-23), ("Vm", "Molar volume (ideal gas)", 0.022413968),
    ("R", "Molar gas constant", 8.3144621), ("C0", "Speed of light", 299792458.0),
    ("C1", "First radiation const", 3.74177153e-16), ("C2", "Second radiation const", 0.01438777),
    ("σ", "Stefan-Boltzmann const", 5.670373e-8), ("ε0", "Electric constant", 8.854187817e-12),
    ("μ0", "Magnetic constant", 1.2566370614e-6), ("Φ0", "Magnetic flux quantum", 2.067833758e-15),
    ("g", "Standard gravity", 9.80665), ("G0", "Conductance quantum", 7.7480917346e-5),
    ("Z0", "Impedance of vacuum", 376.730313461), ("t", "Celsius temperature", 273.15),
    ("G", "Gravitational constant", 6.67384e-11), ("atm", "Standard atmosphere", 101325.0),
]

_PAIRS = [("in", "cm", 2.54), ("ft", "m", 0.3048), ("yd", "m", 0.9144), ("mile", "km", 1.609344),
          ("n mile", "m", 1852), ("acre", "m²", 4046.8564224), ("gal(US)", "ℓ", 3.785411784),
          ("gal(UK)", "ℓ", 4.54609), ("pc", "km", 3.08567782e13), ("km/h", "m/s", 1 / 3.6),
          ("oz", "g", 28.349523125), ("lb", "kg", 0.45359237), ("atm", "Pa", 101325.0),
          ("mmHg", "Pa", 133.322368421), ("hp", "kW", 0.745699872), ("kgf/cm²", "Pa", 98066.5),
          ("kgf·m", "J", 9.80665), ("lbf/in²", "kPa", 6.8947572932)]
CONV = []
for _a, _b, _f in _PAIRS:
    CONV += [(f"{_a} ▸ {_b}", lambda x, f=_f: x * f), (f"{_b} ▸ {_a}", lambda x, f=_f: x / f)]
CONV += [("°F ▸ °C", lambda x: (x - 32) / 1.8), ("°C ▸ °F", lambda x: x * 1.8 + 32),
         ("J ▸ cal", lambda x: x / 4.18605), ("cal ▸ J", lambda x: x * 4.18605)]


# ---------------- linear algebra helpers ----------------
def solve_linear(A, b):
    n = len(A)
    M = [list(map(float, row)) + [float(b[i])] for i, row in enumerate(A)]
    for c in range(n):
        p = max(range(c, n), key=lambda r: abs(M[r][c]))
        if abs(M[p][c]) < 1e-12:
            raise ValueError("singular")
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


class Matrix:
    def __init__(self, rows):
        self.m = [[float(x) for x in r] for r in rows]

    @property
    def shape(self):
        return len(self.m), len(self.m[0])

    def __add__(s, o):
        if not isinstance(o, Matrix) or s.shape != o.shape:
            raise ValueError("Dim")
        return Matrix([[a + b for a, b in zip(r, q)] for r, q in zip(s.m, o.m)])

    def __sub__(s, o):
        return s + (-o)

    def __neg__(s):
        return Matrix([[-a for a in r] for r in s.m])

    def __mul__(s, o):
        if isinstance(o, Matrix):
            if s.shape[1] != o.shape[0]:
                raise ValueError("Dim")
            return Matrix([[sum(a * b for a, b in zip(r, c)) for c in zip(*o.m)] for r in s.m])
        if isinstance(o, (int, float)):
            return Matrix([[a * o for a in r] for r in s.m])
        raise TypeError

    def __rmul__(s, o):
        if isinstance(o, (int, float)):
            return s * o
        raise TypeError

    def __truediv__(s, o):
        return s * (1 / o)

    def det(s):
        r, c = s.shape
        if r != c:
            raise ValueError("Dim")
        M = [row[:] for row in s.m]
        d = 1.0
        for i in range(r):
            p = max(range(i, r), key=lambda k: abs(M[k][i]))
            if abs(M[p][i]) < 1e-14:
                return 0.0
            if p != i:
                M[i], M[p] = M[p], M[i]
                d = -d
            d *= M[i][i]
            for k in range(i + 1, r):
                f = M[k][i] / M[i][i]
                M[k] = [x - f * y for x, y in zip(M[k], M[i])]
        return d

    def trn(s):
        return Matrix([list(c) for c in zip(*s.m)])

    def inv(s):
        n = s.shape[0]
        if s.shape[0] != s.shape[1] or abs(s.det()) < 1e-12:
            raise ValueError("Math")
        cols = [solve_linear(s.m, [1.0 if i == j else 0.0 for i in range(n)]) for j in range(n)]
        return Matrix(cols).trn()

    def __pow__(s, n):
        if n == -1:
            return s.inv()
        if n != int(n) or n < 0 or s.shape[0] != s.shape[1]:
            raise ValueError("Math")
        out = Matrix([[1.0 if i == j else 0.0 for j in range(s.shape[0])] for i in range(s.shape[0])])
        for _ in range(int(n)):
            out = out * s
        return out


class Vec:
    def __init__(self, v):
        self.v = [float(x) for x in v]

    def __add__(s, o):
        if not isinstance(o, Vec) or len(o.v) != len(s.v):
            raise ValueError("Dim")
        return Vec([a + b for a, b in zip(s.v, o.v)])

    def __neg__(s):
        return Vec([-a for a in s.v])

    def __sub__(s, o):
        return s + (-o)

    def __abs__(s):
        return math.sqrt(sum(a * a for a in s.v))

    def __rmul__(s, o):
        if isinstance(o, (int, float)):
            return Vec([a * o for a in s.v])
        raise TypeError

    def __mul__(s, o):
        if isinstance(o, (int, float)):
            return s.__rmul__(o)
        if isinstance(o, Vec) and len(o.v) == len(s.v) and len(s.v) in (2, 3):   # cross product
            a = s.v + [0.0] * (3 - len(s.v))
            b = o.v + [0.0] * (3 - len(o.v))
            return Vec([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])
        raise ValueError("Dim")


# ---------------- polynomial / stats helpers ----------------
def poly_roots(c):
    """c = [a_n ... a_0]  (Durand-Kerner)"""
    n = len(c) - 1
    c = [x / c[0] for x in c]
    r = [(0.4 + 0.9j) ** k for k in range(n)]
    for _ in range(600):
        for i in range(n):
            p = 0
            for coef in c:
                p = p * r[i] + coef
            d = 1
            for j in range(n):
                if j != i:
                    d *= r[i] - r[j]
            r[i] -= p / d
    out = []
    for z in r:
        z = complex(z)
        out.append(z.real if abs(z.imag) < 1e-9 * max(1, abs(z)) else z)
    return out


def quad_roots(a, b, c):
    s = cmath.sqrt(b * b - 4 * a * c)
    out = [(-b + s) / (2 * a), (-b - s) / (2 * a)]
    return [z.real if abs(z.imag) < 1e-12 else z for z in out]


STAT_TYPES = ["1-VAR", "A+BX", "_+CX²", "ln X", "e^X", "A·B^X", "A·X^B", "1/X"]


def quartiles(xs):
    s = sorted(xs)
    n = len(s)

    def med(a):
        m = len(a)
        return a[m // 2] if m % 2 else (a[m // 2 - 1] + a[m // 2]) / 2
    low, up = s[:n // 2], s[(n + 1) // 2:]
    return (med(low) if low else s[0]), med(s), (med(up) if up else s[-1])


def stat_calc(kind, xs, ys):
    n = len(xs)
    if n == 0:
        raise ValueError

    def sd(v):
        m = sum(v) / len(v)
        ss = sum((a - m) ** 2 for a in v)
        return m, math.sqrt(ss / len(v)), (math.sqrt(ss / (len(v) - 1)) if len(v) > 1 else float("nan"))
    mx, px, sx_ = sd(xs)
    items = [("n", n), ("x̄", mx), ("Σx", sum(xs)), ("Σx²", sum(a * a for a in xs)), ("σx", px), ("sx", sx_)]
    if kind == 0:
        q1, md, q3 = quartiles(xs)
        return items + [("minX", min(xs)), ("Q1", q1), ("Med", md), ("Q3", q3), ("maxX", max(xs))], None
    my, py, sy_ = sd(ys)
    items += [("ȳ", my), ("Σy", sum(ys)), ("Σy²", sum(b * b for b in ys)), ("σy", py), ("sy", sy_),
              ("Σxy", sum(a * b for a, b in zip(xs, ys)))]
    if kind == 2:
        def S(p, q=0):
            return sum(u ** p * w ** q for u, w in zip(xs, ys))
        a, b, c = solve_linear([[n, S(1), S(2)], [S(1), S(2), S(3)], [S(2), S(3), S(4)]], [S(0, 1), S(1, 1), S(2, 1)])
        items += [("A", a), ("B", b), ("C", c)]

        def xh(y):
            d = b * b - 4 * c * (a - y)
            if c == 0 or d < 0:
                raise ValueError
            return [(-b + math.sqrt(d)) / (2 * c), (-b - math.sqrt(d)) / (2 * c)]
        return items, (lambda x: [a + b * x + c * x * x], xh)
    idf = lambda t: t
    tf = {1: (idf, idf), 3: (math.log, idf), 4: (idf, math.log), 5: (idf, math.log),
          6: (math.log, math.log), 7: (lambda t: 1 / t, idf)}[kind]
    X, Y = [tf[0](t) for t in xs], [tf[1](t) for t in ys]
    mX, mY = sum(X) / n, sum(Y) / n
    sxx = sum((p - mX) ** 2 for p in X)
    syy = sum((q - mY) ** 2 for q in Y)
    sxy = sum((p - mX) * (q - mY) for p, q in zip(X, Y))
    if sxx == 0:
        raise ValueError
    b = sxy / sxx
    a = mY - b * mX
    r = sxy / math.sqrt(sxx * syy) if syy else float("nan")
    A = math.exp(a) if kind in (4, 5, 6) else a
    B = math.exp(b) if kind == 5 else b
    items += [("A", A), ("B", B), ("r", r)]
    yh = {1: lambda x: A + B * x, 3: lambda x: A + B * math.log(x), 4: lambda x: A * math.exp(B * x),
          5: lambda x: A * B ** x, 6: lambda x: A * x ** B, 7: lambda x: A + B / x}[kind]
    xh = {1: lambda y: (y - A) / B, 3: lambda y: math.exp((y - A) / B), 4: lambda y: math.log(y / A) / B,
          5: lambda y: math.log(y / A) / math.log(B), 6: lambda y: (y / A) ** (1 / B), 7: lambda y: B / (y - A)}[kind]
    return items, (lambda x: [yh(x)], lambda y: [xh(y)])


# ======================= ENGINE =======================
class Engine:
    def __init__(self):
        self.angle = "DEG"
        self.mode = "COMP"          # COMP CMPLX STAT BASE MATRIX VECTOR
        self.base = "DEC"
        self.disp = ("norm", 1)     # ("fix",n) ("sci",n) ("norm",1|2)
        self.polar = False          # complex r∠θ
        self.mixed = False          # fraction a b/c
        self.freq = False           # STAT frequency column
        self.vars = {k: 0.0 for k in VARS}
        self.ans = 0.0
        self.mats = {k: Matrix([[0] * 3 for _ in range(3)]) for k in "ABC"}
        self.vcts = {k: Vec([0, 0, 0]) for k in "ABC"}
        self.matans = self.vctans = None

    def rad(self, x):
        return x * math.pi / 180 if self.angle == "DEG" else x * math.pi / 200 if self.angle == "GRA" else x

    def unrad(self, r):
        return r * 180 / math.pi if self.angle == "DEG" else r * 200 / math.pi if self.angle == "GRA" else r

    def namespace(self):
        r, u = self.rad, self.unrad
        cx = self.mode == "CMPLX"

        def dual(fr, fc):
            if not cx:
                return fr

            def f(x):
                if isinstance(x, complex):
                    return fc(x)
                try:
                    return fr(x)
                except ValueError:
                    return fc(x)
            return f

        def fact(n):
            if n < 0 or n != int(n) or n > 69:
                raise ValueError
            return math.factorial(int(n))

        def cbrt(x):
            return math.copysign(abs(x) ** (1 / 3), x)

        def rtx(n, x):
            if x < 0 and int(n) % 2 == 1:
                return -((-x) ** (1 / n))
            return x ** (1 / n)

        def integ(f, a, b, n=2000):
            h = (b - a) / n
            s = f(a) + f(b)
            for i in range(1, n):
                s += f(a + i * h) * (4 if i % 2 else 2)
            return s * h / 3

        def deriv(f, x):
            h = 1e-5 * max(1, abs(x))
            return (f(x + h) - f(x - h)) / (2 * h)

        def sigma(f, a, b):
            return sum(f(k) for k in range(int(a), int(b) + 1))

        def pol(x, y):
            self.vars["X"] = math.hypot(x, y)
            self.vars["Y"] = u(math.atan2(y, x))
            return self.vars["X"]

        def rec(rr, t):
            self.vars["X"] = rr * math.cos(r(t))
            self.vars["Y"] = rr * math.sin(r(t))
            return self.vars["X"]

        def deg(x):
            return x if self.angle == "DEG" else x * math.pi / 180 if self.angle == "RAD" else x * 200 / 180

        def dot(a, b):
            if len(a.v) != len(b.v):
                raise ValueError
            return sum(p * q for p, q in zip(a.v, b.v))

        def angle(a, b):
            return u(math.acos(max(-1, min(1, dot(a, b) / (abs(a) * abs(b))))))

        def unitv(a):
            return a * (1 / abs(a))

        def ident(n):
            return Matrix([[1 if i == j else 0 for j in range(int(n))] for i in range(int(n))])

        ns = {
            "__builtins__": {},
            "sin": dual(lambda x: math.sin(r(x)), lambda x: cmath.sin(r(x))),
            "cos": dual(lambda x: math.cos(r(x)), lambda x: cmath.cos(r(x))),
            "tan": dual(lambda x: math.tan(r(x)), lambda x: cmath.tan(r(x))),
            "asin": lambda x: u(math.asin(x)), "acos": lambda x: u(math.acos(x)),
            "atan": lambda x: u(math.atan(x)),
            "sinh": math.sinh, "cosh": math.cosh, "tanh": math.tanh,
            "asinh": math.asinh, "acosh": math.acosh, "atanh": math.atanh,
            "log": dual(math.log10, cmath.log10), "ln": dual(math.log, cmath.log),
            "logb": lambda b, x: math.log(x, b),
            "sqrt": dual(math.sqrt, cmath.sqrt), "cbrt": cbrt, "rtx": rtx, "abs": abs,
            "fact": fact, "npr": lambda n, k: math.perm(int(n), int(k)),
            "ncr": lambda n, k: math.comb(int(n), int(k)),
            "rnd": lambda x: float(f"{x:.9g}"), "ran": lambda: round(random.random(), 3),
            "integ": integ, "deriv": deriv, "sigma": sigma, "pol": pol, "rec": rec, "deg": deg,
            "arg": lambda z: u(cmath.phase(z)), "conjg": lambda z: complex(z).conjugate(),
            "cre": lambda z: complex(z).real, "cim": lambda z: complex(z).imag,
            "det": lambda m: m.det(), "trn": lambda m: m.trn(), "ident": ident,
            "dot": dot, "angle": angle, "unitv": unitv,
            "pi": math.pi, "e": math.e, "Ans": self.ans, "MatAns": self.matans, "VctAns": self.vctans,
        }
        if cx:
            ns["i"] = 1j
        for k in "ABC":
            ns["Mat" + k] = self.mats[k]
            ns["Vct" + k] = self.vcts[k]
        ns.update(self.vars)
        return ns

    @staticmethod
    def fix_fact(s):
        while "!" in s:
            i = s.index("!")
            if i > 0 and s[i - 1] == ")":
                depth, k = 0, i - 1
                while k >= 0:
                    if s[k] == ")":
                        depth += 1
                    elif s[k] == "(":
                        depth -= 1
                        if depth == 0:
                            break
                    k -= 1
                start = max(k, 0)
            else:
                k = i
                while k > 0 and (s[k - 1].isalnum() or s[k - 1] == "."):
                    k -= 1
                start = k
            s = s[:start] + "fact(" + s[start:i] + ")" + s[i + 1:]
        return s

    def evaluate_tokens(self, tokens):
        if self.mode == "BASE":
            return self.eval_base(tokens)
        return self.evaluate("".join(tokens))

    def evaluate(self, expr, extra=None):
        s = expr
        if s.count("(") > s.count(")"):
            s += ")" * (s.count("(") - s.count(")"))
        s = re.sub(r"(\d+\.?\d*)°", r"deg(\1)", s).replace("°", "")
        s = self.fix_fact(s)
        for name, py in (("nPr", "npr"), ("nCr", "ncr")):
            s = re.sub(r"(\w+|\([^()]*\))" + name + r"(\w+|\([^()]*\))", py + r"(\1,\2)", s)
        # Casio style:  n ˣ√ ( x   ->  rtx(n, x)
        s = re.sub(r"(\d+\.?\d*|[A-FXYM]|Ans|π|ℯ|\([^()]*\))ˣ√\(", r"rtx(\1,", s)
        for a, b in REPL:
            s = s.replace(a, b)
        s = re.sub(r"(?<![\w.])(\d+\.?\d*)(?=[A-Za-z_(])", r"\1*", s)
        s = re.sub(r"\b(pi|e|Ans|i|[A-FXYM])\b(?=\s*[\w(])", r"\1*", s)
        s = re.sub(r"\)(?=\s*[\w(])", ")*", s)
        ns = self.namespace()
        if extra:
            ns.update(extra)
        v = eval(s, ns)
        if isinstance(v, (Matrix, Vec)):
            return v
        if isinstance(v, complex):
            if self.mode != "CMPLX":
                raise ValueError
            return v
        return float(v)

    def eval_base(self, tokens):
        B = {"DEC": 10, "HEX": 16, "BIN": 2, "OCT": 8}[self.base]
        PB = {"d": 10, "h": 16, "b": 2, "o": 8}
        out, num, nb = [], "", B
        sym = {"×": "*", "÷": "//", "and": "&", "or": "|", "xor": "^", "xnor": "^~", "Not(": "~(", "Neg(": "-("}

        def flush():
            nonlocal num, nb
            if num:
                out.append(str(int(num, nb)))
            num, nb = "", B
        for t in tokens:
            t = t.strip()
            if t in PB:
                flush()
                nb = PB[t]
            elif len(t) == 1 and t in "0123456789ABCDEF":
                num += t
            else:
                flush()
                out.append(sym.get(t, t))
        flush()
        s = "".join(out)
        s += ")" * (s.count("(") - s.count(")"))
        try:
            a = int(self.ans)
        except (TypeError, ValueError, OverflowError):
            a = 0
        return int(eval(s, {"__builtins__": {}, "Ans": a}))

    # ---------- formatting ----------
    def fmt_base(self, v):
        bits = 16 if self.base == "BIN" else 32
        if not -(1 << (bits - 1)) <= v < (1 << (bits - 1)):
            raise ArithmeticError
        if self.base == "DEC":
            return str(v)
        n = v & ((1 << bits) - 1)
        return {"HEX": format(n, "X"), "BIN": format(n, "b"), "OCT": format(n, "o")}[self.base]

    def fmt_c(self, z):
        z = complex(z)
        tol = 1e-12 * max(1, abs(z))
        re_ = 0.0 if abs(z.real) < tol else z.real
        im_ = 0.0 if abs(z.imag) < tol else z.imag
        if im_ == 0:
            return self.fmt(re_)
        if self.polar:
            return f"{self.fmt(abs(z))}∠{self.fmt(self.unrad(cmath.phase(z)))}"
        ims = "" if abs(abs(im_) - 1) < 1e-12 else self.fmt(abs(im_))
        sign = "-" if im_ < 0 else ("+" if re_ != 0 else "")
        return f"{self.fmt(re_) if re_ != 0 else ''}{sign}{ims}i"

    def fmt(self, v):
        if self.mode == "BASE" and isinstance(v, int):
            return self.fmt_base(v)
        if isinstance(v, complex):
            return self.fmt_c(v)
        v = float(v)
        if v != v or v in (math.inf, -math.inf):
            raise ArithmeticError
        if abs(v) < 1e-13:
            v = 0.0
        kind, n = self.disp
        if kind == "fix":
            return f"{v:.{n}f}"
        if kind == "sci":
            m, e = f"{v:.{(n or 10) - 1}e}".split("e")
            return f"{m}×10^{int(e)}"
        lo = 1e-2 if n == 1 else 1e-9
        if v != 0 and (abs(v) < lo or abs(v) >= 1e10):
            m, e = f"{v:.9e}".split("e")
            m = m.rstrip("0").rstrip(".")
            return f"{m}×10^{int(e)}"
        return f"{v:.10g}"

    def exact(self, v):
        """S<=>D"""
        if isinstance(v, (complex, Matrix, Vec)) or v == 0 or abs(v) > 1e9:
            return None
        fr = Fraction(v).limit_denominator(10000)
        if fr.denominator != 1 and abs(float(fr) - v) < 1e-10 * max(1, abs(v)):
            if self.mixed and abs(fr) > 1:
                w = int(abs(fr))
                rem = abs(fr) - w
                return f"{'-' if fr < 0 else ''}{w} {rem.numerator}/{rem.denominator}"
            return f"{fr.numerator}/{fr.denominator}"
        r = Fraction(v / math.pi).limit_denominator(1000)
        if abs(float(r) * math.pi - v) < 1e-9:
            n, d = r.numerator, r.denominator
            top = "π" if abs(n) == 1 else f"{abs(n)}π"
            return ("-" if n < 0 else "") + top + (f"/{d}" if d != 1 else "")
        n = round(v * v)
        if 0 < n < 100000 and abs(v * v - n) < 1e-8:
            k = 1
            for i in range(int(n ** 0.5), 1, -1):
                if n % (i * i) == 0:
                    k = i
                    break
            m = n // (k * k)
            if m > 1:
                return ("-" if v < 0 else "") + (str(k) if k > 1 else "") + f"√{m}"
        return None

    def eng(self, v):
        if isinstance(v, complex) or v == 0:
            return self.fmt(v)
        e = int(math.floor(math.log10(abs(v)) / 3)) * 3
        return self.fmt(v / 10 ** e) + (f"×10^{e}" if e else "")

    def solve(self, expr, guess):
        if "=" in expr:
            l, r = expr.split("=", 1)
            expr = f"({l})-({r})"
        f = lambda x: self.evaluate(expr, {"X": x})
        x = guess
        for _ in range(200):
            fx = f(x)
            if abs(fx) < 1e-12:
                break
            h = 1e-6 * max(1, abs(x))
            d = (f(x + h) - f(x - h)) / (2 * h)
            if d == 0:
                raise ValueError
            nx = x - fx / d
            if abs(nx - x) < 1e-13 * max(1, abs(x)):
                x = nx
                break
            x = nx
        if abs(f(x)) > 1e-7:
            raise ValueError
        return x


# ======================= GUI =======================
# Premium dark UI palette — independent branding, same calculator engine
BG, LCD_BG = "#0b1220", "#111c2e"
SHIFT_C, ALPHA_C, GREEN_C = "#fbbf24", "#fb7185", "#5eead4"
COL = {"fn": ("#1e293b", "#f1f5f9"), "num": ("#263449", "#f8fafc"), "op": ("#164e63", "#ecfeff"),
       "del": ("#0f766e", "#ffffff"), "sh": ("#123047", SHIFT_C), "al": ("#33213a", ALPHA_C),
       "mod": ("#334155", "#e2e8f0")}
MODES = ["1: COMP", "2: CMPLX", "3: STAT", "4: BASE-N", "5: EQN", "6: MATRIX", "7: TABLE", "8: VECTOR"]
MODE_IDS = ["COMP", "CMPLX", "STAT", "BASE", "EQN", "MATRIX", "TABLE", "VECTOR"]
BASE_KEYS = {"²": "DEC", "^": "HEX", "log(": "BIN", "ln(": "OCT"}


def K(label, sh=None, al=None, kind="fn", n=None, s=None, a=None, gr=None):
    return (label, sh, al, kind, n, s, a, gr)


def trig(name):
    al = {"sin": "D", "cos": "E", "tan": "F"}[name]
    return K(name, name + "⁻¹", al, "fn", name + "(", name + "⁻¹(", al)


TOP_L = [K("SHIFT", kind="sh", n="@shift"), K("ALPHA", kind="al", n="@alpha")]
TOP_R = [K("MODE", "SETUP", None, "mod", "@mode", "@setup"), K("ON", "OFF", None, "mod", "@ac", "@off")]
ROW2_L = [K("CALC", "SOLVE", "=", "fn", "@calc", "@solve", "="), K("∫▫", "d/dx", ":", "fn", "∫(", "d/dx(", ":")]
ROW2_R = [K("x⁻¹", "x!", None, "fn", "⁻¹", "!"), K("log▫", "Σ", None, "fn", "logb(", "Σ(")]
ROWS6 = [
    [K("a b/c", None, None, "fn", "÷"), K("√", "∛", None, "fn", "√(", "∛("),
     K("x²", "x³", None, "fn", "²", "³", None, "DEC"), K("x■", "ˣ√", None, "fn", "^", "ˣ√(", None, "HEX"),
     K("log", "10ˣ", None, "fn", "log(", "10^(", None, "BIN"), K("ln", "eˣ", None, "fn", "ln(", "ℯ^(", None, "OCT")],
    [K("(−)", None, "A", "fn", "-", None, "A"), K("° ' \"", None, "B", "fn", "°", None, "B"),
     K("hyp", "Abs", "C", "fn", "@do_hyp", "Abs(", "C"), trig("sin"), trig("cos"), trig("tan")],
    [K("RCL", "STO", None, "fn", "@rcl", "@sto"), K("ENG", "i", None, "fn", "@do_eng", "i"),
     K("(", None, None, "fn", "("), K(")", ",", "X", "fn", ")", ",", "X"),
     K("S⇔D", None, "Y", "fn", "@sd", None, "Y"), K("M+", "M−", "M", "fn", "@mplus", "@mminus", "M")],
]
ROWS5 = [
    [K("7", "CONST", None, "num", "7", "@const"), K("8", "CONV", None, "num", "8", "@conv"),
     K("9", "CLR", None, "num", "9", "@clr"), K("DEL", "INS", None, "del", "@del", "@ins"),
     K("AC", "OFF", None, "del", "@ac", "@off")],
    [K("4", "MATRIX", None, "num", "4", "@matrix_menu"), K("5", "VECTOR", None, "num", "5", "@vector_menu"),
     K("6", None, None, "num", "6"), K("×", "nPr", None, "op", "×", "nPr"), K("÷", "nCr", None, "op", "÷", "nCr")],
    [K("1", "STAT", None, "num", "1", "@stat"), K("2", "CMPLX", None, "num", "2", "@cmplx_menu"),
     K("3", "BASE", None, "num", "3", "@base_menu"), K("+", "Pol", None, "op", "+", "Pol("),
     K("−", "Rec", None, "op", "-", "Rec(")],
    [K("0", "Rnd", None, "num", "0", "Rnd("), K("•", "Ran#", None, "num", ".", "Ran#"),
     K("×10ˣ", "π", "e", "num", "×10^", "π", "ℯ"), K("Ans", "DRG▸", None, "num", "Ans", "@drg"),
     K("=", None, None, "op", "@equals")],
]
AFTER_RESULT_OPS = {"+", "-", "×", "÷", "^", "²", "³", "⁻¹", "!", "nPr", "nCr", "×10^"}


def numtoken(v):
    s = repr(float(v))
    if "e" in s:
        m, e = s.split("e")
        s = f"{m}×10^{int(e)}"
    return "(" + s + ")"


class App:
    def __init__(self, root):
        self.root = root
        self.eng = Engine()
        self.tokens, self.pos = [], 0
        self.result, self.error = "", False
        self.val, self.show_frac = None, False
        self.shift = self.alpha = self.hyp = False
        self.pending = None
        self.just_evaled = False
        self.overwrite = False
        self.history, self.hidx = [], 0
        self.cursor_on = True
        self.stat_type, self.stat_text = 0, ""
        self.build()
        root.bind("<Key>", self.on_keyboard)
        self.blink()
        self.refresh()

    # ---------- UI ----------
    def build(self):
        r = self.root
        r.title("My Calculator | Scientific")
        r.configure(bg=BG)
        r.resizable(False, False)
        r.option_add("*Font", ("Segoe UI", 10))
        head = tk.Frame(r, bg=BG)
        head.pack(fill="x", padx=16, pady=(12, 4))
        mark = tk.Label(head, text="∑", fg="#5eead4", bg="#123047",
                        font=("Segoe UI Symbol", 18, "bold"), width=2, pady=2)
        mark.pack(side="left", padx=(0, 10))
        title_box = tk.Frame(head, bg=BG)
        title_box.pack(side="left", fill="y")
        tk.Label(title_box, text="MY CALCULATOR", fg="#f8fafc", bg=BG,
                 font=("Segoe UI", 15, "bold")).pack(anchor="w")
        tk.Label(title_box, text="SCIENTIFIC  •  PRECISION", fg="#94a3b8", bg=BG,
                 font=("Segoe UI", 8, "bold")).pack(anchor="w", pady=(1, 0))
        tk.Frame(r, bg="#243449", height=1).pack(fill="x", padx=14, pady=(6, 8))
        lcd = tk.Frame(r, bg="#263449", bd=0, highlightthickness=1, highlightbackground="#334155")
        lcd.pack(padx=14, pady=(0, 10), fill="x")
        inner = tk.Frame(lcd, bg=LCD_BG)
        inner.pack(fill="x", padx=1, pady=1)
        self.ind = tk.Label(inner, bg=LCD_BG, fg="#5eead4", font=("Consolas", 9, "bold"), anchor="w")
        self.ind.pack(fill="x", padx=9, pady=(6, 0))
        self.top = tk.Label(inner, bg=LCD_BG, fg="#cbd5e1", font=("Consolas", 15), anchor="w", width=26)
        self.top.pack(fill="x", padx=9, pady=(6, 0))
        self.bot = tk.Label(inner, bg=LCD_BG, fg="#f8fafc", font=("Consolas", 25, "bold"), anchor="e")
        self.bot.pack(fill="x", padx=9, pady=(0, 8))

        body = tk.Frame(r, bg=BG)
        body.pack(padx=10, pady=4)
        top = tk.Frame(body, bg=BG)
        top.pack(fill="x")
        for c in range(6):
            top.columnconfigure(c, weight=1, uniform="k")
        for i, sp in enumerate(TOP_L):
            self.key(top, sp, 0, i, w=4)
        for i, sp in enumerate(TOP_R):
            self.key(top, sp, 0, 4 + i, w=4)
        for i, sp in enumerate(ROW2_L):
            self.key(top, sp, 1, i, w=4)
        for i, sp in enumerate(ROW2_R):
            self.key(top, sp, 1, 4 + i, w=4)
        pad = tk.Frame(top, bg=BG)
        pad.grid(row=0, column=2, columnspan=2, rowspan=2)
        for txt, rr, cc, act in (("▲", 0, 1, "@up"), ("◀", 1, 0, "@left"), ("▶", 1, 2, "@right"), ("▼", 2, 1, "@down")):
            b = tk.Label(pad, text=txt, bg="#334155", fg="#e2e8f0", width=3, relief="flat", bd=0, cursor="hand2", font=("Segoe UI", 9, "bold"))
            b.grid(row=rr, column=cc, padx=2, pady=2, ipady=2)
            b.bind("<Button-1>", lambda e, a=act: self.press(K("", n=a)))
        for row in ROWS6:
            f = tk.Frame(body, bg=BG)
            f.pack(fill="x")
            for c in range(6):
                f.columnconfigure(c, weight=1, uniform="k")
            for i, sp in enumerate(row):
                self.key(f, sp, 0, i, w=4)
        for row in ROWS5:
            f = tk.Frame(body, bg=BG)
            f.pack(fill="x")
            for c in range(5):
                f.columnconfigure(c, weight=1, uniform="k")
            for i, sp in enumerate(row):
                self.key(f, sp, 0, i, w=6, big=True)

    def key(self, parent, spec, row, col, w=5, big=False):
        label, sh, al, kind, gr = spec[0], spec[1], spec[2], spec[3], spec[7]
        cell = tk.Frame(parent, bg=BG)
        cell.grid(row=row, column=col, padx=3, pady=2, sticky="ew")
        hint = tk.Frame(cell, bg=BG)
        hint.pack()
        tk.Label(hint, text=sh or " ", fg=SHIFT_C, bg=BG, font=("Arial", 7, "bold")).pack(side="left")
        if gr:
            tk.Label(hint, text=" " + gr, fg=GREEN_C, bg=BG, font=("Arial", 7, "bold")).pack(side="left")
        if al:
            tk.Label(hint, text=" " + al, fg=ALPHA_C, bg=BG, font=("Arial", 7, "bold")).pack(side="left")
        bg, fg = COL[kind]
        b = tk.Label(cell, text=label, bg=bg, fg=fg, width=w, relief="flat", bd=0, cursor="hand2",
                     font=("Segoe UI", 12 if big else 10, "bold"),
                     highlightthickness=1, highlightbackground="#334155", highlightcolor="#5eead4")
        b.pack(fill="x", ipady=6 if big else 4)
        b.bind("<Button-1>", lambda e, s=spec, w=b: (w.config(relief="sunken", highlightbackground="#5eead4"), self.press(s)))
        b.bind("<ButtonRelease-1>", lambda e, w=b: w.config(relief="flat", highlightbackground="#334155"))

    def blink(self):
        self.cursor_on = not self.cursor_on
        self.refresh()
        self.root.after(550, self.blink)

    def refresh(self):
        left, right = "".join(self.tokens[:self.pos]), "".join(self.tokens[self.pos:])
        cur = ("▮" if self.overwrite else "│") if self.cursor_on and not self.just_evaled else " "
        s = left + cur + right
        width = 26
        if len(s) > width:
            start = max(0, min(len(left) - width + 6, len(s) - width))
            s = s[start:start + width]
        self.top.config(text=s)
        self.bot.config(text=self.result)
        e = self.eng
        parts = []
        if self.shift: parts.append("S")
        if self.alpha: parts.append("A")
        if self.hyp: parts.append("hyp")
        if e.vars["M"]: parts.append("M")
        if self.pending: parts.append(self.pending + "▸ALPHA+var")
        if e.mode != "COMP":
            parts.append({"CMPLX": "CMPLX", "STAT": "STAT", "BASE": "BASE-" + e.base[0],
                          "MATRIX": "MAT", "VECTOR": "VCT"}[e.mode])
        parts.append({"DEG": "D", "RAD": "R", "GRA": "G"}[e.angle])
        if e.disp[0] != "norm": parts.append(e.disp[0].capitalize() + str(e.disp[1]))
        self.ind.config(text="  ".join(parts))

    # ---------- generic dialogs ----------
    def pick(self, title, items, start=0):
        top = tk.Toplevel(self.root)
        top.title(title)
        top.transient(self.root)
        fr = tk.Frame(top)
        fr.pack(padx=6, pady=6)
        lb = tk.Listbox(fr, width=42, height=min(len(items), 16), font=("Consolas", 11), exportselection=False)
        sb = tk.Scrollbar(fr, command=lb.yview)
        lb.config(yscrollcommand=sb.set)
        lb.pack(side="left")
        sb.pack(side="left", fill="y")
        for it in items:
            lb.insert("end", it)
        lb.selection_set(start)
        lb.see(start)
        res = [None]

        def ok(e=None):
            s = lb.curselection()
            if s:
                res[0] = s[0]
                top.destroy()
        lb.bind("<Double-Button-1>", ok)
        lb.bind("<Return>", ok)
        top.bind("<Escape>", lambda e: top.destroy())
        tk.Button(top, text="OK", command=ok).pack(pady=(0, 6))
        lb.focus_set()
        top.grab_set()
        self.root.wait_window(top)
        return res[0]

    def edit_grid(self, title, data, vector=False):
        """Editable grid. Returns list-of-lists (floats) or None."""
        top = tk.Toplevel(self.root)
        top.title(title)
        top.transient(self.root)
        rv, cv = tk.IntVar(value=len(data)), tk.IntVar(value=1 if vector else len(data[0]))
        entries, res = [], [None]
        ctl = tk.Frame(top)
        ctl.pack(padx=8, pady=6)
        box = tk.Frame(top)
        box.pack(padx=8, pady=4)

        def build():
            old = [[e.get() for e in row] for row in entries]
            for w in box.winfo_children():
                w.destroy()
            entries.clear()
            r, c = rv.get(), (1 if vector else cv.get())
            for i in range(r):
                row = []
                for j in range(c):
                    e = tk.Entry(box, width=9, justify="center")
                    e.grid(row=i, column=j, padx=2, pady=2)
                    if i < len(old) and j < len(old[i]):
                        e.insert(0, old[i][j])
                    elif not old:
                        e.insert(0, format(data[i][j], "g"))
                    else:
                        e.insert(0, "0")
                    row.append(e)
                entries.append(row)
        tk.Label(ctl, text="Dim:" if vector else "Rows:").pack(side="left")
        tk.Spinbox(ctl, from_=1, to=3, width=3, textvariable=rv, command=build).pack(side="left")
        if not vector:
            tk.Label(ctl, text="  Cols:").pack(side="left")
            tk.Spinbox(ctl, from_=1, to=3, width=3, textvariable=cv, command=build).pack(side="left")
        build()

        def save():
            try:
                res[0] = [[self.eng.evaluate(e.get()) for e in row] for row in entries]
            except Exception:
                messagebox.showerror(title, "Kisi cell mein galat value hai.")
                return
            top.destroy()
        tk.Button(top, text="Save", command=save).pack(pady=6)
        top.grab_set()
        self.root.wait_window(top)
        return res[0]

    def show_matrix(self, title, rows):
        top = tk.Toplevel(self.root)
        top.title(title)
        box = tk.Frame(top, padx=12, pady=10)
        box.pack()
        for i, row in enumerate(rows):
            for j, x in enumerate(row):
                try:
                    t = self.eng.fmt(x)
                except Exception:
                    t = "ERR"
                tk.Label(box, text=t, width=12, relief="groove", font=("Consolas", 11)).grid(row=i, column=j, padx=1, pady=1)

    def info_window(self, title, lines):
        top = tk.Toplevel(self.root)
        top.title(title)
        lb = tk.Listbox(top, width=44, height=min(max(len(lines), 3), 20), font=("Consolas", 11))
        lb.pack(padx=8, pady=8)
        for l in lines:
            lb.insert("end", l)

    # ---------- editing ----------
    def clear_error(self):
        if self.error:
            self.error, self.result = False, ""

    def insert(self, t):
        if self.just_evaled:
            self.tokens = ["Ans"] if t in AFTER_RESULT_OPS else []
            self.pos = len(self.tokens)
            self.just_evaled, self.result = False, ""
        if self.overwrite and self.pos < len(self.tokens):
            self.tokens[self.pos] = t
        else:
            self.tokens.insert(self.pos, t)
        self.pos += 1

    def press(self, spec):
        n, s, a = spec[4], spec[5], spec[6]
        if n == "@shift":
            self.shift, self.alpha = not self.shift, False
            return self.refresh()
        if n == "@alpha":
            self.alpha, self.shift = not self.alpha, False
            return self.refresh()
        act = s if self.shift else a if self.alpha else n
        self.shift = self.alpha = False
        self.clear_error()
        if act:
            self.run(act)
        self.refresh()

    def run(self, act):
        e = self.eng
        if self.pending:
            mode, self.pending = self.pending, None
            if act in VARS:
                if mode == "STO":
                    e.vars[act] = e.ans
                else:
                    self.insert(act)
                return
        if act.startswith("@"):
            return getattr(self, act[1:])()
        if e.mode == "BASE":
            if act in BASE_KEYS:
                return self.set_base(BASE_KEYS[act])
            valid = {"BIN": "01", "OCT": "01234567", "DEC": "0123456789", "HEX": "0123456789ABCDEF"}[e.base]
            if (len(act) == 1 and act in "0123456789ABCDEF" and act not in valid) or act == ".":
                return
        if self.hyp and act[:3] in ("sin", "cos", "tan") and act.endswith("("):
            act = act[:3] + "h" + ("⁻¹" if "⁻¹" in act else "") + "("
        self.hyp = False
        self.insert(act)

    def on_keyboard(self, ev):
        ks, ch = ev.keysym, ev.char
        self.clear_error()
        if ks in ("Return", "KP_Enter") or ch == "=":
            self.equals()
        elif ks == "BackSpace":
            self.delete_()
        elif ks == "Escape":
            self.ac()
        elif ks in ("Left", "Right", "Up", "Down"):
            getattr(self, ks.lower())()
        elif ch and ch in "0123456789.+-()!^":
            self.run(ch) if self.eng.mode == "BASE" else self.insert(ch)
        elif ch == "*":
            self.insert("×")
        elif ch == "/":
            self.insert("÷")
        elif ch and self.eng.mode == "BASE" and ch.upper() in "ABCDEF":
            self.run(ch.upper())
        elif ch == "i" and self.eng.mode == "CMPLX":
            self.insert("i")
        self.refresh()

    # ---------- simple actions ----------
    def delete_(self):
        if self.just_evaled:
            self.just_evaled, self.result, self.pos = False, "", len(self.tokens)
        if self.pos > 0:
            del self.tokens[self.pos - 1]
            self.pos -= 1

    def ac(self):
        self.tokens, self.pos, self.result = [], 0, ""
        self.just_evaled = self.error = False
        self.hidx = len(self.history)

    def off(self):
        self.root.destroy()

    def ins(self):
        self.overwrite = not self.overwrite

    def clr(self):
        if messagebox.askyesno("CLR", "Saari memory (A-F, X, Y, M, Ans, Mat, Vct) aur history clear karein?"):
            e = self.eng
            e.vars = {k: 0.0 for k in VARS}
            e.ans, self.val, self.history = 0.0, None, []
            e.mats = {k: Matrix([[0] * 3 for _ in range(3)]) for k in "ABC"}
            e.vcts = {k: Vec([0, 0, 0]) for k in "ABC"}
            self.ac()

    def leave_result(self):
        if self.just_evaled:
            self.just_evaled, self.result, self.pos = False, "", len(self.tokens)

    def left(self):
        self.leave_result()
        self.pos = max(0, self.pos - 1)

    def right(self):
        self.leave_result()
        self.pos = min(len(self.tokens), self.pos + 1)

    def up(self):
        if self.hidx > 0:
            self.hidx -= 1
            self.load()

    def down(self):
        if self.hidx < len(self.history) - 1:
            self.hidx += 1
            self.load()

    def load(self):
        self.tokens = list(self.history[self.hidx])
        self.pos, self.result, self.just_evaled = len(self.tokens), "", False

    def do_hyp(self):
        self.hyp = True

    def drg(self):
        order = ["DEG", "RAD", "GRA"]
        self.eng.angle = order[(order.index(self.eng.angle) + 1) % 3]

    def sto(self):
        self.pending = "STO"

    def rcl(self):
        self.pending = "RCL"

    def mplus(self):
        self.eng.vars["M"] += self.eng.ans

    def mminus(self):
        self.eng.vars["M"] -= self.eng.ans

    def fail(self, ex):
        self.error = True
        self.result = "Syntax ERROR" if isinstance(ex, (SyntaxError, NameError, TypeError)) else "Math ERROR"

    def add_history(self):
        if not self.history or self.history[-1] != self.tokens:
            self.history.append(list(self.tokens))
        self.hidx = len(self.history)

    def equals(self):
        if not self.tokens:
            return
        e = self.eng
        try:
            v = e.evaluate_tokens(self.tokens)
            if isinstance(v, (Matrix, Vec)):
                return self.show_obj(v)
            res = e.fmt(v)
        except Exception as ex:
            return self.fail(ex)
        self.val = e.ans = v
        self.result, self.show_frac, self.just_evaled = res, False, True
        self.add_history()

    def show_obj(self, v):
        if isinstance(v, Matrix):
            self.eng.matans, name, rows = v, "MatAns", v.m
        else:
            self.eng.vctans, name, rows = v, "VctAns", [[x] for x in v.v]
        self.result, self.val, self.just_evaled = name, None, True
        self.add_history()
        self.show_matrix(name, rows)

    def sd(self):
        if self.val is None or self.error or self.eng.mode == "BASE":
            return
        if self.show_frac:
            self.result, self.show_frac = self.eng.fmt(self.val), False
        else:
            ex = self.eng.exact(self.val)
            if ex:
                self.result, self.show_frac = ex, True

    def do_eng(self):
        if self.val is not None and not self.error and self.eng.mode != "BASE":
            self.result = self.eng.eng(self.val)

    def calc(self):
        expr = "".join(self.tokens)
        for n in sorted(set(re.findall(r"\b[A-FXYM]\b", expr))):
            v = simpledialog.askfloat("CALC", f"{n} = ?", initialvalue=self.eng.vars[n], parent=self.root)
            if v is None:
                return
            self.eng.vars[n] = v
        self.equals()

    def solve(self):
        expr = "".join(self.tokens)
        if not expr:
            return
        g = simpledialog.askfloat("SOLVE", "X ka initial value (guess):", initialvalue=self.eng.vars["X"], parent=self.root)
        if g is None:
            return
        try:
            x = self.eng.solve(expr, g)
        except Exception:
            self.error, self.result = True, "Can't Solve"
            return
        self.eng.vars["X"] = self.val = self.eng.ans = x
        self.result, self.just_evaled = "X=" + self.eng.fmt(x), True

    # ---------- MODE / SETUP ----------
    def mode(self):
        i = self.pick("MODE", MODES, MODE_IDS.index(self.eng.mode) if self.eng.mode in MODE_IDS else 0)
        if i is None:
            return
        name = MODE_IDS[i]
        if name == "EQN":
            return self.eqn()
        if name == "TABLE":
            return self.table()
        self.eng.mode = name
        self.ac()
        if name == "STAT":
            self.stat(ask=True)

    def setup(self):
        e = self.eng
        i = self.pick("SETUP", ["1: Angle unit", "2: Number format (Fix/Sci/Norm)", "3: Complex format",
                                "4: Fraction result", "5: STAT frequency"])
        if i is None:
            return
        if i == 0:
            j = self.pick("Angle unit", ["1: Deg", "2: Rad", "3: Gra"])
            if j is not None:
                e.angle = ("DEG", "RAD", "GRA")[j]
        elif i == 1:
            j = self.pick("Number format", ["1: Fix", "2: Sci", "3: Norm 1", "4: Norm 2"])
            if j == 0:
                n = simpledialog.askinteger("Fix", "Decimal places (0-9):", parent=self.root, minvalue=0, maxvalue=9)
                e.disp = ("fix", n) if n is not None else e.disp
            elif j == 1:
                n = simpledialog.askinteger("Sci", "Significant digits (1-10):", parent=self.root, minvalue=1, maxvalue=10)
                e.disp = ("sci", n) if n is not None else e.disp
            elif j in (2, 3):
                e.disp = ("norm", j - 1)
        elif i == 2:
            j = self.pick("Complex format", ["1: a+bi", "2: r∠θ"])
            if j is not None:
                e.polar = bool(j)
        elif i == 3:
            j = self.pick("Fraction result", ["1: d/c (improper)", "2: ab/c (mixed)"])
            if j is not None:
                e.mixed = bool(j)
        elif i == 4:
            j = self.pick("STAT frequency", ["1: ON", "2: OFF"])
            if j is not None:
                e.freq = (j == 0)
        self.rerender()

    def rerender(self):
        if self.just_evaled and self.val is not None and not self.error:
            try:
                self.result = self.eng.fmt(self.val)
            except Exception:
                pass

    # ---------- CONST / CONV ----------
    def const(self):
        i = self.pick("CONST", [f"{k + 1:02d}: {s:<4} {n}" for k, (s, n, v) in enumerate(CONSTS)])
        if i is not None:
            self.insert(numtoken(CONSTS[i][2]))

    def conv(self):
        i = self.pick("CONV", [f"{k + 1:02d}: {n}" for k, (n, f) in enumerate(CONV)])
        if i is None:
            return
        e = self.eng
        try:
            x = e.evaluate_tokens(self.tokens) if self.tokens else e.ans
            y = CONV[i][1](float(x))
            res = e.fmt(y)
        except Exception as ex:
            return self.fail(ex)
        self.val = e.ans = y
        self.result, self.show_frac, self.just_evaled = res, False, True

    # ---------- CMPLX / BASE menus ----------
    def cmplx_menu(self):
        i = self.pick("CMPLX", ["1: Arg(", "2: Conjg(", "3: Re(", "4: Im(", "5: ▸r∠θ", "6: ▸a+bi"])
        if i is None:
            return
        if i < 4:
            self.insert(("Arg(", "Conjg(", "Re(", "Im(")[i])
        else:
            self.eng.polar = (i == 4)
            self.rerender()

    def base_menu(self):
        items = ["d", "h", "b", "o", "Neg(", "Not(", " and ", " or ", " xor ", " xnor "]
        i = self.pick("BASE-N", ["1: d  (decimal prefix)", "2: h  (hex prefix)", "3: b  (binary prefix)",
                                 "4: o  (octal prefix)", "5: Neg(", "6: Not(", "7: and", "8: or", "9: xor", "10: xnor"])
        if i is not None:
            self.insert(items[i])

    def set_base(self, b):
        self.eng.base = b
        if self.just_evaled and isinstance(self.val, int):
            try:
                self.result = self.eng.fmt(self.val)
            except Exception:
                self.fail(ArithmeticError())

    # ---------- MATRIX / VECTOR ----------
    def matrix_menu(self):
        items = ["Edit MatA", "Edit MatB", "Edit MatC", "MatA", "MatB", "MatC", "MatAns", "det(", "Trn(", "Ident("]
        i = self.pick("MATRIX", [f"{k + 1}: {t}" for k, t in enumerate(items)])
        if i is None:
            return
        if i < 3:
            k = "ABC"[i]
            d = self.edit_grid("Mat" + k, self.eng.mats[k].m)
            if d:
                self.eng.mats[k] = Matrix(d)
        else:
            self.insert(items[i])

    def vector_menu(self):
        items = ["Edit VctA", "Edit VctB", "Edit VctC", "VctA", "VctB", "VctC", "VctAns", "Dot(", "Angle(", "UnitV("]
        i = self.pick("VECTOR", [f"{k + 1}: {t}" for k, t in enumerate(items)])
        if i is None:
            return
        if i < 3:
            k = "ABC"[i]
            d = self.edit_grid("Vct" + k, [[x] for x in self.eng.vcts[k].v], vector=True)
            if d:
                self.eng.vcts[k] = Vec([r[0] for r in d])
        else:
            self.insert(items[i])

    # ---------- STAT ----------
    def stat(self, ask=False):
        if ask:
            i = self.pick("STAT type", [f"{k + 1}: {n}" for k, n in enumerate(STAT_TYPES)], self.stat_type)
            if i is None:
                return
            self.stat_type, self.stat_text = i, ""
        kind, fr, e = self.stat_type, self.eng.freq, self.eng
        two = kind != 0
        cols = "X" + (" Y" if two else "") + (" FREQ" if fr else "")
        top = tk.Toplevel(self.root)
        top.title("STAT - " + STAT_TYPES[kind])
        tk.Label(top, text=f"Har line mein: {cols}   (space ya comma se alag)").pack(padx=8, pady=(8, 2))
        mid = tk.Frame(top)
        mid.pack(padx=8)
        txt = tk.Text(mid, width=24, height=14, font=("Consolas", 11))
        txt.pack(side="left")
        txt.insert("1.0", self.stat_text)
        res = tk.Listbox(mid, width=26, height=14, font=("Consolas", 11), exportselection=False)
        res.pack(side="left", padx=(8, 0))
        state = {"items": [], "pred": None}
        out = tk.Label(top, text="", font=("Consolas", 11))

        def compute():
            self.stat_text = txt.get("1.0", "end").strip()
            xs, ys = [], []
            try:
                for line in self.stat_text.splitlines():
                    p = [e.evaluate(t) for t in re.split(r"[,\s]+", line.strip()) if t]
                    if not p:
                        continue
                    if len(p) != 1 + two + fr:
                        raise ValueError
                    cnt = int(p[-1]) if fr else 1
                    xs += [p[0]] * cnt
                    if two:
                        ys += [p[1]] * cnt
                state["items"], state["pred"] = stat_calc(kind, xs, ys)
            except Exception:
                messagebox.showerror("STAT", "Data ya calculation mein error (format check karo).")
                return
            res.delete(0, "end")
            for n, v in state["items"]:
                try:
                    s = e.fmt(v)
                except Exception:
                    s = "Math ERROR"
                res.insert("end", f"{n:<5}= {s}")

        def insert_sel():
            s = res.curselection()
            if s:
                self.insert(numtoken(state["items"][s[0]][1]))
                self.refresh()

        def predict(which):
            if not state["pred"]:
                return
            try:
                v = e.evaluate(ent.get())
                r = state["pred"][which](v)
                out.config(text=("ŷ = " if which == 0 else "x̂ = ") + ",  ".join(e.fmt(x) for x in r))
            except Exception:
                out.config(text="Math ERROR")
        btn = tk.Frame(top)
        btn.pack(pady=6)
        tk.Button(btn, text="Calculate", command=compute).pack(side="left", padx=4)
        tk.Button(btn, text="Insert selected ▸ expression", command=insert_sel).pack(side="left", padx=4)
        if two:
            pr = tk.Frame(top)
            pr.pack()
            ent = tk.Entry(pr, width=10)
            ent.pack(side="left")
            tk.Button(pr, text="ŷ (x se)", command=lambda: predict(0)).pack(side="left", padx=3)
            tk.Button(pr, text="x̂ (y se)", command=lambda: predict(1)).pack(side="left", padx=3)
        out.pack(pady=4)

    # ---------- EQN ----------
    def eqn(self):
        i = self.pick("EQN", ["1: aX+bY=c", "2: aX+bY+cZ=d", "3: aX²+bX+c=0", "4: aX³+bX²+cX+d=0"])
        if i is None:
            return
        rows, cols, names = [(2, 3, "abc"), (3, 4, "abcd"), (1, 3, "abc"), (1, 4, "abcd")][i]
        top = tk.Toplevel(self.root)
        top.title("EQN")
        grid = tk.Frame(top)
        grid.pack(padx=10, pady=8)
        ents = []
        for r in range(rows):
            row = []
            for c in range(cols):
                tk.Label(grid, text=names[c] + (str(r + 1) if rows > 1 else "") + " =").grid(row=r, column=2 * c)
                en = tk.Entry(grid, width=7, justify="center")
                en.insert(0, "0")
                en.grid(row=r, column=2 * c + 1, padx=2, pady=2)
                row.append(en)
            ents.append(row)
        box = tk.Listbox(top, width=34, height=5, font=("Consolas", 11))
        e = self.eng

        def go():
            try:
                v = [[e.evaluate(x.get()) for x in row] for row in ents]
                if i < 2:
                    n = i + 2
                    sol = solve_linear([r[:n] for r in v], [r[n] for r in v])
                    lines = [f"{'XYZ'[k]} = {e.fmt(sol[k])}" for k in range(n)]
                else:
                    c = v[0]
                    if c[0] == 0:
                        raise ValueError
                    roots = quad_roots(*c) if i == 2 else poly_roots(c)
                    if i == 3:
                        roots = sorted(roots, key=lambda z: (z.real, z.imag) if isinstance(z, complex) else (z, 0))
                    lines = [f"X{k + 1} = {e.fmt(z)}" for k, z in enumerate(roots)]
            except Exception:
                lines = ["Math ERROR (no unique solution)"]
            box.delete(0, "end")
            for l in lines:
                box.insert("end", l)
        tk.Button(top, text="Solve  (=)", command=go).pack(pady=4)
        box.pack(padx=10, pady=(0, 10))

    # ---------- TABLE ----------
    def table(self):
        top = tk.Toplevel(self.root)
        top.title("TABLE  f(X)")
        e = self.eng
        form = tk.Frame(top)
        form.pack(padx=10, pady=8)
        fields = {}
        for r, (lab, dv) in enumerate([("f(X) =", "X^2"), ("Start =", "1"), ("End =", "5"), ("Step =", "1")]):
            tk.Label(form, text=lab).grid(row=r, column=0, sticky="e")
            en = tk.Entry(form, width=24)
            en.insert(0, dv)
            en.grid(row=r, column=1, pady=2)
            fields[lab] = en
        box = tk.Listbox(top, width=40, height=12, font=("Consolas", 11))

        def go():
            box.delete(0, "end")
            try:
                a, b, s = (e.evaluate(fields[k].get()) for k in ("Start =", "End =", "Step ="))
                if s <= 0 or (b - a) / s > 99:
                    raise ValueError
            except Exception:
                return box.insert("end", "Start/End/Step galat (max 100 rows)")
            f, x = fields["f(X) ="].get(), a
            while x <= b + 1e-12:
                try:
                    y = e.fmt(e.evaluate(f, {"X": x}))
                except Exception:
                    y = "ERROR"
                box.insert("end", f"X={e.fmt(x):<10} f(X)={y}")
                x += s
        tk.Button(top, text="Generate table  (=)", command=go).pack()
        box.pack(padx=10, pady=8)


if __name__ == "__main__":
    root = tk.Tk()
    App(root)
    root.mainloop()
