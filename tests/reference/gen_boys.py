#!/usr/bin/env python3
"""High-precision reference values of the Boys function for milestone M1.

    F_n(x) = integral_0^1 t^(2n) exp(-x t^2) dt

Evaluated with 60-digit decimal arithmetic (Python standard library) from the
series with positive terms only,

    F_n(x) = exp(-x) * sum_{k>=0} (2x)^k / [(2n+1)(2n+3)...(2n+2k+1)],

which follows from repeated integration by parts,
F_n(x) = [2x F_{n+1}(x) + exp(-x)] / (2n+1). No cancellation occurs for any
x >= 0. The values are then rounded to the nearest double.

Self-checks of the reference before it is written (the script aborts if one
fails):
  * F_0(x) against sqrt(pi/(4x)) erf(sqrt(x)) (math.erf, double precision),
  * the recursion F_{n+1} = [(2n+1) F_n - exp(-x)] / (2x) at 60 digits,
  * the closed form at x = 0, F_n(0) = 1/(2n+1),
  * the large-x limit F_n(x) -> (2n-1)!! / 2^(n+1) * sqrt(pi / x^(2n+1)).

Writes data/m1/boys.json. Does not use PySCF.

Usage: python tests/reference/gen_boys.py
"""

import json
import math
import pathlib
import platform
from decimal import Decimal, getcontext

getcontext().prec = 60
EPS = Decimal(10) ** -58

N_MAX = 32  # acceptance criterion needs n <= 24; a margin is included


def boys(n: int, x: float) -> Decimal:
    """F_n(x) with about 55 correct digits. `x` is taken as the exact double."""
    xd = Decimal(x)  # exact binary value of the double
    term = Decimal(1) / (2 * n + 1)
    total = term
    two_x = 2 * xd
    k = 0
    while True:
        k += 1
        term = term * two_x / (2 * n + 2 * k + 1)
        total += term
        if term < EPS * total and 2 * n + 2 * k + 1 > two_x:
            break
    return (-xd).exp() * total


def x_values() -> list[float]:
    xs = [0.0]
    # very small arguments
    xs += [10.0 ** e for e in range(-16, 0)]
    # dense in the range where implementations switch between methods
    xs += [0.1 * k for k in range(1, 400)]
    # irregular points up to 500 (not on any tabulation grid)
    xs += [40.0 + 460.0 * (k / 150.0) ** 1.5 + 0.0137 * math.sin(7.0 * k) for k in range(151)]
    # exact ends and a few points beyond the required range
    xs += [500.0, 750.0, 1000.0, 5000.0, 1.0e5]
    return sorted(set(xs))


def self_check(xs: list[float], table: dict) -> None:
    pi = Decimal("3.14159265358979323846264338327950288419716939937510582097494")
    for x in xs:
        row = table[x]
        xd = Decimal(x)
        if x == 0.0:
            for n in range(N_MAX + 1):
                assert row[n] == Decimal(1) / (2 * n + 1)
            continue
        # F_0 against the error function (double precision check)
        f0 = math.sqrt(math.pi / (4.0 * x)) * math.erf(math.sqrt(x))
        assert abs(float(row[0]) - f0) <= 4e-15 * f0, (x, float(row[0]), f0)
        # upward recursion at 60 digits; it loses digits for small x, so the
        # tolerance is scaled with the amplification factor (2n+1)/(2x)
        ex = (-xd).exp()
        for n in range(N_MAX):
            up = ((2 * n + 1) * row[n] - ex) / (2 * xd)
            amplification = max(Decimal(1), Decimal(2 * n + 1) / (2 * xd))
            assert abs(up - row[n + 1]) <= Decimal(10) ** -50 * amplification * row[n], (x, n)
        # Large-x limit. The difference is the tail integral_1^inf t^(2n)
        # exp(-x t^2) dt, which is below exp(-x) for x > n. That is far below
        # the 60-digit working precision, so the working precision (relative
        # 1e-55) is the effective bound here.
        if x >= 200.0:
            for n in range(N_MAX + 1):
                double_factorial = Decimal(1)
                for m in range(1, 2 * n, 2):
                    double_factorial *= m
                limit = double_factorial / Decimal(2) ** (n + 1) * (pi / xd ** (2 * n + 1)).sqrt()
                assert abs(row[n] - limit) <= ex + Decimal(10) ** -55 * limit, (x, n)


def main() -> None:
    xs = x_values()
    table = {x: [boys(n, x) for n in range(N_MAX + 1)] for x in xs}
    self_check(xs, table)

    out = {
        "meta": {
            "generated_by": "tests/reference/gen_boys.py",
            "method": "positive series in 60-digit decimal arithmetic, rounded to double",
            "python_version": platform.python_version(),
            "definition": "F_n(x) = integral_0^1 t^(2n) exp(-x t^2) dt",
            "n_max": N_MAX,
            "layout": "values[i][n] = F_n(x[i])",
        },
        "x": xs,
        # float(str(Decimal)) is the correctly rounded double
        "values": [[float(str(v)) for v in table[x]] for x in xs],
    }
    path = pathlib.Path(__file__).resolve().parent / "data" / "m1" / "boys.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["{", ' "meta": ' + json.dumps(out["meta"]) + ",", ' "x": ' + json.dumps(xs) + ",",
             ' "values": [']
    lines += ["  " + json.dumps(row) + ("," if i + 1 < len(xs) else "")
              for i, row in enumerate(out["values"])]
    lines += [" ]", "}"]
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path.name}: {len(xs)} arguments, n = 0..{N_MAX}  ({path.stat().st_size / 1024:.0f} kB)")


if __name__ == "__main__":
    main()
