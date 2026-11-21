"""Deterministic checks of the proposed stencil, not production simulations.

Run: python3 scripts/validate_scalar_regulator.py
Uses only the Python standard library. Checks finite-difference derivatives,
shift symmetry, the free Fourier eigenvalue, and finite-lattice inequalities.
"""
import itertools
import json
import math
import random


def run():
    n = 4
    points = list(itertools.product(range(n), repeat=4))
    index = {p: i for i, p in enumerate(points)}
    neighbors = []
    for p in points:
        row = []
        for mu in range(4):
            for step in (-1, 1):
                q = list(p)
                q[mu] = (q[mu] + step) % n
                row.append(index[tuple(q)])
        neighbors.append(row)
    v = len(points)

    def lap(s):
        return [sum(s[j] - s[i] for j in row) for i, row in enumerate(neighbors)]

    def quantities(s, coupling):
        d = lap(s)
        q = [0.5 * sum((s[j] - s[i]) ** 2 for j in row)
             for i, row in enumerate(neighbors)]
        b = [x + coupling * y for x, y in zip(d, q)]
        return 0.5 * sum(x * x for x in b), d, q, b

    def force(s, coupling):
        _, _, _, b = quantities(s, coupling)
        db = lap(b)
        return [db[i] + coupling * sum((b[i] + b[j]) * (s[i] - s[j]) for j in row)
                for i, row in enumerate(neighbors)]

    rng = random.Random(20261120)
    field = [rng.gauss(0, 0.3) for _ in points]
    mean = sum(field) / v
    field = [x - mean for x in field]
    kmin = 4 * math.sin(math.pi / n) ** 2
    derivative_errors, bound_ratios, shift_errors = [], [], []
    for coupling in (0.0, 0.1, -0.1, 0.8):
        action, delta, q, _ = quantities(field, coupling)
        norm2 = sum(x * x for x in field)
        bound = coupling**2 * kmin**2 * norm2**2 / (2 * v) if coupling else kmin**2 * norm2 / 2
        assert action >= bound * (1 - 1e-12)
        assert abs(sum(delta)) < 1e-10
        assert abs(sum(q) + sum(x*y for x,y in zip(field,delta))) < 1e-10
        assert sum(q) >= kmin * norm2 * (1 - 1e-12)
        bound_ratios.append(action / bound)
        grad = force(field, coupling)
        assert abs(sum(grad)) < 1e-9
        shifted = quantities([x + 2.3 for x in field], coupling)[0]
        shift_errors.append(abs(shifted-action)/max(1,abs(action)))
        for _ in range(8):
            direction = [rng.uniform(-1,1) for _ in points]
            mean = sum(direction)/v
            direction = [x-mean for x in direction]
            length = math.sqrt(sum(x*x for x in direction))
            direction = [x/length for x in direction]
            h = 1e-5
            plus = quantities([x+h*y for x,y in zip(field,direction)],coupling)[0]
            minus = quantities([x-h*y for x,y in zip(field,direction)],coupling)[0]
            numerical = (plus-minus)/(2*h)
            analytic = sum(x*y for x,y in zip(grad,direction))
            derivative_errors.append(abs(numerical-analytic)/max(1,abs(analytic)))
    wave = [math.cos(2*math.pi*(p[0]+p[1])/n) for p in points]
    eigenvalue = (8*math.sin(math.pi/n)**2)**2
    applied = lap(lap(wave))
    eigen_error = max(abs(x-eigenvalue*y) for x,y in zip(applied,wave))
    assert max(derivative_errors)<1e-6
    assert max(shift_errors)<1e-12
    assert eigen_error<1e-10
    print(json.dumps({
        'status':'passed','scope':'deterministic stencil checks only; no interacting ensemble',
        'N':n,'dimensions':4,'directions_checked':len(derivative_errors),
        'maximum_relative_force_error':max(derivative_errors),
        'maximum_relative_shift_error':max(shift_errors),
        'free_bilaplacian_eigenvalue_error':eigen_error,
        'minimum_action_to_bound_ratio':min(bound_ratios)
    }, indent=2))


if __name__ == '__main__':
    run()
