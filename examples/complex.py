"""Benchmark polytope: a randomly cut box in 10 dimensions.

The polytope is the box [-2, 2]^10 intersected with 28 pseudo-random
halfspaces a_i . x <= b_i, where the coefficients a_ij are integers in
[-999, 999] drawn from a small linear congruential generator and
b_i = sum_j |a_ij| (the cutting plane sits halfway between the origin and
the deepest box corner, so every cut slices the box, creating many
vertices). On the reference machine (2026 laptop, Windows, lrslib v5.1
with lrsmp arithmetic) enumeration takes roughly 8 seconds and yields
40317 vertices.

Everything is deterministic and integer-valued, so the exact same
polytope can be rebuilt in other environments for benchmarking — e.g. in
the R package vertexenum, or against a newer lrslib. To port the
generator, only the LCG below is needed:

    x_{k+1} = (1103515245 * x_k + 12345) mod 2^31,  x_0 = 20260719
    value_k = (x_{k+1} mod 1999) - 999

(the classic glibc rand() constants; all arithmetic fits in 64-bit
integers). Row i of A is values dim*i .. dim*i+dim-1, read row-major.

Run with --ine FILE to also write the system in lrs .ine format, for
feeding directly to a standalone lrs binary.
"""

import argparse
import time

from vertexenum import enumerate_vertices, lrs_version

DIM = 10        # number of variables
CUTS = 28       # number of pseudo-random cutting planes
BOX = 2         # box is [-BOX, BOX]^DIM
SEED = 20260719

LCG_MULT = 1103515245
LCG_INC = 12345
LCG_MOD = 2**31


def lcg_ints(count, seed=SEED):
    """First `count` values of the deterministic coefficient stream."""
    x = seed
    out = []
    for _ in range(count):
        x = (LCG_MULT * x + LCG_INC) % LCG_MOD
        out.append(x % 1999 - 999)
    return out


def build_polytope(dim=DIM, cuts=CUTS, box=BOX):
    """Return (A, b) for the benchmark polytope, all entries integers."""
    vals = lcg_ints(cuts * dim)
    A = []
    b = []
    for i in range(cuts):
        row = vals[i * dim:(i + 1) * dim]
        if all(v == 0 for v in row):  # never happens for the default seed
            row[0] = 1
        A.append(row)
        b.append(sum(abs(v) for v in row))
    for j in range(dim):
        for sign in (1, -1):
            row = [0] * dim
            row[j] = sign
            A.append(row)
            b.append(box)
    return A, b


def write_ine(path, A, b):
    """Write the system in lrs .ine (H-representation) format."""
    with open(path, "w") as f:
        f.write("random-cut-polytope\n")
        f.write("H-representation\n")
        f.write("begin\n")
        f.write(f"{len(A)} {len(A[0]) + 1} rational\n")
        for row, rhs in zip(A, b):
            f.write(" ".join([str(rhs)] + [str(-v) for v in row]) + "\n")
        f.write("end\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ine", metavar="FILE",
                        help="also write the system in lrs .ine format")
    args = parser.parse_args()

    A, b = build_polytope()
    print(f"{lrs_version()}, polytope: dim={DIM}, {len(A)} inequalities "
          f"({CUTS} cuts + {2 * DIM} box)")
    if args.ine:
        write_ine(args.ine, A, b)
        print(f"wrote {args.ine}")

    t0 = time.perf_counter()
    vertices = enumerate_vertices(A, b)
    elapsed = time.perf_counter() - t0

    # order-independent figures for cross-implementation comparison
    print(f"{vertices.shape[0]} vertices in {elapsed:.2f} s")
    print(f"checksum: coordinate sum = {vertices.sum():.6f}, "
          f"abs sum = {abs(vertices).sum():.6f}")


if __name__ == "__main__":
    main()
