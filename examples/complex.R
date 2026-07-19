## Benchmark polytope for the R package vertexenum -- the exact counterpart of
## examples/complex.py in this repository.
##
## Builds the identical polytope (the box [-2, 2]^10 intersected with 28
## pseudo-random integer cutting planes), enumerates its vertices with the R
## package, verifies the result against the checksums produced by the Python
## implementation, and reports the elapsed time for benchmarking.
##
## Requirements: the archived CRAN package `vertexenum` (and its dependency
## `numbers`). Install one of these ways (Rtools needed on Windows):
##   install.packages("numbers")
##   install.packages("r_package", repos = NULL, type = "source")   # from repo root
## or from the CRAN archive:
##   install.packages("https://cran.r-project.org/src/contrib/Archive/vertexenum/vertexenum_1.0.2.tar.gz",
##                    repos = NULL, type = "source")
##
## Run:  Rscript examples/complex.R

suppressPackageStartupMessages(library(vertexenum))

DIM <- 10L        # number of variables
CUTS <- 28L       # number of pseudo-random cutting planes
BOX <- 2L         # box is [-BOX, BOX]^DIM
SEED <- 20260719

## Reference values from examples/complex.py (Python + lrslib v5.1):
EXPECTED_VERTICES <- 40317L
EXPECTED_COORD_SUM <- -8522.372864
EXPECTED_ABS_SUM <- 678450.682632
CHECKSUM_TOL <- 1e-3

## Deterministic coefficient stream, identical to complex.py:
##   x_{k+1} = (1103515245 * x_k + 12345) mod 2^31,  x_0 = SEED
##   value_k = (x_{k+1} mod 1999) - 999
## R has no 64-bit integers and 1103515245 * x overflows the 2^53 exact-double
## range, so the multiplication is split into 16-bit halves; every
## intermediate stays below 2^48 and the arithmetic remains exact.
lcg_ints <- function(count, seed = SEED) {
  MULT <- 1103515245
  INC <- 12345
  MOD <- 2^31
  x <- seed
  out <- integer(count)
  for (k in seq_len(count)) {
    hi <- x %/% 65536            # < 2^15
    lo <- x %% 65536             # < 2^16
    x <- ((MULT * hi) %% MOD * 65536 + MULT * lo + INC) %% MOD
    out[k] <- x %% 1999 - 999
  }
  out
}

## Same construction as complex.py: CUTS rows of LCG coefficients (row-major)
## with b_i = sum_j |a_ij|, then the box constraints +-x_j <= BOX.
build_polytope <- function(dim = DIM, cuts = CUTS, box = BOX) {
  vals <- lcg_ints(cuts * dim)
  m <- cuts + 2L * dim
  A <- matrix(0, nrow = m, ncol = dim)
  b <- numeric(m)
  for (i in seq_len(cuts)) {
    row <- vals[((i - 1L) * dim + 1L):(i * dim)]
    if (all(row == 0)) row[1L] <- 1
    A[i, ] <- row
    b[i] <- sum(abs(row))
  }
  r <- cuts
  for (j in seq_len(dim)) {
    for (s in c(1, -1)) {
      r <- r + 1L
      A[r, j] <- s
      b[r] <- box
    }
  }
  list(A = A, b = b)
}

check <- function(label, ok, detail) {
  cat(sprintf("%-14s %s  (%s)\n", label, if (ok) "OK" else "FAIL", detail))
  ok
}

p <- build_polytope()
cat(sprintf("R %s, vertexenum %s, polytope: dim=%d, %d inequalities (%d cuts + %d box)\n",
            getRversion(), as.character(packageVersion("vertexenum")),
            DIM, nrow(p$A), CUTS, 2L * DIM))

timing <- system.time(V <- enumerate.vertices(p$A, p$b))
elapsed <- timing[["elapsed"]]

n_vertices <- nrow(V)
coord_sum <- sum(V)
abs_sum <- sum(abs(V))

cat(sprintf("%d vertices in %.2f s\n", n_vertices, elapsed))
cat(sprintf("checksum: coordinate sum = %.6f, abs sum = %.6f\n",
            coord_sum, abs_sum))

all_ok <- all(
  check("vertex count", n_vertices == EXPECTED_VERTICES,
        sprintf("%d vs %d", n_vertices, EXPECTED_VERTICES)),
  check("coord sum", abs(coord_sum - EXPECTED_COORD_SUM) < CHECKSUM_TOL,
        sprintf("%.6f vs %.6f", coord_sum, EXPECTED_COORD_SUM)),
  check("abs sum", abs(abs_sum - EXPECTED_ABS_SUM) < CHECKSUM_TOL,
        sprintf("%.6f vs %.6f", abs_sum, EXPECTED_ABS_SUM))
)

cat(sprintf("benchmark: R/vertexenum elapsed = %.2f s (Python/lrslib v5.1 reference: ~8 s)\n",
            elapsed))

if (!all_ok) {
  cat("checksum mismatch: R result differs from the Python reference\n")
  quit(status = 1L)
}
