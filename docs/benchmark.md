# Benchmark

Running complex.py and complex.R to compare performance of the Python and R implementations for the same polytope.

## Python version
```python

lrslib v.5.1 2015.1.28, polytope: dim=10, 48 inequalities (28 cuts + 20 box)
40317 vertices in 7.90 s
checksum: coordinate sum = -8522.372864, abs sum = 678450.682632

Process finished with exit code 0
```

## R version
- requires RTools
- requires 'numbers' CRAN library
  - ```install.packages("numbers")  # dependency, still on CRAN```
- requires vertexenum package (no longer on CRAN)
  - ```install.packages("vertexenum_1.0.2.tar.gz", repos = NULL, type = "source")```

```r
C:\dev\python\vertexenum\examples>Rscript complex.R
R 4.6.1, vertexenum 1.0.2, polytope: dim=10, 48 inequalities (28 cuts + 20 box)
40317 vertices in 10.98 s
checksum: coordinate sum = -8522.372864, abs sum = 678450.682632
vertex count   OK  (40317 vs 40317)
coord sum      OK  (-8522.372864 vs -8522.372864)
abs sum        OK  (678450.682632 vs 678450.682632)
benchmark: R/vertexenum elapsed = 10.98 s (Python/lrslib v5.1 reference: ~8 s)
```

## Conclusion
- Python implementation is faster than R implementation for this benchmark (7.90 s vs 10.98 s).