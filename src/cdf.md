Use `cumulativefunc` for normalization integrals.
Instead of one integral it computes two now, `CF(x2)-CF(x1)`, it does not impact the time significantly.

The new benchmark:
```
31892.662783336527
  197.800 μs (10 allocations: 158.41 KiB)  # @btime nll(1.1)
  201.700 μs (19 allocations: 158.72 KiB)  # @btime nll(1.1, [2.4, 1.3])
  198.600 μs (10 allocations: 158.41 KiB)  # @btime nll(1.1; p=(μ=2.4, σ=1.3))
  198.700 μs (10 allocations: 158.41 KiB)  # @btime nll(1.1; p=(σ=1.3, μ=2.4))
```
to be compared to the previous one
```
31892.66278333851
  191.500 μs (7 allocations: 156.72 KiB)   # @btime nll(1.1)
  196.000 μs (16 allocations: 157.03 KiB)  # @btime nll(1.1, [2.4, 1.3])
  194.800 μs (7 allocations: 156.72 KiB)   # @btime nll(1.1; p=(μ=2.4, σ=1.3))
  194.400 μs (7 allocations: 156.72 KiB)   # @btime nll(1.1; p=(σ=1.3, μ=2.4))
```


It allows dispatch on the function type. Usage of the analytic integrals speed up the calculations drastically.
```
  73.200 μs  # @btime gauss.(randn(1_000))
```
to be compared to the numerical evaluation of the integral
```
  6.978 ms  # @btime gauss.(randn(1_000))
```

