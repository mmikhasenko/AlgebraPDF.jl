"""
    nt(s::Symbol, v = 0.0)

Creates a named tuple `(s=v,)` where `s` is a provided symbol, and `v` is the value.
"""
nt(s::Symbol, v = 0.0) = NamedTuple{(s,)}([v])

inrange(x, r) = r[1] ≤ x ≤ r[2]
