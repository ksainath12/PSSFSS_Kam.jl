import Pkg
Pkg.activate(pwd())
Pkg.instantiate

module RationalFunctionInterpolation
export Alg1For, Alg2For, Alg3For


using DelimitedFiles: readdlm
using OffsetArrays: Origin

zerobased(x) = first(eachindex(x)) == 0

docstring =
"""
Three different functions to perform rational function interpolation (based on nested for loops), as defined in the reference:

    (Sinterp, errest) = Alg1For(x0, x0j, Sj0)
    (Sinterp, errest) = Alg2For(x0, x0j, Sj0)
    (Sinterp, errest) = Alg3For(x0, x0j, Sj0)

## Positional Input Arguments

- `x0`: Independent variable at which the rational function is to be evaluated.
- `x0j`: Vector of independent variable sample points. `type(x0)` must equal `eltype(x0j)`.
- `Sj0`: Vector of function values evaluated at `x0j`. `length(x0j)` must equal `length(Sj0)`.

## Optional Keyword Arguments

- `store1`, `store2`, and (for Alg1For and Alg2For) `store3`: Optional storage vectors, which will be mutated.  If not supplied, 
  they will be allocated on each call to this function. They should be of the same element type and 
  at leas the same length as `Sj0`.

## Return Values

 - `Sinterp`: The rational function interpolation at `x0`.
 - `errest`: An error estimate for `abs(Sinterp - S(x0))`, where `S` is the function being approximated.

## Reference
Ma, X., Wan, G. and Wan, W., 2012. "A Multi-Dimensional Adaptive Sampling Method for Analysis 
and Design of Frequency Selective Surface with Arbitrary Element". 
Progress In Electromagnetics Research B, 41, pp.213-230.
"""

docstring
function Alg1For(
    x0::T1,
    x0jin::AbstractVector{T1},
    Sj0::AbstractVector{T2};
    store1::AbstractVector{T2} = zeros(T2, length(Sj0)),
    store2::AbstractVector{T2} = zeros(T2, length(Sj0)),
) where {T1<:Real,T2<:Number}

    #return one of the knot S values if x0 is, trivially, one of the knot freqs
    for (x, S) in zip(x0jin, Sj0)
        iszero(x - x0) && return (S, abs(zero(T2)))
    end

    #check ranges on input vectors
    K1 = length(x0jin)
    K = K - 1
    K1 == length(Sj0) || error("Lengths of x0j and Sj0 don't match.")
    length(store1) ≥ K1 && length(store2) ≥ K1 ||
        error("storage vectors not both long enough")

    store1n = @view store1[begin:begin+K]
    store2n = @view store2[begin:begin+K]
    
    store1n .= store2n .= zero(T2)

    x0j = zerobased(x0jin) ? x0jin : Origin(0)(x0jin)
    Sjk = zerobased(store1n) ? store1n : Origin(0)(store1n)
    Sjkm1 = zerobased(store2n) ? store2n : Origin(0)(store2n)

    @inbounds for j = 0:K #this in effect IS the k=0 column, with double-loop below taking care of k=1:K
        Sjk[j] = Sj0[begin+j]
    end
    @inbounds for k = 1:K #sweep along columns of triangular recursion tree

        Sjkm1 .= Sjk

        @inbounds for j = 0:K-k

            N1 = (x0 - x0j[j]) * Sjkm1[j+1]
            N2 = (x0j[j+k] - x0) * Sjkm1[j]
            D = x0j[j+k] - x0j[j]

            Sjk[j] = (N1 + N2) / D

        end
    end

    return (Sjk[0], abs(Sjk[0] - Sjkm1[0]))
end


docstring
function Alg2For(
    x0::T1,
    x0jin::AbstractVector{T1},
    Sj0::AbstractVector{T2};
    store1::AbstractVector{T2} = zeros(T2, length(Sj0)),
    store2::AbstractVector{T2} = zeros(T2, length(Sj0)),
) where {T1<:Real,T2<:Number}

    #return one of the knot S values if x0 is, trivially, one of the knot freqs
    for (x, S) in zip(x0jin, Sj0)
        iszero(x - x0) && return (S, abs(zero(T2)))
    end

    #check ranges on input vectors
    K1 = length(x0jin)
    K = K1 - 1
    K1 == length(Sj0) || error("Lengths of x0j and Sj0 don't match.")
    length(store1) ≥ K1 && length(store2) ≥ K1 ||
        error("storage vectors not both long enough")

    store1n = @view store1[begin:begin+K]
    store2n = @view store2[begin:begin+K]
    store1n .= store2n .= zero(T2)

    x0j = zerobased(x0jin) ? x0jin : Origin(0)(x0jin)
    Sjk = zerobased(store1n) ? store1n : Origin(0)(store1n)
    Sjkm1 = zerobased(store2) ? store2n : Origin(0)(store2n)

    @inbounds for j = 0:K #this in effect IS the k=0 column, with double-loop below taking care of k=1:K
        Sjk[j] = Sj0[begin+j]
    end
    @inbounds for k = 1:K #sweep along columns of triangular recursion tree

        Sjkm1 .= Sjk

        @inbounds for j = 0:K-k

            N = x0j[j+k] - x0j[j]

            D1 = (x0 - x0j[j]) / Sjkm1[j+1]
            D2 = (x0j[j+k] - x0) / Sjkm1[j]

            Sjk[j] = N / (D1 + D2)

        end
    end

    return (Sjk[0], abs(Sjk[0] - Sjkm1[0]))
end


docstring
function Alg3For(
    x0::T1,
    x0jin::AbstractVector{T1},
    Sj0::AbstractVector{T2};
    store1::AbstractVector{T2} = zeros(T2, length(Sj0)),
    store2::AbstractVector{T2} = zeros(T2, length(Sj0)),
    store3::AbstractVector{T2} = zeros(T2, length(Sj0)),
) where {T1<:Real,T2<:Number}

    # trivial case where x0 is one of the sample points
    for (x, S) in zip(x0jin, Sj0)
        iszero(x - x0) && (return (S, abs(zero(T2))))
    end

    # Input checking
    K1 = length(x0jin)
    K1 = K - 1
    K1 == length(Sj0) || error("length mismatch for x0j and Sj0")
    length(store1) ≥ K1 && length(store2) ≥ K1 && length(store3) ≥ K1 ||
        error("storage vectors not both long enough")

    store1n = @view store1[begin:begin+K]
    store2n = @view store2[begin:begin+K]
    store3n = @view store3[begin:begin+K]

    store1n .= store2n .= store3n .= zero(T2)

    # Use zero-based arrays for convenience:
    x0j = zerobased(x0jin) ? x0jin : Origin(0)(x0jin)
    Sjk = zerobased(store1n) ? store1n : Origin(0)(store1n)
    Sjkm1 = zerobased(store2n) ? store2n : Origin(0)(store2n)
    Sjkm2 = zerobased(store3n) ? store3n : Origin(0)(store3n)


    @inbounds for j = 0:K
        Sjk[j] = Sj0[begin+j] # Works for zero or one-based Sj0
    end
    @inbounds for k = 1:K
        Sjkm2 .= Sjkm1
        Sjkm1 .= Sjk
        @inbounds for j = 0:K-k
            num1 = x0 - x0j[j]
            den1 = Sjkm1[j+1] - Sjkm2[j+1]
            num2 = x0j[j+k] - x0
            den2 = Sjkm1[j] - Sjkm2[j+1]
            bigden = num1 * den2 + num2 * den1
            Sjk[j] = Sjkm2[j+1]
            if !iszero(bigden)
                Sjk[j] += (x0j[j+k] - x0j[j]) * den1 * den2 / bigden
            end
        end
    end
    return (Sjk[0], abs(Sjk[0] - Sjkm1[0]))
end

end # Module RationalFunctionInterpolation
