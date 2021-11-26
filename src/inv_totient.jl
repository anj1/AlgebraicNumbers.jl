function all_divisors(factors::Nemo.Fac{T}) where T <: Integer
    divisors = []

    function generate_divisors(curIndex::Int, curDivisor::T, factors::Nemo.Fac)
        result = iterate(factors, curIndex)
        if isnothing(result)
            push!(divisors, curDivisor)
            return
        else
            (pf, nextIndex) = result
        end

        for i = 0:pf.second
            generate_divisors(nextIndex, curDivisor, factors);
            curDivisor *= pf.first
        end
    end

    generate_divisors(1, T(1), factors)
    divisors
end 

# Find all inverses of the Euler Totient function.
# This implements the algorithm of Contini, Croot, and Shparlinski (2006)
function inv_totient(x::T) where T <: Integer
    invs = Set{T}()

    function totient_reps(x::T, factor_list::Vector{Tuple{T,Int}})
        if x<1
            return
        end 

        if x==1
            m = prod([p[1]^(p[2]+1) for p in factor_list[2:end]])
            push!(invs, m)
        end

        # TODO: Don't need to recalcuate this each time.
        factors = Nemo.factor(x)

        # Find all divisors of x of the form (p^0)*(p-1), where p is a prime
        for divisor in all_divisors(factors)
            d = divisor+1
            if Nemo.isprime(d) && d > factor_list[end][1]
                pair = (d, 0)
                totient_reps(div(x, divisor), cat(factor_list, pair, dims=1))
            end
        end

        # Find all divisors of x of the form (p^γ)*(p-1), γ>=1, where p is a prime
        for factor in factors
            for γ = 1:factor.second
                d = (factor.first^γ)*(factor.first-1)
                if divrem(x, d)[2]==0 && factor.first > factor_list[end][1]
                    pair = (factor.first, γ)
                    totient_reps(div(x, d), cat(factor_list, pair, dims=1))
                end
            end 
        end
    end 

    totient_reps(x, [(T(0),0)])

    return invs
end