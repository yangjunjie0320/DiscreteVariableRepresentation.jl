struct Sinc <: AbstractBasis
    xg::Vector
    dx::Real
    ng::Int

    function Sinc(xg::Vector{T}) where {T}
        ng = length(xg)
        dx = xg[2] - xg[1]

        for i in 2:ng
            xi = xg[i]; xj = xg[i-1]
            @assert abs(xi - xj - dx) < 1e-10
        end

        new(xg, dx, ng)
    end
end

function kinet(dvr_obj::DVR{T, Sinc, V}) where {T, V}
    xg = dvr_obj.b.xg
    dx = dvr_obj.b.dx
    ng = dvr_obj.b.ng
    pi2 = pi^2
    
    tm = zeros(T, ng, ng)
    for i in 1:ng
        for j in 1:ng
            imj = i - j
            tm[i, j] = i == j ? pi2 / 3 : (-1)^imj / imj^2 * 2.0
        end
    end

    return tm / 2 / dx^2
end
