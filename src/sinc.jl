struct SincBasis <: AbstractBasis end

function kinetic(dvr_obj::DVR{T, SincBasis}) where {T}
    xg = dvr_obj.xg
    dx = dvr_obj.dx
    ng = length(xg)
    
    tm = zeros(T, ng, ng)
    for i in 1:ng
        for j in 1:ng
            if i == j
                tm[i, j] = π^2 / 3 / dx^2
            else
                tm[i, j] = (-1)^(i-j) / (i-j)^2 / dx^2 * 2.0
            end
        end
    end
    return tm / 2
end
