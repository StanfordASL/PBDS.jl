struct SO3 <: Manifold end

dim(::Type{SO3}) = 3
embdim(::Type{SO3}) = 9

abstract type ExpMap <: ChartID end
isglobal(::Chart{<:ExpMap,SO3}) = false
struct ExpMap0 <: ExpMap end
struct ExpMapX <: ExpMap end
struct ExpMapY <: ExpMap end
struct ExpMapZ <: ExpMap end

skew(w) = SA[0.   -w[3]  w[2]
             w[3]  0.   -w[1]
            -w[2]  w[1]  0.]
unskew(ŵ) = SA[ŵ[3,2], ŵ[1,3], ŵ[2,1]]

function logSO3(R)
    trR = tr(R)
    if trR < 3
        θ = acos((trR-1)/2)
        # ŵ = (R - R')/(2*sin(θ))
        # w = unskew(ŵ)
        w = SA[R[3,2] - R[2,3], R[1,3] - R[3,1], R[2,1] - R[1,2]]/(2*sin(θ))
        return θ*w
    else
        return @SVector zeros(3)
    end
end

function expSO3(w)
    w == (@SVector zeros(3)) && (@show "expSO3 singularity"; return SMatrix{3,3,eltype(w)}(I)) 
    ŵ = skew(w)
    wnorm = norm(w)
    SMatrix{3,3,eltype(w)}(I) + ŵ*sin(wnorm)/wnorm + ŵ^2*(1-cos(wnorm))/wnorm^2
end

tr_vec(pe::SVector{9,S}) where S = pe[1] + pe[5] + pe[9]
chart_boundary_test(pe, ::Chart{ExpMap0,SO3}) = tr_vec(pe)
chart_boundary_test(pe, ::Chart{ExpMapX,SO3}) = tr(RotX(π/2)*reshape(pe, Size(3,3)))
chart_boundary_test(pe, ::Chart{ExpMapY,SO3}) = tr(RotY(π/2)*reshape(pe, Size(3,3)))
chart_boundary_test(pe, ::Chart{ExpMapZ,SO3}) = tr(RotZ(π/2)*reshape(pe, Size(3,3)))

chart_to_emb(p0, ::Chart{ExpMap0,SO3}) = reshape(expSO3(p0), Size(9))
chart_to_emb(px, ::Chart{ExpMapX,SO3}) = reshape(RotX(-π/2)*expSO3(px), Size(9))
chart_to_emb(py, ::Chart{ExpMapY,SO3}) = reshape(RotY(-π/2)*expSO3(py), Size(9))
chart_to_emb(pz, ::Chart{ExpMapZ,SO3}) = reshape(RotZ(-π/2)*expSO3(pz), Size(9))

emb_to_chart(pe, ::Chart{ExpMap0,SO3}) = logSO3(reshape(pe, Size(3,3)))
emb_to_chart(pe, ::Chart{ExpMapX,SO3}) = logSO3(RotX(π/2)*reshape(pe, Size(3,3)))
emb_to_chart(pe, ::Chart{ExpMapY,SO3}) = logSO3(RotY(π/2)*reshape(pe, Size(3,3)))
emb_to_chart(pe, ::Chart{ExpMapZ,SO3}) = logSO3(RotZ(π/2)*reshape(pe, Size(3,3)))

function zero_jacobian_component(::Chart{ExpMap0,SO3})
    SA[0 0  0  0 0 1 0 -1 0
       0 0 -1  0 0 0 1  0 0
       0 1  0 -1 0 0 0  0 0.]
end

function chart_to_emb_jacobian(p, C::Chart{ExpMap0,SO3})
    if norm(p - (@SVector zeros(3))) < 1e-12
        return zero_jacobian_component(C)'
    else
        return ForwardDiff.jacobian(p -> chart_to_emb(p, C), p)
    end
end

function emb_to_chart_jacobian(pe, C::Chart{ExpMap0,SO3})
    if (chart_boundary_test(pe, C) - 3) < 1e-12
        return zero_jacobian_component(C)/2
    else
        return ForwardDiff.jacobian(pe -> emb_to_chart(pe, C), pe)
    end
end

chart_choice_coord_rep(::Chart{<:ExpMap,SO3}) = EmbRep()
function choose_chart_emb(::EmbRep, pe, Cn::Chart{<:ExpMap,SO3})
    C0, Cx = Chart{ExpMap0,SO3}(), Chart{ExpMapX,SO3}()
    Cy, Cz = Chart{ExpMapY,SO3}(), Chart{ExpMapZ,SO3}()

    lo_bound, hi_bound = -0.5, 2.5
    val = chart_boundary_test(pe, Cn)
    if val < lo_bound || val > hi_bound
        val0 = chart_boundary_test(pe, C0)
        valx = chart_boundary_test(pe, Cx)
        valy = chart_boundary_test(pe, Cy)
        valz = chart_boundary_test(pe, Cz)

        val0 = val0 > hi_bound ? -1. : val0
        valx = valx > hi_bound ? -1. : valx
        valy = valy > hi_bound ? -1. : valy
        valz = valz > hi_bound ? -1. : valz

        val_max, ind_max = findmax(SA[val0, valx, valy, valz])
        return (C0, Cx, Cy, Cz)[ind_max]
    else
        return Cn
    end
end

function apply_rotation_emb(xme, R)
    Rtarget = reshape(xme, Size(3,3))
    reshape(R*Rtarget, Size(9))
end

function apply_inverse_rotation_jacobian_emb(xme)
    R = reshape(xme, Size(3,3))
    Rinv = R'
    ForwardDiff.jacobian(xme -> apply_rotation_emb(xme, Rinv), xme)
end

function apply_rotation_chart(xm, R, Ci::Chart{<:ExpMap,SO3})
    xme = chart_to_emb(xm, Ci)
    xme_out = apply_rotation_emb(xme, R)
    xm_out = emb_to_chart(xme_out, Ci)
end

function apply_rotation_chart_emb(xm, R, Ci::Chart{<:ExpMap,SO3})
    xme = chart_to_emb(xm, Ci)
    xme_out = apply_rotation_emb(xme, R)
end

function apply_inverse_rotation_jacobian_chart(xm, Ci::Chart{<:ExpMap,SO3})
    xme = chart_to_emb(xm, Ci)
    R = reshape(xme, Size(3,3))
    Rinv = R'
    xme_out = reshape(SMatrix{3,3,Float64}(I), Size(9))
    J1 = ForwardDiff.jacobian(xm -> apply_rotation_chart_emb(xm, Rinv, Ci), xm)
    J2 = emb_to_chart_jacobian(xme_out, Ci)
    J2*J1
end

function default_metric(pn, task::Task{<:TaskMap{M,SO3}}, CN::Chart{<:ExpMap,SO3}) where M
    C0 = Chart{ExpMap0,SO3}()
    JL = apply_inverse_rotation_jacobian_chart(pn, CN)
    Jϕ = chart_transition_jacobian(pn, CN, C0)
    J = Jϕ*JL
    J'*SMatrix{3,3,Float64}(I)*J
end

default_weight_metric(pn, wn, task::Task{<:TaskMap{M,SO3}}, CN::Chart{<:ExpMap,SO3}) where M =
    default_metric(pn, task, CN)