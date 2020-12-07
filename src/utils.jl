# Approaches v/norm(v) for large v, but approaches zero smoothly as v -> 0
function smooth_normalization(v, α=1.e3)
    h = norm(v) + 1/α * log(1 + exp(-2*α*norm(v)))
    return v/h
end

# Smooth approximation of Euclidean norm
function smooth_norm(x, α=1.e-3)
    sqrt(sum(x.^2) + α^2) - α
end

# 1D Barrier functions
# Finite exponential barrier with value specification
# Inputs:
#   x : Eval location
#   x_lim : Barrier location 
#   x_ctl : Control location
#   f_lim : Value at x_lim
#   ctl_ratio : Control ratio (value at x_ctl = f_lim*ctl_ratio)
function exp_barrier_value(x, x_lim, x_ctl, f_lim, ctl_rat)
    σ = (x_ctl - x_lim)/log(1/ctl_rat)
    f_lim*exp(-(x - x_lim)/σ)
end

# Finite exponential barrier with slope specification
# Works for both min and max barrier
# Inputs:
#   x : Eval location
#   x_lim : Barrier location 
#   x_ctl : Control location
#   fdot_lim : Slope at x_lim
#   ctl_ratio : Control ratio (slope at x_γ = fdot_lim*ctl_ratio)
function exp_barrier_slope(x, x_lim, x_ctl, fdot_lim, ctl_ratio)
    σ = (x_lim - x_ctl)/log(1/ctl_ratio)
    fdot_lim*σ*exp((x - x_lim)/σ)
end

# Softplus with slope specification
# Inputs:
# Works for both min and max barrier
#   x : Eval location
#   x_lim : Barrier location
#   x_ctl : Control location
#   fdot_lim : Slope at x_lim
#   ctl_ratio : Control ratio (slope at x_ctl = fdot_lim*ctl_ratio)
function softplus_slope(x, x_lim, x_ctl, fdot_lim, γ)
    σ = (x_ctl - x_lim)*log(1 - 2*ctl_ratio)
    2*γ*fdot_lim*σ*log(1 + exp((x - x_lim)/σ))
end

# Smooth approximation of max and min functions
smooth_max(x::AbstractVector, α=1.e2) = sum(@. x*exp(α*x))/sum(@. exp(α*x))
smooth_max(x, y, α=1.e2) = (x*exp(α*x) + y*exp(α*y))/(exp(α*x) + exp(α*y))

smooth_min(x::AbstractVector, α=1.e2) = -smooth_max(-x, α)
smooth_min(x, y, α=1.e2) = -smooth_max(-x, -y, α)

# Smooth approximation of clamp function
# Inputs:
#  x : Eval location
#  x_trans : Clamp transition start point
#  x_lim : Clamp transition ending point
#  slope_lim : Slope at x_lim
function smooth_clamp(x, x_trans, x_lim, slope_lim=1.e-3)
    if x < -x_trans
        c = 2. - 4. / slope_lim
        d = (-c + sqrt(c^2 - 4.))/2.
        k = -log(d)/(-x_lim + x_trans)
        L = 4. / k
        h = -x_trans - 0.5*L
        
        return L/(1+exp(-k*(x + x_trans))) + h
        
    elseif -x_trans <= x < x_trans
        return x
        
    else
        c = 2. - 4. / slope_lim
        d = (-c + sqrt(c^2 - 4.))/2.
        k = -log(d)/(x_lim - x_trans)
        L = 4. / k
        h = x_trans - 0.5*L

        return L/(1+exp(-k*(x - x_trans))) + h
    end
end

# Smooth absolute value function
function smooth_abs(x, ε=1.e-3)
    sqrt(x^2 + ε)
end

# Logistic function with value specification
# Inputs:
#  x : Eval location
#  xmin : Location of specified near-minimum value
#  xmax : Location of specified near-maximum value
#  fmin : Minimum value at limit
#  fmax : Maximum value at limit
#  ε : Gap from fmin/fmax at xmin/xmax
function logistic_val(x, xmin, xmax, fmin, fmax, ε)
    h = fmin
    L = fmax - fmin
    k = log(((fmin-fmax+ε)/ε)^2)/(xmax - xmin)
    x0 = 1. / k * log((fmax-fmin)/ε - 1.) + xmin
    
    L/(1. + exp(-k*(x-x0))) + h
end

# Angle wrapping functions
wrapTo2Pi(x) = mod(2*π + mod(x, 2*π), 2*π)
wrapToPi(x) = -π + wrapTo2Pi(x + π)

# Distance to line segment
# Inputs:
#  x : Eval location
#  p1 : Endpoint 1 of line segment
#  p2 : Endpoint 2 of line segment
function line_segment_distance(x, p1, p2)
    p2p1 = p2 - p1
    xp1 = x - p1
    xp2 = x - p2
    if dot(xp1, p2p1) > 0.
        if dot(xp2, p2p1) < 0.
            normal_vec = xp1 - p2p1*dot(xp1, p2p1)/norm(p2p1)^2
            return norm(normal_vec)
        else
            return norm(xp2)
        end
    else
        return smooth_norm(xp1)
    end
end

# Distance from point with positive x to box bordering y axis
# Inputs:
#  x : Eval location
#  p_lo : Low corner location
#  p_hi : High corner location
#  α : smoothing parameter
function aligned_box_distance(x, p_lo, p_hi, α=1.e2)
    if x[1] <= p_lo[1]
        if p_lo[2] <= x[2] <= p_hi[2]
            return 0.
        else
            return smooth_max(smooth_max(0., p_lo[2] - x[2], α), smooth_max(0., x[2] - p_hi[2], α), α)
        end
    else
        return line_segment_distance(x, p_lo, p_hi)
    end
end

## Plotting utiliies
function html_video(filename)
    base64_video = Base64.base64encode(open(read, filename))
    """<video controls src="data:video/x-m4v;base64,$base64_video">"""
end