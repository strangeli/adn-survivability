using PowerDynamics: @DynamicNode

@DynamicNode ExtDampedGrid(H, D, Γ, Ω, V) begin
    @assert D >= 0 "damping (D) should be >=0"
    @assert H > 0 "inertia (H) should be >0"
    @assert abs(V) >= 0 "slack voltage (V) should be >=0"
    Ω_H = Ω / (2*H)
end [[ω, dω]] begin
    v = abs(u)
    dv = - Γ * (v - V)
    dω = - Ω_H * D * ω
    du = u * (im * ω + dv / v)
end
