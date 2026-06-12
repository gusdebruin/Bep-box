function solution(geometry::Geometry.AbstractGeometry)
    ϵ > 0.0 || throw(ArgumentError("ϵ should be greater than 0."))

    s(p) = ϵ / ((p - 0.5)^2 + ϵ)
    u(x, y) = 16 * s(4 * x * y) * x * (1 - x) * y * (1 - y)
    u_expr(x) = [@. u(x[:, 1], x[:, 2])]

    return Forms.AnalyticalFormField(0, u_expr, geometry, "u")
end
