function forcing(geometry::Geometry.AbstractGeometry)
    ϵ > 0.0 || throw(ArgumentError("ϵ should be greater than 0."))

    dx2(x, y) = -(128 * ϵ * (y - 1) * y * (16 * y * (4 * (y * x * (x * (16 * y * (4 * y - 1)
    * x + 12 * ϵ + 3) - 12 * ϵ - 3) + ϵ) + 1) - (4 * ϵ + 1)^2)) / (16 * y * x * (4 * y * x -
    1) + 4 * ϵ + 1)^3

    dy2(x, y) = -(128 * ϵ * (x - 1) * x * (16 * x * (4 * (x * y * (y * (16 * x * (4 * x - 1)
    * y + 12 * ϵ + 3) - 12 * ϵ - 3) + ϵ) + 1) - (4 * ϵ + 1)^2)) / (16 * x * y * (4 * x * y -
    1) + 4 * ϵ + 1)^3

    f_expr(x) = [@. -(dx2(x[:, 1], x[:, 2]) + dy2(x[:, 1], x[:, 2]))]

    return Forms.AnalyticalFormField(0, f_expr, geometry, "f")
end
