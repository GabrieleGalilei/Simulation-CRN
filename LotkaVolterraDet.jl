function LotkaVolterraDet(c, num0, time0, time_max)
    function f(u, c, t)
        dx = c[1]*u[1] - c[2]*u[1]*u[2];
        dy = c[3]*u[1]*u[2] - c[4]*u[2];
        return [dx, dy];
    end
    u0 = num0;
    tspan = (time0, time_max);
    prob = ODEProblem(f, u0, tspan, c);
    sol = solve(prob);

    return sol;
end