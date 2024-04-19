function LotkaVolterraConst(alpha, beta, delta, gamma)

    # X     --> 2X (alpha)
    # X + Y --> Y  (beta - delta)
    # X + Y --> 2Y (delta)
    # Y     --> \o (gamma)

    N = 2;
    M = 4;

    g = [2.0, 2.0];
    zeta = [[1, 0], [-1, 0], [-1, 1], [0, -1]];

    function propensities(num)
        #return [alpha*num[1], beta*num[1]*num[2], delta*num[1]*num[2], gamma*num[2]];
        return [alpha*num[1], (beta - delta)*num[1]*num[2], delta*num[1]*num[2], gamma*num[2]];
    end

    return N, M, g, zeta, propensities;
end