function TauLeapPostCheckNano2(
    sz::Int64,
    eps::Float64,
    p::Float64,
    pstar::Float64,
    q::Float64,
    V,
    mins::Int16,
    nuc::Float64,
    grw::Float64,
    repetition::Int16,
    slope::Float64,
    intercept::Float64,
    threshold::Int64 = -10
    )

    sequenza = Int32.(0:(mins-1));

    if threshold < 0
        threshold = ceil(Int64, 1.5 * sqrt(V));
    end
    comptimes::Array{Float64} = zeros(Float64, repetition);

    if slope > 2.6 #Alternative
        time25 = 0.2 * 4.49;
        time5 = 0.4 * 4.49;
        time75 = 0.6 * 4.49;
        time10 = 0.8 * 4.49;
        #threshold = V+1;
    else #Classical
        time25 = 0.2 * 4.49;
        time5 = 0.4 * 4.49;
        time75 = 0.6 * 4.49;
        time10 = 0.8 * 4.49;
        threshold = ceil(Int64, 2 * sqrt(sqrt(V)));
    end

    taus::Matrix{Float64} = zeros(Float64, repetition, 100000);
    out::Matrix{Int64} = zeros(Int64, repetition, threshold);
    rec25::Matrix{Int64} = zeros(Int64, threshold, repetition);
    rec5::Matrix{Int64} = zeros(Int64, threshold, repetition);
    rec75::Matrix{Int64} = zeros(Int64, threshold, repetition);
    rec10::Matrix{Int64} = zeros(Int64, threshold, repetition);

    Threads.@threads for i = 1:repetition
        tick()
        num::Array{Int64} = zeros(Int64, threshold);
        monomeri::Int32 = V;
        rateNUC::Float64 = 0.0;
        ratesGR::Array{Int64} = zeros(Int64, threshold);
        largest::Int64 = mins+1;
        taken25::Bool = false;
        taken5::Bool = false;
        taken75::Bool = false;
        taken10::Bool = false;
        log_num::Float64 = 0.0;
        log_den::Float64 = 0.0;
        r::Float64 = 0.0;
        accepted::Bool = true;
        linea::Int64 = 0;
        temp::Float64 = 0.0;
        tempmat::Matrix{Union{Float64, Int64}} = Matrix{Union{Float64, Int64}}(undef, (sz, 2))
        contatore = 0;

        newpar::Float64 = nuc/grw;

        time::Float64 = 0.0;

        S = [Matrix{Union{Float64, Int64}}(undef, (sz, 2)) for _ = 1:threshold];
        B::Array{Int64} = ones(Int64, threshold);
        T::Array{Float64} = zeros(Float64, threshold);
        C::Array{Int64} = zeros(Int64, threshold);

        global cont = 1;
        for k = 1:largest
            S[k][1,:] = [0.0, 0];
        end

        if monomeri >= mins
            rateNUC = newpar * exp(sum(log.(monomeri .- sequenza)));
        end
        for k = mins:largest-1
            ratesGR[k+1] = 0;
            if num[k+1] > 0
                ratesGR[k+1] = monomeri * num[k+1];
            end
        end

        # CALCOLO DI TAU
        mu::Array{Float64} = zeros(Float64, largest+1);
        sigmasq::Array{Float64} = zeros(Float64, largest+1);

        #nucleazione
        if rateNUC != 0
            mu[1] = ((-1)*mins) * rateNUC;
            sigmasq[1] = (mins^2) * rateNUC;
            mu[mins+1] = rateNUC;
            sigmasq[mins+1] = rateNUC;
        end
        #crescita
        for i = mins:largest-1
            if ratesGR[i+1] != 0
                mu[1] -= ratesGR[i+1];
                sigmasq[1] += ratesGR[i+1];
                mu[i+1] -= ratesGR[i+1];
                sigmasq[i+1] += ratesGR[i+1];
                mu[i+2] += ratesGR[i+1];
                sigmasq[i+2] += ratesGR[i+2];
            end
        end

        mu = grw.*mu;
        sigmasq = grw.*sigmasq;

        tau::Float64 = Inf;

        #monomeri
        temp = log(eps) + log(monomeri) - log(3 + (1/(monomeri-1)) + (2/(monomeri-2)));
        temp = maximum([temp, 0.0]);
        firstterm::Float64 = temp - log(abs(mu[1]));
        secondterm::Float64 = 2.0 * temp - log(sigmasq[1]);
        tau = minimum([firstterm, secondterm, tau]);
        #particelle
        for i = mins:largest-1
            if num[i+1] == 0
                firstterm = - log(abs(mu[i+1]));
                secondterm = - log(abs(sigmasq[i+1]));
                tau = minimum([firstterm, secondterm, tau]);
            else
                temp = log(eps) + log(num[i+1]) - log(2);
                temp = maximum([temp, 0.0]);
                firstterm = temp - log(abs(mu[i+1]));
                secondterm = 2.0 * temp - log(abs(sigmasq[i+1]));
                tau = minimum([firstterm, secondterm, tau]);
            end
        end
        tau = exp(tau);

        while monomeri > 0
            N::Array{Int64} = zeros(Int64, largest);
            row::Array{Int64} = Array{Int64}(undef, largest);
            N_nuc::Int64 = 0;

            rateNUC = 0.0;
            if monomeri >= mins
                rateNUC = newpar * exp(sum(log.(monomeri .- sequenza)));
            end

            for k = mins:(largest-1)
                ratesGR[k+1] = 0;
                if num[k+1] > 0
                    ratesGR[k+1] = monomeri * num[k+1];
                end
            end

            #nucleazione
            if rateNUC == 0.0
                N_nuc = 0;
                row[1] = 1;
            else
                if exp(log(grw) + log(rateNUC) + log(tau)) + T[1] >= S[1][B[1],1]
                    N_nuc = convert(Int64, rand(Poisson( exp(log(grw) + log(rateNUC) + log(tau)) + T[1] - S[1][B[1],1] )) + S[1][B[1],2] - C[1]);
                    row[1] = B[1];
                else
                    linea = 2;
                    while exp(log(grw) + log(rateNUC) + log(tau)) + T[1] >= S[1][linea,1] && linea < B[1]
                        linea += 1;
                    end
                    if S[1][linea,1] - S[1][linea-1,1] < 0.0 
                        N_nuc = 0;
                    else
                        log_num = log(exp( log(grw) + log(rateNUC) + log(tau)) + T[1] - S[1][linea-1, 1]);
                        log_den = log(S[1][linea,1] - S[1][linea-1,1]);
                        r = exp(log_num - log_den);
                        N_nuc = convert(Int64, rand(Binomial( S[1][linea,2] - S[1][linea-1,2] , r )) + S[1][linea-1,2] - C[1]);
                    end
                    row[1] = linea - 1;
                end
            end

            #crescita
            for k = mins:largest-1
                # definiamo row e N_k
                if ratesGR[k+1] == 0
                    N[k+1] = 0;
                    row[k+1] = 1;
                else
                    if exp(log(grw) + log(ratesGR[k+1]) + log(tau)) + T[k+1] >= S[k+1][B[k+1], 1]
                        N[k+1] = convert(Int64, rand(Poisson( exp(log(grw) + log(ratesGR[k+1]) + log(tau)) + T[k+1] - S[k+1][B[k+1],1])) + S[k+1][B[k+1],2] - C[k+1]);
                        row[k+1] = B[k+1];
                    else
                        linea = 2;
                        while exp(log(grw) + log(ratesGR[k+1]) + log(tau)) + T[k+1] >= S[k+1][linea,1] && linea < B[k+1]
                            linea += 1;
                        end
                        if S[k+1][linea,1] - S[k+1][linea-1,1] == 0.0
                            N[k+1] = 0;
                        else
                            log_num = log(exp(log(grw) + log(ratesGR[k+1]) + log(tau)) + T[k+1] - S[k+1][linea-1, 1]);
                            log_den = log(S[k+1][linea,1] - S[k+1][linea-1,1]);
                            r = exp(log_num - log_den);
                            N[k+1] = convert(Int64, rand(Binomial(S[k+1][linea,2] - S[k+1][linea-1,2], r)) + S[k+1][linea-1,2] - C[k+1]);
                        end
                        row[k+1] = linea - 1;
                    end
                end
            end

            #CONTROLLO DELLA LEAP CONDITION
            accepted = true;
            #monomeri
            if N_nuc * mins + sum(N) > maximum([(eps * monomeri)/(3+(1/(monomeri-1))+(2/(monomeri-2))), 1.0])
                accepted = false;
            end
            #crescita
            if accepted && abs(N_nuc * mins - N[mins+1]) > maximum([(eps * monomeri)/2, 1.0]) 
                accepted = false;
            end
            for k = mins+1:largest
                if accepted && abs(N[k-1]-N[k]) > maximum([(eps/2) * num[k+1], 1.0])
                    accepted = false;
                end
            end

            if accepted
                time += tau/(0.7908686 * log(V) - 1.3279955);
                contatore += 1;
                taus[i, cont] = tau;
                #nucleazione
                if row[1] != 0
                    tempmat = S[1][(row[1]+1):B[1],:];
                    S[1][1:B[1], :] = Matrix{Union{Float64, Int64}}(undef, (B[1], 2));
                    if rateNUC != 0
                        T[1] += exp(log(grw) + log(rateNUC) + log(tau));
                        C[1] += N_nuc;
                    end
                    S[1][1,1] = T[1];
                    S[1][1,2] = C[1];
                    S[1][2:(B[1] - row[1] + 1),:] = tempmat;
                    B[1] = B[1] - row[1] + 1;
                end
                #crescita
                for k = mins:largest-1
                    if row[k+1] != 0
                        tempmat = S[k+1][(row[k+1]+1):B[k+1],:];
                        S[k+1][1:B[k+1], :] = Matrix{Union{Float64, Int64}}(undef, (B[k+1], 2));
                        if ratesGR[k+1] != 0
                            T[k+1] += exp(log(grw) + log(ratesGR[k+1]) + log(tau));
                            C[k+1] += N[k+1];
                        end
                        S[k+1][1,1] = T[k+1];
                        S[k+1][1,2] = C[k+1];
                        S[k+1][2:(B[k+1] - row[k+1] + 1),:] = tempmat;
                        B[k+1] = B[k+1] - row[k+1] + 1;
                    end
                end
            else
                #nucleazione
                if row[1] + 1 > B[1]
                    if rateNUC != 0
                        S[1][B[1]+1, :] = [T[1] + exp(log(grw) + log(rateNUC) + log(tau)), C[1] + N_nuc];
                    else
                        S[1][B[1]+1, :] = [T[1], C[1] + N_nuc];
                    end
                else
                    S[1][(row[1]+2):(B[1]+2), :] = S[1][(row[1]+1):(B[1]+1),:];
                    if rateNUC != 0
                        S[1][row[1]+1, :] = [T[1] + exp(log(grw) + log(rateNUC) + log(tau)), C[1] + N_nuc];
                    else
                        S[1][row[1]+1, :] = [T[1], C[1] + N_nuc];
                    end
                end
                B[1] += 1;
                #crescita
                for k = mins:largest-1
                    if row[k+1] + 1 > B[k+1]
                        if ratesGR[k+1] != 0
                            S[k+1][B[k+1]+1, :] = [T[k+1] + exp(log(grw) + log(ratesGR[k+1]) + log(tau)), C[k+1] + N[k+1]];
                        else
                            S[k+1][B[k+1]+1, :] = [T[k+1], C[k+1] + N[k+1]];
                        end
                    else
                        S[k+1][(row[k+1]+2):(B[k+1]+2), :] = S[k+1][(row[k+1]+1):(B[k+1]+1),:];
                        if ratesGR[k+1] != 0
                            S[k+1][row[k+1]+1, :] = [T[k+1] + exp(log(grw) + log(ratesGR[k+1]) + log(tau)), C[k+1] + N[k+1]];
                        else
                            S[k+1][row[k+1]+1, :] = [T[k+1], C[k+1] + N[k+1]];
                        end
                    end
                    B[k+1] += 1;
                end
            end

            if accepted
                #CONTROLLO DELLA LEAP CONDITION
                accepted = true;
                #monomeri
                if N_nuc * mins + sum(N) > maximum([(0.75 * eps * monomeri)/(3+(1/(monomeri-1))+(2/(monomeri-2))), 1.0])
                    accepted = false;
                end
                #crescita
                if accepted && abs(N_nuc*mins - N[mins+1]) > maximum([(0.75 * eps * monomeri)/2, 1.0]) 
                    accepted = false;
                end
                for k = mins+1:largest-1
                    if accepted && abs(N[k]-N[k+1]) > maximum([(0.75 * eps/2) * num[k+1], 1.0])
                        accepted = false;
                    end
                end
                
                if N_nuc > 0
                    monomeri -= N_nuc * mins;
                    num[mins+1] += N_nuc;
                end
    
                for i = mins:largest-1
                    if N[i+1] > 0
                        monomeri -= N[i+1];
                        num[i+1] -= N[i+1];
                        num[i+2] += N[i+1];
                    end
                end

                if num[largest] > 0
                    largest += 1;
                end

                if accepted
                    tau = exp(q * log(tau));
                else
                    tau = exp(log(pstar) + log(tau));
                end
                cont += 1;

                if !taken25
                    if monomeri > 0.8 * V
                        rec25[1:threshold, i] = num[1:threshold];
                        rec25[1, i] = monomeri;
                    else
                        taken25 = true;
                    end
                end
                if !taken5 && taken25
                    if monomeri > 0.6 * V
                        rec5[1:threshold, i] = num[1:threshold];
                        rec5[1,i] = monomeri;
                    else
                        taken5 = true;
                    end
                end
                if !taken75 && taken5
                    if monomeri > 0.4 * V
                        rec75[1:threshold, i] = num[1:threshold];
                        rec75[1,i] = monomeri;
                    else
                        taken75 = true;
                    end
                end
                if !taken10 && taken75
                    if monomeri > 0.2 * V
                        rec10[1:threshold, i] = num[1:threshold];
                        rec10[1,i] = monomeri;
                    else
                        taken10 = true;
                    end
                end
            else
                tau = p * tau;
            end
        end
        comptimes[i] = tok();
        out[i,:] = Int64.(num);
        out[i, 1] = monomeri;
        println(time)
    end

    return out, comptimes, rec25, rec5, rec75, rec10;
end