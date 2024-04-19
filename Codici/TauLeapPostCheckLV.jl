function TauLeapPostCheckLV(
    TLmeth::String,
    sz::Int64,
    time_max::Float64,
    V,
    cV,
    eps::Float64,
    p::Float64,
    pstar::Float64,
    q::Float64,
    repetition::Int16)

    out::Array{Matrix{Int64}} = [zeros(Int64, 2, sz) for _= 1:repetition];
    comptimes::Array{Float64} = zeros(Float64, repetition);
    times::Array{Array{Float64}} = [zeros(Float64, sz) for _= 1:repetition];
    conts::Array{Int64} = zeros(Int64, repetition);
    rec25::Matrix{Int64} = zeros(Int64, 2, repetition);
    rec5::Matrix{Int64} = zeros(Int64, 2, repetition);
    rec75::Matrix{Int64} = zeros(Int64, 2, repetition);
    rec10::Matrix{Int64} = zeros(Int64, 2, repetition);

    Threads.@threads for i = 1:repetition
        tick()
        preys::Int64 = V;
        preds::Int64 = V;
        taken25::Bool = false;
        taken5::Bool = false;
        taken75::Bool = false;
        taken10::Bool = false;

        times[i][1] = 0.0;
        time = 0.0;

        out[i][:,1] = [preys, preds];

        S = [Matrix{Union{Float64, Int64}}(undef, (sz, 2)) for _ = 1:4];
        B = ones(Int64, 4);
        T = zeros(Float64, 4);
        C = zeros(Int64, 4);

        global cont = 1;
        for k = 1:4
            S[k][1,:] = [0.0, 0];
        end

        mult::Int64 = preys*preds;
        props::Array{Float64} = [cV[1]*preys, (cV[2]-cV[3])*mult, cV[3]*mult, cV[4]*preds];
        #preys
        mupreys::Float64 = props[1] - props[2] - props[3];
        sigmasqpreys::Float64 = props[1] + props[2] + props[3];
        #predators
        mupreds::Float64 = props[2] + props[3] - props[4];
        sigmasqpreds::Float64 = props[2] + props[3] + props[4];

        tau = Inf;

        temp = log(eps) + log(preys) - log(2);
        temp = maximum([temp, 0.0]);
        firstterm = temp - log(abs(mupreys));
        secondterm = 2.0*temp - log(sigmasqpreys);
        tau = minimum([firstterm, secondterm, tau]);
        temp = log(eps)+ log(preds) - log(2);
        temp = maximum([temp, 0.0]);
        firstterm = temp - log(abs(mupreds));
        secondterm = 2.0*temp - log(sigmasqpreds);
        tau = minimum([firstterm, secondterm, tau]);
        tau = exp(tau);

        while time < time_max
            N = Array{Int64}(undef, 4);
            row = Array{Int64}(undef, 4);

            props = [cV[1]*preys, (cV[2]-cV[3])*preys*preds, cV[3]*preys*preds, cV[4]*preds];

            if cmp(TLmeth, "MidPoint") == 0
                rhopreys = preys;
                rhopreds = preds;
                addpreys = props[1] - props[2] - props[3];
                addpreds = props[3] - props[4];
                rhopreys = preys + 0.5 * tau * addpreys;
                rhopreds = preds + 0.5 * tau * addpreds;
                props = [cV[1]*rhopreys, (cV[2]-cV[3])*rhopreys*rhopreds, cV[3]*rhopreys*rhopreds, cV[4]*rhopreds];
            end

            for k = 1:4
                # definiamo row e N_k
                if props[k] == 0.0
                    N[k] = 0;
                    row[k] = 0;
                else
                    if props[k] * tau + T[k] >= S[k][B[k], 1]
                        N[k] = rand(Poisson( props[k] * tau + T[k] - S[k][B[k],1])) + S[k][B[k],2] - C[k];
                        row[k] = B[k];
                    else
                        linea = 2;
                        while props[k] * tau + T[k] >= S[k][linea,1] && linea < B[k]
                            linea += 1;
                        end
                        if S[k][linea,1] - S[k][linea-1,1] == 0.0
                            N[k] = 0;
                        else
                            log_num = log(props[k] * tau + T[k] - S[k][linea-1, 1]);
                            log_den = log(S[k][linea,1] - S[k][linea-1,1]);
                            log_r = log_num - log_den;
                            r = exp(log_r);
                            N[k] = convert(Int64, rand(Binomial(S[k][linea,2] - S[k][linea-1,2], r))) + S[k][linea-1,2] - C[k];
                        end
                        row[k] = linea - 1;
                    end
                end
            end

            #CONTROLLO DELLA LEAP CONDITION
            accepted = true;
            #prede
            if abs(N[1]-N[2]-N[3]) > maximum([(eps/2) * preys, 1.0])
                accepted = false;
            end
            #predatori
            if abs(N[3]-N[4]) > maximum([(eps/2) * preds, 1.0]) && accepted
                accepted = false;
            end

            if accepted
                for k = 1:4
                    temp = S[k][(row[k]+1):B[k],:];
                    S[k][1:B[k], :] = Matrix{Union{Float64, Int64}}(undef, (B[k], 2));
                    log_term = log(props[k]) + log(tau);
                    T[k] += exp(log_term);
                    C[k] += N[k];
                    S[k][1,1] = T[k];
                    S[k][1,2] = C[k];
                    S[k][2:(B[k] - row[k] + 1),:] = temp;
                    B[k] = B[k] - row[k] + 1;
                end
            else
                for k = 1:4
                    log_term = log(props[k]) + log(tau);
                    if row[k] + 1 > B[k]
                        S[k][B[k]+1, :] = [T[k] + exp(log_term), C[k] + N[k]];
                    else
                        S[k][(row[k]+2):(B[k]+2), :] = S[k][(row[k]+1):(B[k]+1),:];
                        S[k][row[k]+1, :] = [T[k] + exp(log_term), C[k] + N[k]];
                    end
                    B[k] += 1;
                end
            end

            if accepted
                time += tau;
                #CONTROLLO DELLA LEAP CONDITION
                accepted = true;
                #prede
                if abs(N[1]-N[2]-N[3]) > maximum([(0.75*eps/2) * preys, 1.0])
                    accepted = false;
                end
                #predatori
                if abs(N[3]-N[4]) > maximum([(0.75*eps/2) * preds, 1.0]) && accepted
                    accepted = false;
                end

                if accepted
                    tau = exp(q*log(tau));
                else
                    tau = exp(log(pstar) + log(tau));
                end
                cont += 1;
                times[i][cont] = time;
                preys = preys + N[1] - N[2] - N[3];
                preds = preds + N[3] - N[4];
                out[i][:, cont] = [preys, preds];
            else
                tau = exp(log(p) + log(tau));
            end

            if time > 2.5 - eps && time < 2.5 + eps && !taken25
                taken25 = true;
                rec25[:,i] = [preys, preds];
            end
            if time > 5.0 - eps && time < 5.0 + eps && !taken5
                taken5 = true;
                rec5[:,i] = [preys, preds];
            end
            if time > 7.5 - eps && time < 7.5 + eps && !taken75
                taken75 = true;
                rec75[:,i] = [preys, preds];
            end
            if time > 10.0 - eps && time < 10.0 + eps && !taken10
                taken10 = true;
                rec10[:,i] = [preys, preds];
            end
        end
        comptimes[i] = tok();
        conts[i] = cont;
    end

    return out, times, conts, comptimes, rec25, rec5, rec75, rec10;
end