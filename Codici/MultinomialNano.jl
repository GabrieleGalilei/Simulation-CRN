function MultinomialNano(
    V,
    mins::Int16,
    nuc::Float64,
    grw::Float64,
    tau::Float64,
    repetition::Int16,
    threshold::Int64 = -10
)
    sequenza = Int32.(0:(mins-1));
    monomeri::Int32 = V;

    if threshold < 0
        threshold = ceil(Int64, sqrt(sqrt(monomeri)));
        #threshold = ceil(Int64, 1.2 * sqrt(monomeri));
        #threshold = V+1;
    end

    comptimes::Array{Float64} = zeros(Float64, repetition);

    out::Matrix{Int64} = zeros(Int64, repetition, threshold);

    rec25::Matrix{Int64} = zeros(Int64, threshold, repetition);
    rec5::Matrix{Int64} = zeros(Int64, threshold, repetition);
    rec75::Matrix{Int64} = zeros(Int64, threshold, repetition);
    rec10::Matrix{Int64} = zeros(Int64, threshold, repetition);

    Threads.@threads for i = 1:repetition
        tick()
        num = zeros(Int64, threshold);
        monomeri = V;
        rateNUC::Float64 = 0;
        largest::Int64 = mins+1;
        newpar::Float64 = nuc/grw;
        taken25::Bool = false;
        taken5::Bool = false;
        taken75::Bool = false;
        taken10::Bool = false;

        time = 0.0;

        while monomeri > 0
            if monomeri >= mins
                rateNUC = newpar * exp(sum(log.(monomeri .- sequenza)));
            else 
                rateNUC = 0.0;
            end
            
            ratesGR = zeros(Int64, largest);
            for k = mins:(largest-1)
                if num[k+1] > 0 
                    ratesGR[k+1] = monomeri * num[k+1];
                end
            end

            if rateNUC > 0
                N_nuc = rand( Binomial(floor(monomeri/mins), 1 - exp( (-1) * grw * rateNUC * tau) ));
                #N_nuc = rand( Binomial(floor(monomeri/mins), ((grw * rateNUC * tau)^(grw * rateNUC * tau))*(exp( (-1) * grw * rateNUC * tau)/factorial(floor(Int64, grw * rateNUC * tau) ))));
            else
                N_nuc = 0;
            end

            prob = zeros(Float64, largest);
            for k = mins:(largest-1)
                prob[k+1] = 1 - exp((-1) * grw * ratesGR[k+1] * tau);
                #prob[k+1] = ((grw * ratesGR[k+1] * tau)^(grw * ratesGR[k+1] * tau))*(exp((-1)*grw * ratesGR[k+1] * tau)/ factorial(floor(Int64, grw * ratesGR[k+1] * tau)))
            end
            prob_null = 1 - sum(prob);

            #=
            found = false;
            while !found
                Crescite = rand(Multinomial(convert(Int64, monomeri), [prob; prob_null]));
                accept = true;
                k = mins;
                while accept
                    if Crescite[k+1] > num[k+1]
                        accept = false;
                    end
                    k += 1;
                end
                if accept
                    found = true;
                end
            end
            =#

            Crescite = rand(Multinomial(convert(Int64, monomeri - N_nuc*mins), [prob; prob_null]));
            for k = mins:largest-1
                Crescite[k+1] = min(Crescite[k+1], num[k+1]);
            end

            if N_nuc > 0
                monomeri -= N_nuc * mins;
                num[mins+1] += N_nuc;
            end

            for k = mins:largest-1
                if Crescite[k+1] > 0
                    monomeri -= convert(Int64, Crescite[k+1]);
                    num[k+1] -= convert(Int64, Crescite[k+1]);
                    num[k+2] += convert(Int64, Crescite[k+1]);
                end
            end

            if Crescite[largest] > 0
                largest += 1; # largest is the dimension of the biggest particle + 1 
            end

            time += tau;

            if !taken25
                rec25[1,i] = monomeri;
                if monomeri > 0.8 * V
                    rec25[1:soffitto, i] = num[1:soffitto];
                    rec25[1, i] = monomeri;
                else
                    taken25 = true;
                end
            end
            if !taken5 && taken25
                if monomeri > 0.6 * V
                    rec5[1:soffitto, i] = num[1:soffitto];
                    rec5[1,i] = monomeri;
                else
                    taken5 = true;
                end
            end
            if !taken75 && taken5
                if monomeri > 0.4 * V
                    rec75[1:soffitto, i] = num[1:soffitto];
                    rec75[1,i] = monomeri;
                else
                    taken75 = true;
                end
            end
            if !taken10 && taken75
                if monomeri > 0.2 * V
                    rec10[1:soffitto, i] = num[1:soffitto];
                    rec10[1,i] = monomeri;
                else
                    taken10 = true;
                end
            end
            
        end
        comptimes[i] = tok();
        out[i, :] = Int64.(num);
        out[i, 1] = monomeri;
    end

    return out, comptimes, rec25, rec5, rec75, rec10;
end





