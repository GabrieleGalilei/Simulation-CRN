function TauLeapNano(
    V,
    mins::Int16,
    nuc::Float64,
    grw::Float64,
    eps::Float64,
    repetition::Int16,
    threshold::Int64 = -10
)

    sequenza = Int32.(0:(mins-1));
    monomeri::Int32 = V;

    if threshold < 0
        threshold = ceil(Int64, 3 * sqrt(monomeri));
    end
    comptimes::Array{Float64} = zeros(Float64, repetition);

    out::Matrix{Int64} = zeros(Int64, repetition, threshold);

    Threads.@threads for i = 1:repetition
        tick()
        num = zeros(Int64, threshold);
        monomeri = V;
        rateNUC::Float64 = 0;
        ratesGR = zeros(Int64, threshold);
        largest::Int64 = mins + 1; #indica la massima reazione che può avvenire

        newpar::Float64 = nuc/grw;

        while monomeri > 0

            #rate[1] is the rate of the nucleation
            #rate[i+1] is the rate of  growth of nanoparticles of size i
            #and since 
            #num[i+1] is the number of nanoparticles of size i
            #then
            #rate[i+1] is the rate of  growth of nanoparticles whose amount is  in num[i+1]

            if monomeri >= mins
                rateNUC = newpar * exp(sum(log.(monomeri .- sequenza)));
            end
            for k = mins:(largest-1)
                ratesGR[k+1] = monomeri * num[k+1];
            end

            
            # CALCOLO DI TAU
            mu = zeros(Float64, largest+1);
            sigmasq = zeros(Float64, largest+1);

            #nucleazione
            if rateNUC != 0
                mu[1] = ((-1)*mins) * rateNUC;
                sigmasq[1] = (mins^2) * rateNUC;
                mu[mins+1] = rateNUC;
                sigmasq[mins+1] = rateNUC;
            end
            #crescita
            for k = mins:largest-1
                if ratesGR[i+1] != 0
                    mu[1] -= ratesGR[k+1];
                    sigmasq[1] += ratesGR[k+1];
                    mu[k+1] -= ratesGR[k+1];
                    sigmasq[k+1] += ratesGR[k+1];
                    mu[k+2] += ratesGR[k+1];
                    sigmasq[k+2] += ratesGR[k+2];
                end
            end

            mu = grw.*mu;
            sigmasq = grw.*sigmasq;

            tau = Inf;

            #monomeri
            temp = log(eps) + log(monomeri) - log(3 + (1/(monomeri-1)) + (2/(monomeri-2)));
            temp = maximum([temp, 0.0]);
            firstterm = temp - log(abs(mu[1]));
            secondterm = 2.0 * temp - log(sigmasq[1]);
            tau = minimum([firstterm, secondterm, tau]);
            #particelle
            for k = mins:(largest-1)
                if num[k+1] == 0
                    firstterm = - log(abs(mu[k+1]));
                    secondterm = - log(abs(sigmasq[k+1]));
                    tau = minimum([firstterm, secondterm, tau]);
                else
                    temp = log(eps) + log(num[k+1]) - log(2);
                    temp = maximum([temp, 0.0]);
                    firstterm = temp - log(abs(mu[k+1]));
                    secondterm = 2.0 * temp - log(abs(sigmasq[k+1]));
                    tau = minimum([firstterm, secondterm, tau]);
                end
            end
            tau = exp(tau);
            
            #tau::Float64 = 0.00001;
            #eps = tau;
    
            N_grw::Array{Int64} = zeros(Int64, threshold);
            N_nuc::Int64 = 0;

            mult = log(grw) + log(tau);

            if monomeri >= mins
                N_nuc = rand(Poisson(exp(log(rateNUC) + mult)));
                println("nucleazione: ", N_nuc)
            end #NB lavoro con tutti i rate divisi per growth, anche la nucleazione

            for k = mins:(largest-1)
                if ratesGR[k+1] != 0
                    N_grw[k+1] = rand(Poisson(exp(mult + log(ratesGR[k+1]))));
                    println("crescita di P", i, " : ", N_grw[k+1])
                end
            end

            if N_grw[largest] != 0
                largest += 1;
            end

            if largest == threshold
                error("Dimension of nanoparticles is exceding threshold")
            end

            if N_nuc != 0
                monomeri -= N_nuc * mins;
                num[mins+1] += N_nuc;
            end

            for k = mins:(largest-1)
                if N_grw[k+1] != 0
                    monomeri -= N_grw[k+1];
                    num[k+1] -= N_grw[k+1];
                    num[k+2] += N_grw[k+1];
                end
            end
            
            if any(x -> x<0, num) || monomeri < 0
                error("Something went negative")
            end
        end
        comptimes[i] = tok();
        out[i,:] = num;
        out[i, 1] = monomeri;
    end

    return out, comptimes;
end