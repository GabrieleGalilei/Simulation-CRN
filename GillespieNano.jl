function GillespieNano(
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
    monomeri::Int32 = V;

    if threshold < 0
        #threshold = ceil(Int64, 1.2 * sqrt(monomeri));
        threshold = V+1;
    end
    comptimes::Array{Float64} = zeros(Float64, repetition);

    if slope > 2.6 #Alternative
        time25 = 0.2 * 4.49;
        time5 = 0.4 * 4.49;
        time75 = 0.6 * 4.49;
        time10 = 0.8 * 4.49;
        soffitto = threshold;
    else #Classical
        time25 = 0.2 * 4.49;
        time5 = 0.4 * 4.49;
        time75 = 0.6 * 4.49;
        time10 = 0.8 * 4.49;
        threshold = ceil(Int64, 2 * sqrt(sqrt(monomeri)));
        #threshold = V+1;
        soffitto = threshold;
    end

    out::Matrix{Int64} = zeros(Int64, repetition, threshold);

    rec25::Matrix{Int64} = zeros(Int64, soffitto, repetition);
    rec5::Matrix{Int64} = zeros(Int64, soffitto, repetition);
    rec75::Matrix{Int64} = zeros(Int64, soffitto, repetition);
    rec10::Matrix{Int64} = zeros(Int64, soffitto, repetition);

    Threads.@threads for i = 1:repetition
        tick()
        num = zeros(Int64, threshold);
        monomeri = V;
        rateNUC::Float64 = 0;
        ratesGR = zeros(Int64, threshold);
        prob = zeros(Float64, threshold);
        largest::Int64 = mins+1;
        taken25::Bool = false;
        taken5::Bool = false;
        taken75::Bool = false;
        taken10::Bool = false;

        newpar::Float64 = nuc/grw;

        time = 0.0;

        while monomeri > 0
            #rate[1] is the rate of the nucleation
            #rate[i+1] is the rate of  growth of nanoparticles of size i
            #and since 
            #num[i+1] is the number of nanoparticles of size i
            #then
            #rate[i+1] is the rate of  growth of nanoparticles whose amount is in num[i+1]
            if monomeri >= mins
                rateNUC = newpar * exp(sum(log.(monomeri .- sequenza)));
            else 
                rateNUC = 0.0;
            end #NB lavoro con tutti i rate divisi per growth, anche la nucleazione

            for k = mins:(largest-1)
                ratesGR[k+1] = monomeri * num[k+1]; #* (i+1)^pow;
            end

            sumratesGR::Int64 = sum(ratesGR);
            sumrates = sumratesGR + rateNUC;
                
            prob[1] = exp(log(rateNUC) - log(sumrates)); 

            for k = (mins+1):largest
                
                if ratesGR[k]>0;
                    prob[k] = exp(log(ratesGR[k])-log(sumrates)); 
                else
                    prob[k] = 0.0
                end
            end
            fires::Int64 = rand(Categorical(prob));

            if fires == threshold
                error("Dimension of nanoparticles is exceding threshold")
            end

            if fires == 1 # && monomeri>=mins
                monomeri = monomeri - mins;
                num[mins+1] = num[mins+1] + 1;
            else
                if fires == largest
                    largest = largest+1; # largest is the dimension of the biggest particle + 1 
                end
                monomeri = monomeri - 1;
                num[fires] = num[fires] - 1;
                num[fires+1] = num[fires+1] + 1;
            end

            time += rand(Exponential(grw * sumrates * 10 / (V^2)));
            #time += rand(Exponential(grw * sumrates));

            if !taken25
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
        out[i,:] = Int64.(num);
        out[i, 1] = monomeri;
        println(time)
    end

    return out, comptimes, rec25, rec5, rec75, rec10;

end