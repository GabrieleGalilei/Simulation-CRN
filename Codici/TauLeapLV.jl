function TauLeapLV(
    V,
    time_max::Float64,
    cV::NTuple{4, Float64},
    TLmeth::String,
    eps::Float64,
    repetition::Int16,
    sz::Int64 = 10^7
    )

    comptimes::Array{Float64} = zeros(Float64, repetition);
    out::Array{Matrix{Int64}} = [zeros(Int64, 2, sz) for _= 1:repetition];
    times::Array{Array{Float64}} = [zeros(Float64, sz) for _= 1:repetition];
    conts::Array{Int64} = zeros(Int64, repetition);
    rec25::Matrix{Int64} = zeros(Int64, 2, repetition);
    rec5::Matrix{Int64} = zeros(Int64, 2, repetition);
    rec75::Matrix{Int64} = zeros(Int64, 2, repetition);
    rec10::Matrix{Int64} = zeros(Int64, 2, repetition);

    for i = 1:repetition 
        out[i][:,1] = [V, V];
        times[i][1] = 0.0;
    end

    Threads.@threads for i = 1:repetition
        tick()
        preys::Int64 = V;
        preds::Int64 = V;
        prob = zeros(Float64, 4);
        time::Float64 = 0.0;
        cont::Int64 = 1;
        timer::Int16 = 0;
        taken25::Bool = false;
        taken5::Bool = false;
        taken75::Bool = false;
        taken10::Bool = false;

        while time < time_max
            mult::Int64 = preys*preds;
            prob::Array{Float64} = [cV[1]*preys, (cV[2]-cV[3])*mult, cV[3]*mult, cV[4]*preds];
            sumprob::Float64 = sum(prob);
            
            #tau = 0.001;
            
            #preys
            mupreys::Float64 = prob[1] - prob[2] - prob[3];
            sigmasqpreys::Float64 = prob[1] + prob[2] + prob[3];
            #predators
            mupreds::Float64 = prob[2] + prob[3] - prob[4];
            sigmasqpreds::Float64 = prob[2] + prob[3] + prob[4];

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

            time += tau;

            if cmp(TLmeth, "MidPoint") == 0
                rhopreys = preys;
                rhopreds = preds;
                addpreys = prob[1] - prob[2] - prob[3];
                addpreds = prob[3] - prob[4];
                rhopreys = preys + 0.5 * tau * addpreys;
                rhopreds = preds + 0.5 * tau * addpreds;
                prob = [cV[1]*rhopreys, (cV[2]-cV[3])*rhopreys*rhopreds, cV[3]*rhopreys*rhopreds, cV[4]*rhopreds];
            end

            N::Int64 = 0;

            if prob[1] == 0
                N = 0;
            else
                N = rand(Poisson(exp( log(prob[1]) + log(tau) )));
                preys += N;
            end

            if prob[2] == 0
                N = 0;
            else
                N = rand(Poisson(exp( log(prob[2]) + log(tau) )));
                preys -= N;
            end
            
            if prob[3] == 0
                N = 0;
            else
                N = rand(Poisson(exp( log(prob[3]) + log(tau) )));
                preys -= N;
                preds += N;
            end

            if prob[4] == 0
                N = 0;
            else
                N = rand(Poisson(exp( log(prob[4]) + log(tau) )));
                preds -= N;
            end

            if timer == 0
                cont += 1;
                out[i][:,cont] = [preys, preds];
                times[i][cont] = time;
                timer = 0;
            else
                timer -= 1;
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