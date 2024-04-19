function GillespieLV(
    V,
    time_max::Float64,
    cV::NTuple{4, Float64},
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
    
    Threads.@threads for i = 1:repetition
        tick()
        preys::Int64 = V;
        preds::Int64 = V;
        prob = zeros(Float64, 4);
        time::Float64 = 0.0;
        cont::Int64 = 1;
        timer::Int16 = 1000;
        taken25::Bool = false;
        taken5::Bool = false;
        taken75::Bool = false;
        taken10::Bool = false;
        
        for i = 1:repetition 
            out[i][:,1] = [preys, preds];
            times[i][1] = 0.0;
        end

        while time < time_max
            if preys + preds > 0
                mult::Int64 = preys*preds;
                prob::Array{Float64} = [cV[1]*preys, (cV[2]-cV[3])*mult, cV[3]*mult, cV[4]*preds];
                sumprob::Float64 = mult*cV[2] + cV[1]*preys + cV[4]*preds;
                r1 = rand();
                tau::Float64 = (1/sumprob) * log(1.0 / r1);
                prob = prob./sumprob;
                fires::Int64 = rand(Categorical(prob));

                if fires == 1
                    preys += 1;
                elseif fires == 2
                    preys -=1;
                elseif fires == 3
                    preys -= 1;
                    preds += 1;
                else
                    preds -= 1;
                end
            else
                tau = 1.0;
            end

            time += tau;
            if timer == 0
                cont += 1;
                out[i][:,cont] = [preys, preds];
                times[i][cont] = time;
                timer = 100;
            else
                timer -= 1;
            end

            if time > 2.5 - 0.001 && time < 2.5 + 0.001 && !taken25
                taken25 = true;
                rec25[:,i] = [preys, preds];
            end
            if time > 5.0 - 0.001 && time < 5.0 + 0.001 && !taken5
                taken5 = true;
                rec5[:,i] = [preys, preds];
            end
            if time > 7.5 - 0.001 && time < 7.5 + 0.001 && !taken75
                taken75 = true;
                rec75[:,i] = [preys, preds];
            end
            if time > 10.0 - 0.001 && time < 10.0 + 0.001 && !taken10
                taken10 = true;
                rec10[:,i] = [preys, preds];
            end
        end
        comptimes[i] = tok();
        conts[i] = cont;
    end

    return out, times, conts, comptimes, rec25, rec5, rec75, rec10;

end