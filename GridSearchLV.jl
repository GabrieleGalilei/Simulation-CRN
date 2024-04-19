using Random
using PyPlot
using Plots
using TickTock
using Distributions
using LinearAlgebra
using JLD

include("TauLeapPostCheckLV.jl")

esp::Int16 = 5;
V = 10^esp;

TLmeth = "Euler";
#TLmeth = "MidPoint";

alpha = 2.0;
beta = 0.5;
delta = 0.2;
gamma = 2.0;
alpha_scl = alpha;
beta_scl = beta/V;
delta_scl = delta/V;
gamma_scl = gamma;
cV = (alpha_scl, beta_scl, delta_scl, gamma_scl);

address = "C:\\Users\\Gabriele Galilei\\Dropbox (Politecnico Di Torino Studenti)\\PC\\Desktop\\tesi\\risultati\\";

epsvect = [0.05, 0.01, 0.1];
pvect = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
pstarvect = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]; #ricordarsi di non runnare dove p>pstar
qvect = [0.6, 0.7, 0.8, 0.9, 0.95];
sz = 2000000;
time_max = 10.5;

errormatrix = zeros(Float64, (3, 9, 7, 5));
timematrix = zeros(Float64, (3, 9, 7, 5));

sims::Int16 = 8;

for i = 1:3
    for j = 1:9
        for k in 1:7
            for l in 1:5
                if pstarvect[k] > pvect[j]
                    println("[eps: ", epsvect[i], ", p: ", pvect[j], ", pstar: ", pstarvect[k], ", q: ", qvect[l], "]")
                    _, _, _, comptimes, rec25, rec5, rec75, rec10 = TauLeapPostCheckLV(TLmeth, sz, time_max, V, cV, epsvect[i], pvect[j], pstarvect[k], qvect[l], sims);
                    timematrix[i,j,k,l] = sum(comptimes)/sims;
                    errvect = [];
                    for numero in [25, 5, 75, 10]
                        gill = load(address*"records"*"$numero"*"_Gillespie_10e"*"$esp"*".jld")["data"];
                        if numero == 25
                            rec = rec25;
                        elseif numero == 5
                            rec = rec5;
                        elseif numero == 75
                            rec = rec75;
                        else
                            rec = rec10;
                        end
                        diffsq_preys = (gill[1,1:sims] - rec[1,:]).^2;
                        diffsq_preds = (gill[2,1:sims] - rec[2,:]).^2;
                        dist = sqrt.(diffsq_preys + diffsq_preds);
                        stronger = (sum(dist)/sims)/V;
                        append!(errvect, stronger);
                    end
                    errormatrix[i,j,k,l] = maximum(errvect);
                else
                    timematrix[i,j,k,l] = NaN;
                    errormatrix[i,j,k,l] = NaN;
                end
                save(address*"errormatrix_10e"*"$esp"*".jld", "data", errormatrix)
                save(address*"timematrix_10e"*"$esp"*".jld", "data", timematrix)
            end
        end
    end
end

errormat = errormatrix./maximum(filter(!isnan,errormatrix));
timemat = timematrix./maximum(filter(!isnan,timematrix));

display(findmin(filter(!isnan, errormat + timemat)))
display(findmin(filter(!isnan, errormat)))
display(findmin(filter(!isnan, timemat)))