using Random
using PyPlot
using Plots
using TickTock
using Distributions
using LinearAlgebra
using JLD

include("TauLeapPostCheckNano2.jl")

esp::Int16 = 5;
V = 10^esp;

#scaling = "Classical";
scaling = "Alternative";

grw::Float64 = 1.0;
nuc::Float64 = 5.0;
mins::Int16 = 3;

if cmp(scaling, "Classical") == 0
    nuc = exp(log(nuc) + (1-mins) * log(V));
    grw = exp(log(grw) - log(V));
    slope::Float64 = 2.0284;
    intercept::Float64 = -0.4442;
    soffitto = ceil(Int64, 2 * sqrt(sqrt(V)));
elseif cmp(scaling, "Alternative") == 0
    nuc = exp(log(nuc) + (1-mins) * log(V));
    slope::Float64 = 2.601746693;
    intercept::Float64 = -0.656053137;
    soffitto = ceil(Int64, 1.2 * sqrt(V));
end

address = "C:\\Users\\Gabriele Galilei\\Dropbox (Politecnico Di Torino Studenti)\\PC\\Desktop\\tesi\\risultati\\";

epsvect = [0.6, 0.7, 0.8, 0.9];
pvect = [0.3, 0.4, 0.5, 0.6];
pstarvect = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]; #ricordarsi di non runnare dove p>pstar
qvect = [0.7, 0.8, 0.9];
sz = 100000;

errormatrix = zeros(Float64, (4, 4, 6, 3));
timematrix = zeros(Float64, (4, 4, 6, 3));

sims::Int16 = 2;

for i = 1:4
    for j = 1:4
        for k in 1:6
            for l in 1:3
                if pstarvect[k] > pvect[j]
                    println("[eps: ", epsvect[i], ", p: ", pvect[j], ", pstar: ", pstarvect[k], ", q: ", qvect[l], "]")
                    _, comptimes, rec25, rec5, rec75, rec10 = TauLeapPostCheckNano2(sz, epsvect[i], pvect[j], pstarvect[k], qvect[l], V, mins, nuc, grw, sims, slope, intercept);
                    timematrix[i,j,k,l] = sum(comptimes)/sims;
                    errvect = [];
                    for numero in [25, 5, 75, 10]
                        gill = load(address*"Nanoparticles_records"*"$numero"*"_"*scaling*"_Gillespie_10e"*"$esp"*".jld")["data"];
                        if numero == 25
                            rec = rec25;
                        elseif numero == 5
                            rec = rec5;
                        elseif numero == 75
                            rec = rec75;
                        else
                            rec = rec10;
                        end
                        diff = zeros(Float64, sims);
                        for x = 1:sims
                            for y = 1:soffitto
                                diff[x] += (gill[y,x] - rec[y,x])^2;
                            end
                            diff[x] = sqrt(diff[x]);
                        end
                        stronger = (sum(diff)/sims)/V;
                        append!(errvect, stronger);
                    end
                    errormatrix[i,j,k,l] = maximum(errvect);
                else
                    timematrix[i,j,k,l] = NaN;
                    errormatrix[i,j,k,l] = NaN;
                end
                save(address*"Nano_"*scaling*"_errormatrix_10e"*"$esp"*".jld", "data", errormatrix)
                save(address*"Nano_"*scaling*"_timematrix_10e"*"$esp"*".jld", "data", timematrix)
            end
        end
    end
end

errormat = errormatrix./maximum(filter(!isnan,errormatrix));
timemat = timematrix./maximum(filter(!isnan,timematrix));

display(findmin(filter(!isnan, errormat + timemat)))
display(findmin(filter(!isnan, errormatrix)))
display(findmin(filter(!isnan, timematrix)))