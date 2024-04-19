using Random
using PyPlot
using Plots
using Distributions
using LinearAlgebra
using TickTock
using JLD

include("GillespieNano.jl")
include("TauLeapNano.jl")
include("TauLeapPostCheckNano.jl")
include("TauLeapPostCheckNano2.jl")
include("MultinomialNano.jl")

esp::Int16 = 6;
V = 10^esp;

simulations = 8;

#algorithm = "Gillespie";
#algorithm = "TauLeap";
#algorithm = "TauLeapPostLeapCheck";
algorithm = "Multinomial";

scaling = "Classical";
#scaling = "Alternative";

repetition::Int16 = 8;

if cmp(scaling, "Classical") == 0
    soffitto = ceil(Int64, sqrt(sqrt(V)));
    #soffitto = V+1;
elseif cmp(scaling, "Alternative") == 0
    soffitto = ceil(Int64, 1.2 * sqrt(V));
    #soffitto = V+1;
end

sz = 100000;

records25::Matrix{Int64} = zeros(Int64, soffitto, simulations);
records5::Matrix{Int64} = zeros(Int64, soffitto, simulations);
records75::Matrix{Int64} = zeros(Int64, soffitto, simulations);
records10::Matrix{Int64} = zeros(Int64, soffitto, simulations);
outres::Matrix{Int64} = zeros(Int64, soffitto, simulations);

grw::Float64 = 1.0;
nuc::Float64 = 5.0;
mins::Int16 = 3;

if cmp(algorithm, "Multinomial") == 0
    taumon::Float64 = 10^(-7);
end

if cmp(scaling, "Classical") == 0
    nuc = exp(log(nuc) + (1 - mins) * log(V));
    grw = exp(log(grw) - log(V));
    slope::Float64 = 2.0284;
    intercept::Float64 = -0.4442;
elseif cmp(scaling, "Alternative") == 0
    nuc = exp(log(nuc) + (1 - mins) * log(V));
    slope::Float64 = 2.601746693;
    intercept::Float64 = -0.656053137;
end

if cmp(scaling, "Classical") == 0
    eps::Float64 = 0.6;
    p::Float64 = 0.3;
    pstar::Float64 = 0.4;
    q::Float64 = 0.7;
elseif cmp(scaling, "Alternative") == 0
    eps::Float64 = 0.9;
    p::Float64 = 0.3;
    pstar::Float64 = 0.6;
    q::Float64 = 0.9;
end

global comptimes = zeros(Float64, repetition);

address = "C:\\Users\\Gabriele Galilei\\Dropbox (Politecnico Di Torino Studenti)\\PC\\Desktop\\tesi\\risultati\\";

global rep = 0;

while rep < simulations
    if cmp(algorithm, "Gillespie") == 0
        global out, comptimes_add, rec25, rec5, rec75, rec10 = GillespieNano(V, mins, nuc, grw, repetition, slope, intercept);
        global comptimes += comptimes_add;
        records25[:, rep+1:rep+repetition] = rec25;
        records5[:, rep+1:rep+repetition] = rec5;
        records75[:, rep+1:rep+repetition] = rec75;
        records10[:, rep+1:rep+repetition] = rec10;
        outres[:, rep+1:rep+repetition] = out';
        println(((rep+repetition)/simulations)*100," %")
    elseif cmp(algorithm, "TauLeap") == 0
        global out, comptimes_add = TauLeapNano(V, mins, nuc, grw, eps, repetition);
        global comptimes += comptimes_add;
        records25[:, rep+1:rep+repetition] = rec25;
        records5[:, rep+1:rep+repetition] = rec5;
        records75[:, rep+1:rep+repetition] = rec75;
        records10[:, rep+1:rep+repetition] = rec10;
        #outres[:, rep+1:rep+repetition] = out;
        println(((rep+repetition)/simulations) * 100," %")
    elseif cmp(algorithm, "TauLeapPostLeapCheck") == 0
        global out, comptimes_add, rec25, rec5, rec75, rec10 = TauLeapPostCheckNano2(sz, eps, p, pstar, q, V, mins, nuc, grw, repetition, slope, intercept);
        global comptimes += comptimes_add;
        records25[:, rep+1:rep+repetition] = rec25;
        records5[:, rep+1:rep+repetition] = rec5;
        records75[:, rep+1:rep+repetition] = rec75;
        records10[:, rep+1:rep+repetition] = rec10;
        #outres[:, rep+1:rep+repetition] = out';
        println(((rep+repetition)/simulations) * 100," %")
    elseif cmp(algorithm, "Multinomial") == 0
        global out, comptimes_add, rec25, rec5, rec75, rec10 = MultinomialNano(V, mins, nuc, grw, taumon, repetition);
        global comptimes += comptimes_add;
        records25[:, rep+1:rep+repetition] = rec25;
        records5[:, rep+1:rep+repetition] = rec5;
        records75[:, rep+1:rep+repetition] = rec75;
        records10[:, rep+1:rep+repetition] = rec10;
        println(((rep+repetition)/simulations) * 100," %")
    end
    global rep += repetition;
end

plot_array = []; 
for i = 1:repetition
    p1 = Plots.plot();
    Plots.bar!(p1, out[i,:], legend=false, xtickfontsize=5, ytickfontsize=5, titlefontsize= 6, xlabelfontsize=6, ylabelfontsize=6)
    ylims!(0, maximum(out[1,:]))
    #xlims!(1,200)
    xlims!(0.5, maximum(findall(x -> x > 0, out[i,:]))+1)
    xlabel!("size of the nanoparticle")
    ylabel!("no. of nanoparticles")
    push!(plot_array, p1)
end
for i = repetition+1:8
    p1 = Plots.plot();
    xlabel!("size of the nanoparticle")
    ylabel!("no. of nanoparticles")
    push!(plot_array, p1)
end
prep = Plots.plot(plot_array[1], plot_array[2], plot_array[3], plot_array[4], plot_array[5], plot_array[6], plot_array[7], plot_array[8], layout=(2,4))
display(prep)
if cmp(algorithm, "Gillespie") == 0
    #Plots.savefig(prep, address*"Nanoparticles_"*scaling*"_Gillespie_10e"*"$esp"*".png");
    save(address*"Nanoparticles_records25_"*scaling*"_Gillespie_10e"*"$esp"*".jld", "data", records25)
    save(address*"Nanoparticles_records5_"*scaling*"_Gillespie_10e"*"$esp"*".jld", "data", records5)
    save(address*"Nanoparticles_records75_"*scaling*"_Gillespie_10e"*"$esp"*".jld", "data", records75)
    save(address*"Nanoparticles_records10_"*scaling*"_Gillespie_10e"*"$esp"*".jld", "data", records10)
    #save(address*"Nanoparticles_out_"*scaling*"_Gillespie_10e"*"$esp"*".jld", "data", outres)
elseif cmp(algorithm, "TauLeap") == 0 || cmp(algorithm, "TauLeapPostLeapCheck") == 0
    Plots.savefig(prep, address*"Nanoparticles_"*algorithm*"_"*scaling*"_10e"*"$esp"*"_eps="*"$eps"*".png");
    #save(address*"Nanoparticles_records25_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records25)
    #save(address*"Nanoparticles_records5_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records5)
    #save(address*"Nanoparticles_records75_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records75)
    #save(address*"Nanoparticles_records10_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records10)
    #save(address*"Nanoparticles_out_"*scaling*"_"*algorithm*"_10e"*"$esp"*".jld", "data", outres)
elseif cmp(algorithm, "Multinomial") == 0
    Plots.savefig(prep, address*"Nanoparticles_"*algorithm*"_"*scaling*"_10e"*"$esp"*"_tau="*"$taumon"*".png");
    save(address*"Nanoparticles_records25_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_tau="*"$taumon"*".jld", "data", records25)
    save(address*"Nanoparticles_records5_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_tau="*"$taumon"*".jld", "data", records5)
    save(address*"Nanoparticles_records75_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_tau="*"$taumon"*".jld", "data", records75)
    save(address*"Nanoparticles_records10_"*scaling*"_"*algorithm*"_10e"*"$esp"*"_tau="*"$taumon"*".jld", "data", records10)
end
display(sum(comptimes) / rep)

#p3 = Plots.plot(dpi=1000);
#for i = 1:5
#    Plots.plot!(p3, cumsum(taus[i,:]), taus[i,:], linetype=:steppost, lw=2, legend=false);
#end
#xlabel!("time")
#ylabel!("accepted tau (s)")
#display(p3)