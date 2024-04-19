#using Random
#using PyPlot
using Plots
#using DataFrames
#using CSV
using TickTock
#using DelimitedFiles
using Distributions
using LinearAlgebra
using DifferentialEquations
using JLD

include("LotkaVolterraDet.jl")
include("GillespieLV.jl")
include("TauLeapLV.jl")
include("TauLeapPostCheckLV.jl")

#Plots.closeall()

esp::Int16 = 2;
V = 10^esp;

simulations::Int64 = 1;

records25::Matrix{Int64} = zeros(Int64, 2, simulations);
records5::Matrix{Int64} = zeros(Int64, 2, simulations);
records75::Matrix{Int64} = zeros(Int64, 2, simulations);
records10::Matrix{Int64} = zeros(Int64, 2, simulations);

algorithm = "Gillespie";
#algorithm = "TauLeap";
#algorithm = "TauLeapPostLeapCheck";

repetition::Int16 = 1;

TLmeth = "Euler";
#TLmeth = "MidPoint";

sz = 5000000;

# Costanti per Lotka Volterra
alpha::Float64 = 2.0;
beta::Float64 = 0.5;
delta::Float64 = 0.2;
gamma::Float64 = 2.0;
c = (alpha, beta, delta, gamma);

alpha_scl = alpha;
beta_scl = beta/V;
delta_scl = delta/V;
gamma_scl = gamma;
cV = (alpha_scl, beta_scl, delta_scl, gamma_scl);

# Costanti per Post Leap Check
eps::Float64 = 0.01;
p::Float64 = 0.1;
pstar::Float64 = 0.5;
q::Float64 = 0.7;

# Limiti temporali
time_max::Float64 = 20;

address = "C:\\Users\\Gabriele Galilei\\Dropbox (Politecnico Di Torino Studenti)\\PC\\Desktop\\tesi\\risultati\\";

global comptimes = zeros(Float64, repetition);

global rep = 0;
while rep < simulations
    if cmp(algorithm, "Gillespie") == 0
        global out, times, conts, comptimes_add, rec25, rec5, rec75, rec10 = GillespieLV(V, time_max, cV, repetition, sz);
        #records25[:, rep+1:rep+repetition] = rec25;
        #records5[:, rep+1:rep+repetition] = rec5;
        #records75[:, rep+1:rep+repetition] = rec75;
        #records10[:, rep+1:rep+repetition] = rec10;
        println(((rep+repetition)/simulations)*100," %")
        global comptimes += comptimes_add;
    elseif cmp(algorithm, "TauLeap") == 0
        global out, times, conts, comptimes_add, rec25, rec5, rec75, rec10 = TauLeapLV(V, time_max, cV, TLmeth, eps, repetition, sz);
        global comptimes += comptimes_add;
        records25[:, rep+1:rep+repetition] = rec25;
        records5[:, rep+1:rep+repetition] = rec5;
        records75[:, rep+1:rep+repetition] = rec75;
        records10[:, rep+1:rep+repetition] = rec10;
        println(((rep+repetition)/simulations)*100," %")
    elseif cmp(algorithm, "TauLeapPostLeapCheck") == 0
        global out, times, conts, comptimes_add, rec25, rec5, rec75, rec10 = TauLeapPostCheckLV(TLmeth, sz, time_max, V, cV, eps, p, pstar, q, repetition);
        global comptimes += comptimes_add;
        records25[:, rep+1:rep+repetition] = rec25;
        records5[:, rep+1:rep+repetition] = rec5;
        records75[:, rep+1:rep+repetition] = rec75;
        records10[:, rep+1:rep+repetition] = rec10;
        println(((rep+repetition)/simulations)*100," %")
    end
    global rep += repetition;
end

p1 = Plots.plot(dpi=1000);
p2 = Plots.plot(dpi=1000);
if cmp(algorithm, "Gillespie") == 0
    sol = LotkaVolterraDet(c, [1,1], 0.0, time_max);
    rng_x = range(0.0, time_max, 10000);
    rng_y = sol(rng_x);
    #=
    y = [];
    for iter = 1:10000
        append!(y, rng_y[iter][1]);
    end
    Plots.plot!(p1, rng_x, y, lw=4, alpha=0.5, lc="red", legend=false, ls=:dot, dpi=1000)
    xlabel!("time")
    xlims!(0.0, time_max)
    ylabel!("concentration of preys")
    for iter = 1:repetition
        Plots.plot!(p1, times[iter][1:conts[iter]], out[iter][1,1:conts[iter]]./V, lw=1, alpha=1.0, lc="black", legend=false, linetype=:steppre, dpi=1000)
    end
    =#

    p3 = Plots.plot(dpi=1000);
    y = [];
    for iter = 1:10000
        append!(y, rng_y[iter][1]);
    end
    z = [];
    for iter = 1:10000
        append!(z, rng_y[iter][2]);
    end
    Plots.plot!(p3, rng_x, y, lw=3, alpha=0.3, lc="red", label="", ls=:solid, dpi=1000)
    Plots.plot!(p3, rng_x, z, lw=3, alpha=0.3, lc="blue", label="", ls=:solid, dpi=1000)
    Plots.plot!(p3, times[1][1:conts[1]], out[1][1,1:conts[1]]./V, lw=2, alpha=1.0, lc="red", label="Preys", linetype=:steppre, dpi=1000)
    Plots.plot!(p3, times[1][1:conts[1]], out[1][2,1:conts[1]]./V, lw=2, alpha=1.0, lc="blue", label="Predators", linetype=:steppre, dpi=1000)
    xlabel!("time")
    xlims!(0.0, time_max)
    ylabel!("concentration")
    ylims!(0.0, 50.0)
    display(p3)
    Plots.savefig(p3, address*"LotkaVolterra_extinction.png");

    #=
    y = [];
    for iter = 1:10000
        append!(y, rng_y[iter][2]);
    end
    Plots.plot!(p2, rng_x, y, lw=4, alpha=0.5, lc="blue", legend=false, ls=:dot, dpi=1000)
    xlabel!("time")
    xlims!(0.0, time_max)
    ylabel!("concentration of predators")
    for iter = 1:repetition
        Plots.plot!(p2, times[iter][1:conts[iter]], out[iter][2,1:conts[iter]]./V, lw=1, alpha=1.0, lc="black", legend=false, linetype=:steppre)
    end
    display(p1)
    Plots.savefig(p1, address*"LotkaVolterra_Gillespie_10e"*"$esp"*"_"*"prey.png");
    display(p2)
    Plots.savefig(p2, address*"LotkaVolterra_Gillespie_10e"*"$esp"*"_"*"predator.png");
    =#
    save(address*"records25_Gillespie_10e"*"$esp"*".jld", "data", records25)
    save(address*"records5_Gillespie_10e"*"$esp"*".jld", "data", records5)
    save(address*"records75_Gillespie_10e"*"$esp"*".jld", "data", records75)
    save(address*"records10_Gillespie_10e"*"$esp"*".jld", "data", records10)
elseif cmp(algorithm, "TauLeap") == 0 || cmp(algorithm, "TauLeapPostLeapCheck") == 0
    sol = LotkaVolterraDet(c, [1,1], 0.0, time_max);
    rng_x = range(0.0, time_max, 10000);
    rng_y = sol(rng_x);
    #=
    y = [];
    for iter = 1:10000
        append!(y, rng_y[iter][1]);
    end
    Plots.plot!(p1, rng_x, y, lw=4, alpha=0.5, lc="red", legend=false, ls=:dot, dpi=1000)
    xlabel!("time")
    xlims!(0.0, time_max)
    ylabel!("concentration of preys")
    for iter = 1:repetition
        Plots.plot!(p1, times[iter][1:conts[iter]], out[iter][1,1:conts[iter]]./V, lw=1, alpha=1.0, lc="black", legend=false, linetype=:steppre, dpi=1000)
    end
    =#

    p3 = Plots.plot(dpi=1000);
    y = [];
    for iter = 1:10000
        append!(y, rng_y[iter][1]);
    end
    z = [];
    for iter = 1:10000
        append!(z, rng_y[iter][2]);
    end
    Plots.plot!(p3, rng_x, y, lw=4, alpha=0.5, lc="red", label="Preys", ls=:solid, dpi=1000)
    Plots.plot!(p3, rng_x, z, lw=4, alpha=0.5, lc="blue", label="Predators", ls=:solid, dpi=1000)
    Plots.plot!(p3, times[1][1:conts[1]], out[1][1,1:conts[1]]./V, lw=2, alpha=1.0, lc="red", legend=false, linetype=:steppre, dpi=1000)
    Plots.plot!(p3, times[2][1:conts[1]], out[2][1,1:conts[1]]./V, lw=2, alpha=1.0, lc="blue", legend=false, linetype=:steppre, dpi=1000)
    xlabel!("time")
    xlims!(0.0, time_max)
    ylabel!("concentration")
    ylims!(0.0, 50.0)
    display(p3)
    Plots.savefig(p3, address*"LotkaVolterra_extinction.png");
    #=
    y = [];
    for iter = 1:10000
        append!(y, rng_y[iter][2]);
    end
    Plots.plot!(p2, rng_x, y, lw=4, alpha=0.5, lc="blue", legend=false, ls=:dot, dpi=1000)
    xlabel!("time")
    xlims!(0.0, time_max)
    ylabel!("concentration of predators")
    for iter = 1:repetition
        Plots.plot!(p2, times[iter][1:conts[iter]], out[iter][2,1:conts[iter]]./V, lw=1, alpha=1.0, lc="black", legend=false, linetype=:steppre, dpi=1000)
    end
    display(p1)
    #Plots.savefig(p1, address*"LotkaVolterra_"*algorithm*"_"*TLmeth*"_10e"*"$esp"*"_eps="*"$eps"*"prey.png");
    display(p2)
    #Plots.savefig(p2, address*"LotkaVolterra_"*algorithm*"_"*TLmeth*"_10e"*"$esp"*"_eps="*"$eps"*"predator.png");
    =#
    #=
    save(address*"records25_"*algorithm*"_"*TLmeth*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records25)
    save(address*"records5_"*algorithm*"_"*TLmeth*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records5)
    save(address*"records75_"*algorithm*"_"*TLmeth*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records75)
    save(address*"records10_"*algorithm*"_"*TLmeth*"_10e"*"$esp"*"_eps="*"$eps"*".jld", "data", records10)
    =#
end
display(sum(comptimes) / rep)
display(records5)