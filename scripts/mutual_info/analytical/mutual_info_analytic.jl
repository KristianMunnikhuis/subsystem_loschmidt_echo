"""Analytical calculations of TFIM quantities"""

###Data Reading and Writing
base_folder = dirname(@__DIR__);
plot_dir = joinpath(@__DIR__,"Plots")
mkpath(plot_dir)


#TFIM Functions 
include(joinpath(base_folder,"TFIM_Functions/TFIM_functions.jl"));
x=  [sigma_general([1,2,3,4],ti) for ti in t]
y = sqrt.(abs.([sigma_general([1],ti) for ti in t]))

y[10:50] *= -1
plot(y)
plot!(sqrt.(real.(x)))
####PREFACTORS NEED WORK!!!!!
P_n(4,0)
P_n(2,0.5)

#Parameters
L = 40
J = 1
h1 = 0
h2 = 1.5
h_i = [i for i in 0:.01:2]
P10 = [log(P_n(2*4,0.0)/(P_n(4,0.0)^2))]
for hi in h_i[2:end]
    h1 = hi
    push!(P10,log(P_n(8,0.0)/P_n(4,0.0)^2))
end
P10
plot(h_i,abs.(P10))
#Used for plotting
T_final = 5
t = [i for i in 0:0.05:T_final];
k = [2*pi*(n+1/2)/L for n in 0:L-1];

timer = time()
P10 = [[log(P_n(2*ni,ti)/(P_n(ni,ti)^2)) for ti in t] for ni in 1:3];

print(timer-time())
plot(t,real.(P10))

