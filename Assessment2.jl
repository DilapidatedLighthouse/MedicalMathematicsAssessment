using Plots, SciMLBase 

#||||----Variables----||||#

R = 10.0 #Basic reproductive rate
vac = 0.00 #proportion of the population vaccinated

tstart = 0 #initial time
tend = 10 #final time
tspan = (tstart, tend)

intervals = 20 #number of points calculated per unit time
timesteps = (tend-tstart)*intervals #Total number of timesteps

susceptibles = 0.98 #initial fraction of population susceptible

b = 10 #per-capita birth rate
c = 10*b #recovery rate of disease


ics = [susceptibles, 1-susceptibles, 0] #initial conditions for DE problem. u_0, v_0, w_0 respectively



#||--Differential equations--||#


function pandemicDEs!(du,u,p,t) #p = initial conditions, t = timespan in form of tuple (start,end)
    du[1] = (b/(b+c))*(1-vac-u[1]) - R*u[1]*u[2]
    du[2] = (R*u[1] - 1)*u[2]
    du[3] =  (b/(b+c))*vac + (c/(b+c))*u[2]-(b/(b+c))*u[3]
end#Function

using OrdinaryDiffEq

problem = ODEProblem(pandemicDEs!, ics, tspan)

times = LinRange(tstart, tend, Int(timesteps))

sol = solve(problem, Tsit5(), saveat=times) #stored as a 2 dimensional array. Access as sol[i,:]. i=1,2,3 give u,v,w respectively.




#||||----Plotting----|||||#

using LaTeXStrings
using ColorSchemes
using Statistics
using Measures

gr(size=(900,500), xtickfontsize=10, ytickfontsize=10, xguidefontsize=16, yguidefontsize=16, legendfontsize=16, dpi=100, framestyle = :box, grid = :xy);

colorFunc(i, names) = get(ColorSchemes.seaborn_bright,i./length(names))
seriesNames = [L"u",L"v",L"w"] 

myplot = plot(legend =:outertopright, left_margin=8mm, top_margin = 4mm, bottom_margin = 8mm, xlims=[tstart,tend], ylims=[0,1])

for i in 1:length(seriesNames)
    plot!(times, sol[i,:], color = colorFunc(i, seriesNames), linewidth=3, label = seriesNames[i])
end#for


plot!(xlabel="Time (units?)") #L"" declares a latex string.
plot!(ylabel= "Proportion of Population")