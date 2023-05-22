

#||||----Variables----||||#

R = 1.1 #Basic reproductive rate
vac = 0.00 #proportion of the population vaccinated

tstart = 0 #initial time
tend = 400 #final time
tspan = (tstart, tend) #for passing into DE solver

intervals = 5 #number of points calculated per unit time
timesteps = (tend-tstart)*intervals #Total number of timesteps

susceptibles = 0.98 #initial fraction of population susceptible

b = 10 #per-capita birth rate
c = 10*b #recovery rate from disease


ics = [susceptibles, 1-susceptibles, 0] #initial conditions for DE problem. u_0, v_0, w_0 respectively



#||--Differential equations--||#

#Declaring the system of differential equations.
function pandemicDEs!(du,u,p,t) #p = initial conditions, t = timespan in form of tuple (start,end)
    du[1] = (b/(b+c))*(1-vac-u[1]) - R*u[1]*u[2]
    du[2] = (R*u[1] - 1)*u[2]
    du[3] =  (b/(b+c))*vac + (c/(b+c))*u[2]-(b/(b+c))*u[3]
end#Function


using OrdinaryDiffEq #package for solving the differential equations

#Creates the 'problem' from our declared function
problem = ODEProblem(pandemicDEs!, ics, tspan)

#times to save the solution at
times = LinRange(tstart, tend, Int(timesteps))

#The solution
sol = solve(problem, Tsit5(), saveat=times) #stored as a 2 dimensional array. Access as sol[i,:]. i=1,2,3 give u,v,w respectively.




#||||----Plotting----|||||#

using Plots #The plotting package used here

using LaTeXStrings #Allows latex in plot labels, legends, etc
using ColorSchemes #A library of colour pallets
using Measures #Used here for setting libraries

#In order this sets: size of the plot, font size of tick labels, font size of x and y labels, font size in the ledend, and puts a full border around the graph 
gr(size=(900,500), xtickfontsize=10, ytickfontsize=10, xguidefontsize=16, yguidefontsize=16, legendfontsize=16, framestyle = :box);

#Used for asigning colours to lines in plot by sampling evenly spaced points from a gradient
colorFunc(i, names) = get(ColorSchemes.seaborn_bright,i./length(names))

#Wrapping in a function makes swapping out the plots I'm drawing easier.
function plotConcentrationsWithTime(times, sol, colorFunc, seriesNames)
    #Creates a plot. Moves the legend outside the plot, sets the left, right and top margins, sets the x limits and y limits of the graph.
    myplot = plot(legend =:outertopright, left_margin=8mm, top_margin = 4mm, bottom_margin = 8mm, xlims=[tstart,tend], ylims=[0,1])

    #Draws in the plot for each series. sets the colour, linewidth and name
    for i in eachindex(seriesNames)
        plot!(times, sol[i,:], color = colorFunc(i, seriesNames), linewidth=3, label = seriesNames[i])
    end#for

    #sets the x and y labels
    plot!(xlabel="Time (units?)") #L"" declares a latex string.
    plot!(ylabel= "Proportion of Population")
end#function

#Names of the series that will be plotted. 
seriesNames = [L"u",L"v",L"w"] 
println("Final value: ", sol[2,length(sol[2,:])])
myplot1 = plotConcentrationsWithTime(times,sol,colorFunc,seriesNames)



#=
gr(size=(500,500), xtickfontsize=10, ytickfontsize=10, xguidefontsize=16, yguidefontsize=16, legendfontsize=16, framestyle = :box);

myplot2 = plot(legend =:none, left_margin=8mm, top_margin = 8mm, bottom_margin = 8mm, right_margin=8mm, xlims=[0,1], ylims=[0,1])

plot!(sol[1,:], sol[2,:], linewidth=3)
plot!(xlabel = L"u", ylabel = L"v")
=#

#||||----Investigating rate of change----||||#

function RateOfChangeInfected(array)
    diffArray = ones(length(array))
    for i in eachindex(array)
        diffArray[i] = (R*sol[1,i] - 1)*sol[2,i]
    end#for
    return diffArray
end#function

#println("Final value: ", sol[2,length(sol[2,:])])

#myplot = plot(times, sol[2,:])
#plot!(times, RateOfChangeInfected(sol[2,:]))

diffPlot = plot(times, RateOfChangeInfected(sol[2,:]))

rateRatio = zeros(length(sol[2,:]))
for i in eachindex(sol[2,:])
    rateRatio[i] = (RateOfChangeInfected(sol[2,:])[i]/(abs(RateOfChangeInfected(sol[2,:])[i])))*log10( abs(sol[2,i]/RateOfChangeInfected(sol[2,:])[i]))
end#for
myplot2 = plot(times, rateRatio, legend=:none)
