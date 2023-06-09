###||||----Medical Mathematics Assignment 2: Code Re-write----||||###

    #||||----Variables----||||#

saveLocation = "C:/Users/joemu/Desktop/Uni Stuff/Medical Mathematics/MedicalMathematicsAssessment/Plots" #file path of folder where figures will be saved
#--Model Parameters--#
R = 10 #Basic reproductive rate
vac = 0.00 #proportion of the population vaccinated
b = 10 #per-capita birth rate
c = 10*b #recovery rate from disease

#--Solution Parameters--#
tstart = 0 #initial time
tend = 20 #final time
tspan = (tstart, tend) #for passing into DE solver

intervals = 5 #number of points calculated per unit time
timesteps = (tend-tstart)*intervals #Total number of timesteps

#--Initial Conditions--#
susceptibles = 0.98 #initial fraction of population susceptible
ics = [susceptibles, 1-susceptibles, 0] #initial conditions for DE problem. u_0, v_0, w_0 respectively



    #||||----Declaring and Solving----||||#

#Declaring the system of differential equations.
function pandemicDEs!(du,u,p,t) #p = initial conditions, t = timespan in form of tuple (start,end)
    du[1] = (b/(b+c))*(1-vac-u[1]) - R*u[1]*u[2]
    du[2] = (R*u[1] - 1)*u[2]
    du[3] =  (b/(b+c))*vac + (c/(c+b))*u[2]-(b/(b+c))*u[3]
end#Function


using OrdinaryDiffEq #package for solving the differential equations

#Uses the previously declared function to construct a 'problem' for use with the OrdinaryDiffEq library
problem = ODEProblem(pandemicDEs!, ics, tspan)

#times to save the solution at
times = LinRange(tstart, tend, Int(timesteps))

#The solution
sol = solve(problem, Tsit5(), saveat=times) #stored as a 2 dimensional array. Access as sol[i,:]. i=1,2,3 give u,v,w respectively.
#Tsit5 is a solution algorithm. Specifically an efficient variant of Runge-Kutta



    #||||----Plotting----|||||#

using Plots #The plotting package used here

using LaTeXStrings #Allows latex in plot labels, legends, etc
using ColorSchemes #A library of colour pallets
using Measures #Used here for setting libraries

#In order this sets: size of the plot, font size of tick labels, font size of x labels, font size of y labels, font size in the ledend, and puts a full border around the graph, then sets the margins of the graph
gr(xlabel="Time", ylabel= "Proportion of Population", size=(900,500), xtickfontsize=10, ytickfontsize=10, xguidefontsize=16, yguidefontsize=16, legendfontsize=16, framestyle = :box, left_margin=8mm, top_margin = 4mm, bottom_margin = 8mm,);

#Used for asigning colours to lines in plot by sampling evenly spaced points from a gradient
colorFunc(i, names) = get(ColorSchemes.seaborn_bright,i./length(names))

#Wrapping in a function makes swapping out the plots I'm drawing easier.
function plotConcentrationsWithTime(times, sol, colorFunc, seriesNames)
    #Creates a plot. Moves the legend outside the plot, sets the left, right and top margins, sets the x limits and y limits of the graph.
    myplot = plot(legend =:outertopright, left_margin=8mm, top_margin = 4mm, bottom_margin = 8mm, xlims=[tstart,tend], ylims=[0,1])

    #Draws in the plot for each series. sets the colour, linewidth and name
    for i in eachindex(seriesNames)
        plot!(times, sol[i,:], color = colorFunc(i, seriesNames), linewidth=6, label = seriesNames[i], linealpha = 0.6)
    end#for

    #sets the x and y labels
    plot!(xlabel="Time") #L"" declares a latex string.
    plot!(ylabel= "Proportion of Population")
    return myplot
end#function
function plotConcentrationsWithTime!(times, sol, colorFunc, seriesNames)
    #Creates a plot. Moves the legend outside the plot, sets the left, right and top margins, sets the x limits and y limits of the graph.
    plot!(legend =:outertopright, left_margin=8mm, top_margin = 4mm, bottom_margin = 8mm, xlims=[tstart,tend], ylims=[0,1])

    #Draws in the plot for each series. sets the colour, linewidth and name
    for i in eachindex(seriesNames)
        plot!(times, sol[i,:], color = colorFunc(i, seriesNames), linewidth=4, label = seriesNames[i], linestyle=:dot, linealpha=1)
    end#for

    #sets the x and y labels
    plot!(xlabel="Time") #L"" declares a latex string.
    plot!(ylabel= "Proportion of Population")
end#function

#Names of the series that will be plotted. 
seriesNames = [L"u",L"v",L"w"] 
#println("Final values: ", " u = ", sol[1,length(sol[1,:])], " v = ", sol[2,length(sol[2,:])], " w = ", sol[3,length(sol[3,:])]) #bug testing
#myplot1 = plotConcentrationsWithTime(times,sol,colorFunc,seriesNames)
#savefig(string(saveLocation, "/test.png"))

#||--Plots for Assignment--||#
#All methods in this section are also present further up in the code, unless otherwise specified



#--No Vaccination model--#

# R = 10 #
R=10

problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tstart, tend, Int(timesteps))
sol = solve(problem, Tsit5(), saveat=times) 
plot1 =  plotConcentrationsWithTime(times,sol,colorFunc,seriesNames)
savefig(string(saveLocation, "/noVaccR10.png"))

# R = 0.2 #
R=0.2

problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tstart, tend, Int(timesteps))
sol = solve(problem, Tsit5(), saveat=times) 
plot1 =  plotConcentrationsWithTime(times,sol,colorFunc,seriesNames)
savefig(string(saveLocation, "/noVaccR02.png"))


# R = 1 and R = 1.01 # 

#Update R values
R=1

problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tstart, tend, Int(timesteps))
sol1 = solve(problem, Tsit5(), saveat=times)

R=0.99

problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tstart, tend, Int(timesteps))
sol2 = solve(problem, Tsit5(), saveat=times)
seriesNames = [L"u: R=1",L"v: R=1",L"w: R=1"] 
plot1 =  plotConcentrationsWithTime(times,sol1,colorFunc,seriesNames)
colorFunc(i, names) = get(ColorSchemes.glasbey_bw_minc_20_hue_150_280_n256,i./length(names))
seriesNames = [L"u: R=0.99",L"v: R=0.99",L"w: R=0.99"] 
plotConcentrationsWithTime!(times,sol2,colorFunc,seriesNames)
savefig(string(saveLocation, "/noVaccR1R01.png"))



# R = 1 and R = 1.01 close up #
tstart = 0 #initial time
tend = 800 #final time
tspan = (tstart, tend) #for passing into DE solver

seriesNames = [L"v: R=1", L"v: R=1.01", L"v: R=0.99"]

R=1
problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tend-20, tend, Int(timesteps))
sol1 = solve(problem, Tsit5(), saveat=times)

R=1.01
problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tend-20, tend, Int(timesteps))
sol2 = solve(problem, Tsit5(), saveat=times)

R=0.99
problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tend-20, tend, Int(timesteps))
sol3 = solve(problem, Tsit5(), saveat=times)


colorFunc(i, names) = get(ColorSchemes.seaborn_bright,i./length(names))
plot1 =  plot(times, sol1[2,:], color = colorFunc(1, seriesNames), linewidth=4, label = seriesNames[1], linealpha=1)
plot!(times, sol2[2,:], color = colorFunc(2, seriesNames), linewidth=4, label = seriesNames[2], linealpha=1)
plot!(times, sol3[2,:], color = colorFunc(3, seriesNames), linewidth=4, label = seriesNames[3], linealpha=1, legend =:outertopright)
plot!(left_margin = 8mm, bottom_margin = 8mm)
savefig(string(saveLocation, "/noVaccR1R01End.png"))


#--Vaccination Model--#
# set R=3
R = 3 #Basic reproductive rate
b = 10 #per-capita birth rate
c = 10*b #recovery rate from disease

#--Solution Parameters--#
tstart = 0 #initial time
tend = 40 #final time
tspan = (tstart, tend) #for passing into DE solver
# vac = 0.5 #
vac=0.5
seriesNames = [L"u",L"v",L"w"] 
problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tstart, tend, Int(timesteps))
sol = solve(problem, Tsit5(), saveat=times) 
plot1 =  plotConcentrationsWithTime(times,sol,colorFunc,seriesNames)
savefig(string(saveLocation, "/Vacc05R3.png"))

# vac = 0.8 #
vac=0.8
R=10
seriesNames = [L"u",L"v",L"w"] 
problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tstart, tend, Int(timesteps))
sol = solve(problem, Tsit5(), saveat=times) 
plot1 =  plotConcentrationsWithTime(times,sol,colorFunc,seriesNames)
savefig(string(saveLocation, "/Vacc08R10.png"))


#Critical vac values
R=3
tstart = 0 #initial time
tend = 800 #final time
tspan = (tstart, tend) #for passing into DE solver
vac = 0.66
seriesNames = [L"v: p = 0.66",L"v: p = 0.67", L"v: p = 0.68"] 
problem = ODEProblem(pandemicDEs!, ics, tspan)
times = LinRange(tend-20, tend, Int(timesteps))
sol = solve(problem, Tsit5(), saveat=times) 
vac = 0.67
problem = ODEProblem(pandemicDEs!, ics, tspan)
sol2 = solve(problem, Tsit5(), saveat=times) 
vac = 0.68
problem = ODEProblem(pandemicDEs!, ics, tspan)
sol3 = solve(problem, Tsit5(), saveat=times) 
plot1 =  plot(times, sol[2,:], color = colorFunc(1, seriesNames), linewidth = 4, label = seriesNames[1], linealpha = 1, legend = :outertopright)
plot!(times, sol2[2,:], color = colorFunc(2, seriesNames), linewidth = 4, label = seriesNames[2], linealpha = 1)
#savefig(string(saveLocation, "/VaccR066R067End.png"))




#Critical vac values zoomed in

#--Examples--#

#Disease 1

#Disease 2





#=
#||--Plotting u vs v--||#
gr(size=(500,500), xtickfontsize=10, ytickfontsize=10, xguidefontsize=16, yguidefontsize=16, legendfontsize=16, framestyle = :box);

myplot2 = plot(legend =:none, left_margin=8mm, top_margin = 8mm, bottom_margin = 8mm, right_margin=8mm, xlims=[0,1], ylims=[0,1])

plot!(sol[1,:], sol[2,:], linewidth=3)
plot!(xlabel = L"u", ylabel = L"v")
=#

    #||||----Investigating rate of change----||||#
#=
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
=#