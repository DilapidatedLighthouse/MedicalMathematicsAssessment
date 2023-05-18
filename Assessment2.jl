using Plots, SciMLBase 

#||||----Variables----||||#

R = 10.0 #Basic reproductive rate
vac = 0.00 #proportion of the population vaccinated

tstart = 0 #initial time
tend = 10 #final time
tspan = (tstart, tend)

interval = 5.0 #number of points calculated per unit time
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

sol = solve(problem, Tsit5()) #stored as a 2 dimensional array. Access as sol[i,:]. i=1 gives time, i=2,3,4 give u,v,w respectively.

times = LinRange(tstart, tend, length(sol[1,:]))


myplot = plot(times, sol[1,:])
myplot = plot!(times, sol[2,:])