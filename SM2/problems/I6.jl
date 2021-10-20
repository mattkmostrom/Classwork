#using Pkg;Pkg.add("Plots");Pkg.add("PyPlot")

import SpecialFunctions
points = 500000
ints = (points - 1) #intervals between points, the things actually being integrated over
x_min = 0.001
x_max = 1000            #the error seems to go up with how big you make x_max
dx = (x_max - x_min) / ints
range = x_min:dx:x_max  # start:iteration:finish

ge(x) = exp(-2 * x) / x
gs(x) = sin(x) / x

lhr(f,range) = sum(f.(range[1:end-1]) .* step(range)) #f operating on each element of range and evaluating at each range
rhr(f,range) = sum(f.(range[2:end])   .* step(range))

println("")
println("____________________")
println("")
println("Result for Exp from LHR is ", lhr(ge,range))
println("Result for Exp from RHR is ", rhr(ge,range))

av = (lhr(ge,range) + rhr(ge,range)) / 2
println("Trap rule average for Exp(-2x)/x is ", av)
exp_answer = 5.63939        #Ei(x) wolfram
println("Analytical Ei(x): ", exp_answer)
println("")

println("**")
println("Result for Sin from LHR is ", lhr(gs,range))
println("Result for Sin from RHR is ", rhr(gs,range))

av = (lhr(gs,range) + rhr(gs,range)) / 2
println("Trap rule average for sin(x)/x is ", av)
sin_answer = 1.56923
println("Analytical Si(x): ", sin_answer, " which is close to pi/2 (1.57090)")
println("")


#matrix = zeros(2,length(range))
#f = open("I6_exp.dat","w")
#close(f)

#f = open("I6_sin.dat","w")
#close(f)

#using Plots
#plot(rhr(gs,range))
