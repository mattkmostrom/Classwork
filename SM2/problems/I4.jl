import SpecialFunctions
using Printf
using DelimitedFiles
x0 = 0
x_min = 0
x_max = 5
alpha = 3
sigma = 1 / sqrt(2 * alpha)
println("")
println("*******************")
println("This program is an answer problem I4, and will output a 4-column file called 'I4_error.dat', which will graph rectangle rule and trapezoidal rule integration error against dx, respectively.")
println("")
@printf("First Moment: %f \nx-min: %f \nx-max: %f \nalpha : %f \nsigma : %f", x0,x_min,x_max,alpha,sigma)
println("")

g(x) = (1 / sqrt(2 * pi * sigma^2)) * exp(-alpha * (x - x0)^2)            # Change your function here

#f operating on each element of range and evaluating at each range
lhr(f, range) = sum(f.(range[1:end-1]) .* step(range))
rhr(f, range) = sum(f.(range[2:end]) .* step(range))


array_range = 1:1000
err_array = zeros(Float64,length(array_range),3)

for i in array_range
    dx = (x_max - x_min) / i
    range = x_min:dx:x_max

    av = (lhr(g, range) + rhr(g, range)) / 2

    erf = SpecialFunctions.erf
    N = sqrt(pi / alpha) / 2
    answer = N * erf(sqrt(alpha) * x_max)
    #https://www.wolframalpha.com/input/?i=Integrate%5Bexp%28-a*x%5E2%29%5D

    rect_err = ((lhr(g, range) - answer) / answer) * 100
    trap_err = ((av - answer) / answer) * 100

    err_array[i,1] = dx
    err_array[i,2] = rect_err
    err_array[i,3] = trap_err

    println("dx: ",dx)
    println("rect_err: ",rect_err," %")
    println("trap_err: ",trap_err," %")
    println("")
end                 #i loop

#print(err_array)

touch("I4_error.dat")
outfile = "I4_error.dat"
open(outfile, "w") do f
  for i in err_array
    println(f, i)
  end
end
writedlm("I4_error.dat","#dx rectangle_%_error trap._%_error")
writedlm("I4_error.dat",err_array)

# cat I4_error.dat | awk '{print($1,$2)}' > rect_err.dat;cat I4_error.dat | awk '{print($1,$3)}' > trap_err.dat; xmgrace rect_err.dat trap_err.dat
