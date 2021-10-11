import SpecialFunctions
points = 100000
ints = (points - 1)             #intervals between points, the things actually being integrated over
x0 = 0
x_min = 0
x_max = 5
dx = (x_max - x_min) / ints
range = x_min:dx:x_max          # start:iteration:finish

alpha = 3
sigma = sqrt(2*alpha)

g(x) = (1 / sqrt(2 * pi * sigma ^ 2)) * exp(-alpha * (x-x0) ^ 2)            # Change your function here

#k = 2 / sqrt(pi / alpha)
#l(x) = k * exp(-(x-x0) ^ 2)

lhr(f,range) = sum(f.(range[1:end-1]) .* step(range)) #f operating on each element of range and evaluating at each range
rhr(f,range) = sum(f.(range[2:end])   .* step(range))

println("")
println("*************")
println("Result from LHR is ", lhr(g,range))
println("Result from RHR is ", rhr(g,range))

av = (lhr(g,range) + rhr(g,range)) / 2
println("Their average is ", av)
println("")

erf = SpecialFunctions.erf
N = sqrt(pi/alpha) / 2
answer = (N * erf(sqrt(alpha) * x_max) - N * erf(sqrt(alpha) * (x_min)))
#answer = N * (erf(x_max) - erf(x_min))
println("Analytical answer: ", answer) #https://www.wolframalpha.com/input/?i=Integrate%5Bexp%28-a*x%5E2%29%5D
