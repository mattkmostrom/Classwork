using Printf
using Statistics

points = 10
x_min = 0
x_max = 10
dx = (x_max - x_min) / (points - 1)


for n in 1:(points - 1)    #LHR
    suml::Int8 = 0
    x = x_min - dx
    fx = x + dx
    suml = suml + (fx * dx)
end


for n in 1:(points - 1)    #RHR
    sumr::Int8 = 0
    x = x_min + dx
    fx = x + dx
    sumr = sumr + (fx * dx)
end

println("Result from LHR is", suml)
println("Result from RHR is", sumr)

av = (sumr + suml) / 2
println("Their average is", av)
