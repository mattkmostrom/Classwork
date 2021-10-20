points = 10000
ints = (points - 1) #intervals between points, the things actually being integrated over
x_min = 0
x_max = 10
dx = (x_max - x_min) / ints
range = x_min:dx:x_max  # start:iteration:finish

g(x) = x*x*x            # Change your function here

function lhr()
    suml = 0
    x = x_min
    for n in 1:ints    #LHR
        fx = x*x
        suml = suml + (fx * dx)
        x = x + dx
    end
    return suml
end

function rhr()
    sumr = 0
    x = x_min + dx
    for n in 1:ints    #RHR
        fx = x*x
        sumr = sumr + (fx * dx)
        x = x + dx
    end
    return sumr
end


println("")
println("*************")
println("Result from LHR is ", lhr(g,range))
println("Result from RHR is ", rhr(g,range))

av = (lhr(g,range) + rhr(g,range)) / 2
println("Their average is ", av)
println("")

answer = ((x_max*x_max*x_max*x_max)/4) - ((x_min*x_min*x_min*x_min)/4)
println("Analytical answer: ", answer)
