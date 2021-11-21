h = 0.0001 	#the stepsize
x = 0.2 	#the point to evaluate the derivative

alpha = 5.
sigma = 1 / sqrt(2 * alpha)
g(x) = (1 / sqrt(2 * pi * sigma ^ 2)) * exp(-alpha * (x) ^ 2)

#power = 2.
#g(x0,power) = pi * x0 ^ power

#forward derivative
on = 0.
xf = x + h
on = g(xf)
fdx = (on - g(x)) / h


#backward derivative
back = 0.
xb = x - h
back = g(xb)
bdx = (g(x) - back) / h


central_direct = ((on - back) / h) / 2.
central_average = (fdx + bdx) / 2.
println("")
println("***************")

println("Forward: ", fdx)
println("Backward: ", bdx)
println("")

println("Central Direct: ", central_direct)
println("Central Average: ", central_average)
println("")

#ana = power * g( 3. , (power - 1) )
ana = -2 * alpha * x * g(x)
println("Analytical answer is: ", ana)
diff = ((central_direct - ana)/central_direct)
println("")
println("% difference between central methods is: ", abs(diff), "%")
