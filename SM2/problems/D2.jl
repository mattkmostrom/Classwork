import SpecialFunctions
h = 0.0001 	#the stepsize
x = 3. 	#the point to evaluate the derivative
power = 2.
alpha = 1
sigma = sqrt(2*alpha)

erf = SpecialFunctions.erf
g(x) = erf(x)

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
ana = 2 * alpha * exp(- (alpha^2) * (x^2)) / sqrt(pi)
println("Analytical answer is: ", ana)
diff = ((central_direct - ana)/central_direct)
println("")
println("% difference between central methods is: ", abs(diff), "%")
