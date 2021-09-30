h = 0.0000000001 	#the stepsize
x0 = 3. 	#the point to evaluate the derivative
power = 2.

g(x0,power) = pi * x0 ^ power

#forward derivative
on = 0.
xf = x0 + h
on = g(xf,power)
fdx = (on - g(x0,power)) / h


#backward derivative
back = 0.
xb = x0 - h
back = g(xb,power)
bdx = (g(x0,power) - back) / h


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

ana = power * g( 3. , (power - 1) )
println("Analytical answer is: ", ana)
diff = ((central_direct - ana)/central_direct)
println("")
println("% difference between central methods is: ", abs(diff), "%")
