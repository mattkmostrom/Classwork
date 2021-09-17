h = 0.1
x0 = 3

g(x) = pi * x0

on = 0.
function forward(f)
	xf = x0 + h
	on = xf * pi
 	fdx = (on - f(x0)) / h
end

back = 0.
function backward(f)
	xb = x0 - h
	back = xb * pi
	bdx = (f(x0) - back) / h
end

central_direct = ((on - back) * h) / 2.

central_average = (fdx + bdx) / 2.


println("")
println("***************")
println("Forward: ", forward(g))
println("Backward: ", backward(g))


ana = pi
println("Analytical answer is: ", ana)
diff = ((central_direct(g) - ana)/central_direct(g))
println("")
println("% difference is: ", diff)
