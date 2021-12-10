import SpecialFunctions
using Printf
using DelimitedFiles
x = 0.
x0 = 0.            	#enter initial postion in angstroms
v0 = 1.            	#enter velocity in angstroms/ps
t0 = 0.

global k = 3.
global red_mass = 2.      	#reduced mass in amu
global w = sqrt(k/red_mass)
global T = 2. * pi/w       	#time per period
global time = T*100.     	#integration time;
global npts = 100000       	#number of integration points

ener = 10.
dt = 0.01
range = 1:dt:time


println("_________________________\n\n")
println("Harmonic Oscillator MD v1.0\n")

red_mass *= 10. / (6.022169 * 1.380662) #amu to K ps A in units where k_b = 1
fac = sqrt(2. * ener / (red_mass * w ^ 2))
xmax = fac
xmin = -fac

println("The classical turning points are +/-", xmax)
println("Number of integration points/period: ", npts)
println("Timestep size: ",dt)
println("")

function get_acc(x,w,m)
  dudx = m * (w^2) * (x-x0)   #U = (m/2) * (w^2) * (x-x0)^2
  F = -dudx
  acc = F/m

  return acc
end


analytic = zeros(Float64,npts,2)
euler = zeros(Float64,npts,2)
verlet = zeros(Float64,npts,2)
vel_verlet = zeros(Float64,npts,2)

function harmonic(x0,time,range)
  #time = t0
  acc = 0.


  #harmonic analytic variables
  A = x0
  B = v0 / w

  #Euler variables
  euler_x = 0.
  euler_v = 0.
  prev_euler_x = x0
  prev_euler_v = v0

  #Verlet variables
  verlet_x = 0.
  prev_verlet_x = 0. #x0
  prev_prev_verlet_x =  x0 - v0*dt

  #velocity Verlet variables
  vel_verlet_x = 0.
  prev_vel_verlet_x = x0
  prev_vel_verlet_v = v0

  for i in 1:npts

    #harmonic analytical answer ***************************************
    harmonic_analytic = A * cos(w * time) + B * sin(w * time)
    analytic[i,1] = time
    analytic[i,2] = harmonic_analytic #position


    #Euler answer ***************************************
	euler_x = prev_euler_x + prev_euler_v * dt
	euler_v = prev_euler_v + (-w^2) * prev_euler_x * dt
	prev_euler_x = euler_x
	prev_euler_v = euler_v

	euler[i,1] = time
	euler[i,2] = euler_x

    #Verlet answer ***************************************
    acc = get_acc(verlet_x,w,red_mass)
	  #verlet_x = 2 * prev_verlet_x - prev_prev_verlet_x + (-w) * prev_verlet_x * dt^2
    verlet_x = 2 * prev_verlet_x - prev_prev_verlet_x + acc * dt^2
	  prev_prev_verlet_x = prev_verlet_x
	  prev_verlet_x = verlet_x

	  #@printf("Step: %i, verlet_x: %f, prev_verlet_x: %f, prev_prev_verlet_x: %f, acceleration: %f\n",i,verlet_x,prev_verlet_x,prev_prev_verlet_x,acc)
	  verlet[i,1] = time
	  verlet[i,2] = verlet_x

    #Velocity Verlet answer ***************************************
    acc = get_acc(vel_verlet_x,w,red_mass)
	#vel_verlet_x = prev_vel_verlet_x + prev_vel_verlet_v * dt + 0.5 * (-w) * prev_vel_verlet_x * dt^2
    vel_verlet_x = prev_vel_verlet_x + prev_vel_verlet_v * dt + 0.5 * acc * dt^2
	vel_verlet_v = prev_vel_verlet_v + 0.5 * ((-w^2) * vel_verlet_x + (-w^2) * prev_vel_verlet_x) * dt
	prev_vel_verlet_x = vel_verlet_x
	prev_vel_verlet_v = vel_verlet_v

  #@printf("Step: %i, vel_verlet_x: %f, prev_vel_verlet_x: %f, acceleration: %f\n",i,vel_verlet_x,prev_vel_verlet_x,acc)
	vel_verlet[i,1] = time
	vel_verlet[i,2] = vel_verlet_x

    error = verlet_x - vel_verlet_x
    println("verlet error: ",error)

    time += dt
    #x = x + dx
  end
end

println("***************")
println("Running dyanmics...")

harmonic(x0,time,range)

println("... Done!\n")

println("Writing harmonic results to 'analytic.dat'...")
touch("analytic.dat")
outfile = "analytic.dat"
open(outfile, "w") do f
  for i in analytic
    println(f, i)
  end
end
writedlm("analytic.dat",analytic)

println("Writing Euler results to 'euler.dat'...")
touch("euler.dat")
outfile = "euler.dat"
open(outfile, "w") do f
  for i in euler
    println(f, i)
  end
end
writedlm("euler.dat",euler)

println("Writing verlet results to 'verlet.dat'...")
touch("verlet.dat")
outfile = "verlet.dat"
open(outfile, "w") do f
  for i in verlet
    println(f, i)
  end
end
writedlm("verlet.dat",verlet)

println("Writing velocity verlet results to 'vel_verlet.dat'...")
touch("vel_verlet.dat")
outfile = "vel_verlet.dat"
open(outfile, "w") do f
  for i in vel_verlet
    println(f, i)
  end
end
writedlm("vel_verlet.dat",vel_verlet)
println("... Done!")

println("")
println("xmgrace analytic.dat euler.dat verlet.dat vel_verlet.dat")
