import SpecialFunctions
using Printf
using DelimitedFiles
x = 0
x0 = 0.            #enter initial postion in angstroms
v0 = 3.            #enter velocity in angstroms/ps
t0 = 0
k = 3
omega = 7         #harmonic frequency
time = 20         #integration time in (periods)
red_mass = 2      #reduced mass in amu
ener = 10         #total energy of the system
npts = 5000       #number of integration points

println("_________________________\n\n")
println("Harmonic Oscillator MD v1.0\n")

#total_energy(x0,v0,omega,rmass,ener) =
red_mass = red_mass * 10. / (6.022169 * 1.380662) #amu to K ps A in units where k_b = 1
fac = sqrt(2. * ener / (red_mass * omega ^ 2))
xmax = fac
xmin = -fac
println("The classical turning points are +/-", xmax)
println("Number of integration points/period: ", npts)
dt = time / (npts - 1)
println("Timestep size: ",dt)

println("")

analytic = zeros(Float64,npts,2)
euler = zeros(Float64,npts,2)
verlet = zeros(Float64,npts,2)
vel_verlet = zeros(Float64,npts,2)

function harmonic(x,t)
  x = x0
  time = t0

  #harmonic analytic variables
  w = sqrt(k/red_mass)
  A = x0
  B = v0 / w

  #Euler variables
	euler_x = 0.
	euler_v = 0.
	prev_euler_x = x0
	prev_euler_v = v0

	#Verlet variables
	verlet_x = 0.
	verlet_v = 0.
	prev_verlet_x = x0

	#velocity Verlet variables
	prev_vel_verlet_x = 0.
	prev_prev_verlet_x = x0
	prev_vel_verlet_v = v0

  for i in 1:npts

    #harmonic analytical answer ***************************************
    harmonic_analytic = A*cos(w*time) + B*sin(w*time)
    analytic[i,1] = time
    analytic[i,2] = harmonic_analytic

		#print(time," ",euler_x," \n")
		#flush(stdout)

    #Euler answer ***************************************
		euler_x = prev_euler_x + prev_euler_v*dt
		euler_v = prev_euler_v + (-w)*prev_euler_x*dt
		prev_euler_x = euler_x
		prev_euler_v = euler_v

		euler[i,1] = time
		euler[i,2] = euler_x

    #Verlet answer ***************************************
		verlet_x = 2*prev_verlet_x - prev_prev_verlet_x + (-w)*prev_verlet_x*dt*dt;
		prev_prev_verlet_x = prev_verlet_x;
		prev_verlet_x = verlet_x;

		verlet[i,1] = time
		verlet[i,2] = verlet_x

    #Velocity Verlet answer ***************************************
		vel_verlet_x = prev_vel_verlet_x + prev_vel_verlet_v*dt + 0.5*(-w)*prev_vel_verlet_x*dt*dt;
		vel_verlet_v = prev_vel_verlet_v + 0.5*((-w)*vel_verlet_x + (-w)*prev_vel_verlet_x)*dt;
		prev_vel_verlet_x = vel_verlet_x;
		prev_vel_verlet_v = vel_verlet_v;

		vel_verlet[i,1] = time
		vel_verlet[i,2] = vel_verlet_x

		#potential[i,1] = x
    #potential[i,2] = vx
    #print(x," ",vx," ")
    #flush(stdout)
    time = time + dt
    #x = x + dx
  end
end

println("***************")
println("Running dyanmics...")

h = harmonic(x0,time)

println("... Done!\n")
println("Writing harmonic results to 'analytic.dat'...")
println("Writing Euler results to 'euler.dat'...")
println("Writing verlet results to 'verlet.dat'...")
println("Writing velocity verlet results to 'vel_verlet.dat'...")

touch("analytic.dat")
outfile = "analytic.dat"
open(outfile, "w") do f
  for i in analytic
    println(f, i)
  end
end
writedlm("analytic.dat",analytic)

touch("euler.dat")
outfile = "euler.dat"
open(outfile, "w") do f
  for i in euler
    println(f, i)
  end
end
writedlm("euler.dat",euler)

touch("verlet.dat")
outfile = "verlet.dat"
open(outfile, "w") do f
  for i in verlet
    println(f, i)
  end
end
writedlm("verlet.dat",verlet)

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
