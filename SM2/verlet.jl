import SpecialFunctions
using Printf
using DelimitedFiles
x = 0
x0 = 1.            #enter initial postion in angstroms
vx = 0.
v0 = 3.            #enter velocity in angstroms/ps
t0 = 0
omega = 7         #harmonic frequency
time = 7         #integration time in (periods)
red_mass = 2      #reduced mass in amu
ener = 10         #total energy of the system
npts = 50       #number of integration points

println("_________________________\n\n")
println("Harmonic Oscillator MD v0.1\n")

#total_energy(x0,v0,omega,rmass,ener) =
red_mass = red_mass * 10. / (6.022169 * 1.380662) #amu to K ps A in units where k_b = 1
fac = sqrt(2. * ener / (red_mass * omega ^ 2))
xmax = fac
xmin = -fac
println("The classical turning points are +/-", xmax)
println("")
println("***************")

println("Number of integration points/period: ", npts)
dt = time / (npts - 1)
range = [1:dx:npts]
println("Integration stepsize is: ", dx)
println("")

potential = zeros(Float64,npts,2)
function harmonic(x)
  x = x0
  for i in 1:npts
    vx = 0.5*red_mass*omega*omega*x*x
    potential[i,1] = x
    potential[i,2] = vx
    #print(x," ",vx," ")
    #flush(stdout)
    x = x + dx
  end
end

println("Running dyanmics on harmonic oscillator...\n")

h = harmonic(x0)

println("... Done!\n")
println("Writing harmonic results to 'potential.dat'...\n")

touch("potential.dat")
outfile = "potential.dat"
open(outfile, "w") do f
  for i in potential
    println(f, i)
  end
end
writedlm("potential.dat",potential)

println("... Done!")
