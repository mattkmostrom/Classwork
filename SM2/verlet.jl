x0 = 0            #enter initial postion in angstroms
v0 = 3            #enter velocity in angstroms/ps
omega = 7         #harmonic frequency
time = 10         #integration time in (periods)
red_mass = 2      #reduced mass in amu
ener = 10         #total energy of the system
npts = 5000       #number of integration points

#total_energy(x0,v0,omega,rmass,ener) =
red_mass=red_mass * 10. / 6.022169 / 1.380662 #amu to K ps A in units where k_b = 1
fac=sqrt( 2. * ener / red_mass / omega ^ 2)
xmax=fac
xmin=-fac
println("The classical turning points are +/-", xmax)
println("")
println("***************")

println("Number of integration points/period: ", npts)
dx = (xmax - xmin) / (npts - 1)
range = [1:dx:npts]
println("Integration stepsize is: ", dx)
println("")

f = open("potential.dat","w")



#write(f, whatever data)
close(f)
