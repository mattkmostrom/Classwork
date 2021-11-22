import SpecialFunctions
using Printf
using DelimitedFiles
t0 = 0
k = 3
omega = 7         #harmonic frequency
time = 20         #integration time in (periods)
red_mass = 2      #reduced mass in amu
ener = 10         #total energy of the system
npts = 5000       #number of integration points
dx = 3.

println("________________________\n\n")
println("Morse Potential MD v0.1\n")


red_mass = red_mass * 10. / (6.022169 * 1.380662) #amu to K ps A in units where k_b = 1
println("Number of integration points/period: ", npts)
dt = time / (npts - 1)
println("Timestep size: ",dt)

println("")

#****************************
#Make the cubic lattice
#****************************

struct Element      #Make classes a thing in julia
    type::String
    mass::Real
    eps::Real
    sig::Real
    pos::Vector{<:Real}
    vel::Vector{<:Real}
end

Ar(pos::Vector,vel::Vector) = Element("Ar",39.95,128.326802,3.371914,pos,vel) #define some atom parameters

n = 5
dx = 3
atoms = []

for i in 0:n    #Make a cubic lattice
  for j in 0:n
    for k in 0:n
      
      atom = Ar([i,j,k] .* dx,zeros(3)) #Ar is a function of 2 vectors: pos & vel 
      #atom = Ar([i,j,k] .* dx,rand(3)) #again but random velocities

      push!(atoms,atom)
    end #i
  end   #j
end     #k
#print(atoms)

#****************************
#Make the velocities gaussian
#****************************

for i in atoms    
  itemp = 0
  itemp = 



end






println("")
println("***************")
println("Running dynamics...")


println("... Done!\n")

println("Writing Trajectory to 'out.dat'")

#touch("out.dat")
#outfile = "out.dat"
#open(outfile, "w") do f
#  for i in traj
#    println(f, i)
#  end
#end
#writedlm("out.dat",traj)

println("... Done!")


println("")
