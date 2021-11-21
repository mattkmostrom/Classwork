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

abstract type Element end

struct Ar <: Element 
    type::String
    mass::Float64
    eps::Float64
    sig::Float64
    pos_x::Float64
    pos_y::Float64
    pos_z::Float64
    vel_x::Float64
    vel_y::Float64
    vel_z::Float64
end

struct He <: Element  
    type::String
    mass::Float64
    eps::Float64
    sig::Float64
    pos_x::Float64
    pos_y::Float64
    pos_z::Float64
    vel_x::Float64
    vel_y::Float64
    vel_z::Float64
end

type(c::Ar) = "Ar"
mass(c::Ar) = "39.95"
eps(c::Ar) = "128.326802"
sig(c::Ar) = "3.371914"
pos_x(c::Ar) = c.pos_x
pos_y(c::Ar) = c.pos_y
pos_z(c::Ar) = c.pos_z
vel_x(c::Ar) = c.vel_x
vel_y(c::Ar) = c.vel_y
vel_z(c::Ar) = c.vel_z

n = 5
array_range = 1:1:n
atoms = zeros(Float64,n,10)
function cubic_lattice(n,dx)

  for i in 1:n
    for j in 1:n
      for k in 1:n
      atom = Ar(i*dx,j*dx,k*dx,0,0,0) 
      print(atom)
      flush(stdout)
      append!(atoms,atom)
      end #i
    end   #j
  end     #k
  return atoms
  #print(atoms)
  #flush(stdout)
end       #function

cube = cubic_lattice(5,3)

print(atoms)



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
