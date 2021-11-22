#To do:
#1.) Make a box
#2.)
#
import SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
println("________________________\n\n")
println("Morse Potential MD v0.1\n")

global nsteps = 5       #number of integration points
global dt = 0.001
println("Number of integration points/period: ", nsteps)
println("Timestep size: ",dt)
println("")


mutable struct Element      #Make classes a thing in julia
    type::String
    mass::Real
    eps::Real
    sig::Real
    w::Real
    x0::Vector{<:Real}
    pos::Vector{<:Real}
    vel::Vector{<:Real}
    acc::Vector{<:Real}
    acc_old::Vector{<:Real}
end

Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",39.95,128.326802,3.371914,1,pos,pos,vel,acc,acc_old) #define some atom parameters

#****************************
#Make the cubic lattice
#****************************

n = 2
dx = 3.
atoms = []

for i in 1:n    #Make a cubic lattice
  for j in 1:n
    for k in 1:n
      #atom = Ar([i,j,k] .* dx,zeros(3),zeros(3),zeros(3)) #Ar is a function of 2 vectors: pos & vel
      atom = Ar([i,j,k] .* dx,[i,j,k] .* dx,rand(3),zeros(3),zeros(3)) #again but random velocities
      push!(atoms,atom)
    end #i
  end   #j
end     #k


#****************************
#Make the velocities gaussian
#****************************

#target_temp = 25
function gaussify(atoms,target_temp)
  for i in 1:length(atoms)
    #target_temp = 25
    itemp = 0
    itemp *= atoms[i].mass * dot(atoms[i].vel,atoms[i].vel)
    itemp /= 3 * length(atoms) - 3
    atoms[i].vel *= sqrt(target_temp/itemp)
  end
end


function forces(atoms)
  for i in 1:length(atoms)
    atoms[i].acc = -(atoms[i].w)^2 .* (atoms[i].pos - atoms[i].x0)    #calculate forces
  end
end

print(atoms)
traj = zeros(Float64,nsteps,3)
function integrate(atoms,dt)
  for j in 1:nsteps
    posx = atoms[j].pos[1]
    posy = atoms[j].pos[2]
    posz = atoms[j].pos[3]

    #print(posx,", ")
    flush(stdout)

    #traj[j][1] = posx
    #traj[j][2] = posy
    #traj[j][3] = posz

    atoms[j].pos = atoms[j].pos + atoms[j].vel * dt + 0.5 * atoms[j].acc * dt^2   #Integrate
    atoms[j].acc_old = atoms[j].acc

    forces(atoms)         #Calculate forces

    atoms[j].vel = atoms[j].vel + 0.5 * (atoms[j].acc + atoms[j].acc_old) * dt  #Update velocities
  end
end

#****************************
#Initialize some stuff
#****************************

println("")
println("***************")
println("Running dynamics...")


gaussify(atoms,25)
forces(atoms)


#****************************
#Write trajectories
#****************************
function write_output(filename,atoms)
  touch("traj.dat")
  outfile = "traj.dat"
  open(outfile, "w") do f
    for i in traj
      println(f, i)
    end
  end
  writedlm("traj.dat",traj)
end





#****************************
#Run the dynamics
#****************************

for n in 1:nsteps
  integrate(atoms,dt)

end

println("")
println("... Done!\n")
println("")
