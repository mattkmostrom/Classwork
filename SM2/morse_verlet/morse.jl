#To do:
#1.) Make a box
import SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
println("________________________\n\n")
println("Morse Potential MD v0.2\n")

global nsteps = 400       #number of integration points
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

#He(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("He",mass,eps,sigma,w,x0,pos,vel,acc,acc_old)
Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",39.95,128.326802,3.371914,1,x0,pos,vel,acc,acc_old) #define some atom parameters

#****************************
#Make the cubic lattice
#****************************

n = 2
dx = 5.
atoms = []

for i in 1:n    #Make a cubic lattice
  for j in 1:n
    for k in 1:n
      x0 = [i,j,k] .* dx
      #atom = Ar([i,j,k] .* dx,[i,j,k] .* dx,rand(3),zeros(3),zeros(3))
      atom = Ar(x0,x0,rand(3),zeros(3),zeros(3)) #Ar is a function of 5 vectors: x0,pos,vel,acc,acc_old
      push!(atoms,atom)
    end #i
  end   #j
end     #k


#****************************
#Make the velocities gaussian
#****************************

function gaussify(atoms,target_temp)
  itemp = 0
  for i in 1:length(atoms)
    atoms[i].vel[1] -= 0.5
    atoms[i].vel[2] -= 0.5
    atoms[i].vel[3] -= 0.5
    itemp += atoms[i].mass * dot(atoms[i].vel,atoms[i].vel)
  end

  itemp /= 3 * length(atoms) #- 3      #Need to comment out the '-3' if there's only 1 atom
  for i in length(atoms)
    atoms[i].vel *= sqrt(target_temp/itemp)
  end
end


function forces(atoms)      #calculate forces
  for i in 1:length(atoms)
    atoms[i].acc = -(atoms[i].w) ^ 2 .* (atoms[i].pos - atoms[i].x0)    #F/m = -dU/dr, U=k(x-x0)^2
  end
end

#print(atoms)
#println("")
function integrate(atoms,dt)
  traj = zeros(Float64,length(atoms),3)

  for j in 1:length(atoms)
    traj[j,1] = atoms[j].pos[1]
    traj[j,2] = atoms[j].pos[2]
    traj[j,3] = atoms[j].pos[3]

    println("")
    print(traj)

    atoms[j].pos = atoms[j].pos + atoms[j].vel * dt + 0.5 * atoms[j].acc * dt^2   #Integrate
    atoms[j].acc_old = atoms[j].acc

    forces(atoms)         #Calculate forces

    atoms[j].vel = atoms[j].vel + 0.5 * (atoms[j].acc + atoms[j].acc_old) * dt  #Update velocities
    return traj
  end
end

#****************************
#Initialize some stuff
#****************************

println("")
println("***************")
println("Running dynamics...")

gaussify(atoms,250)
forces(atoms)

#****************************
#Run the dynamics
#****************************
for n in 1:nsteps
  integrate(atoms,dt)
  println("Steps completed:",n,"/",nsteps)
end
print(traj)

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




println("")
println("Done!\n")
println("")
