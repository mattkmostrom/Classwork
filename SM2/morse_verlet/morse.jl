#To do:
#1.) Make a box
import SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
println("________________________\n\n")
println("Morse Potential MD v0.2\n")

global nsteps = 2000       #number of integration points
global dt = 0.0001         #timestep size
global n = 3                     #number of atoms to be simulated
T = 10
dx = 4.                   #spacing between adjacent atoms in cube
println("Number of steps: ", nsteps)
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
#Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",39.95,10,3.371914,1,x0,pos,vel,acc,acc_old) #define some atom parameters


#****************************
#Make the cubic lattice
#****************************
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

  itemp /= 3 * length(atoms) - 3      #Need to comment out the '-3' if there's only 1 atom
  for i in length(atoms)
    atoms[i].vel *= sqrt(target_temp/itemp)
  end
end

#****************************
#Calculate forces (update potentials)
#****************************
function forces(atoms)      #calculate forces
  for i in 1:length(atoms)
    for j in 1:length(atoms)
      if atoms[i] != atoms[j]
        dx = atoms[i].pos - atoms[j].pos
        r = sqrt(dot(dx,dx))
        #println("dx: ",dx)
        dx = dx / r
        #println("r: ",r)
        #println("dx2: ",dx)

        #Lennard-Jones
        #force = 24. * atoms[i].eps / (r ^ 2) * ( 2 * (atoms[i].sig / r) ^ 12 - (atoms[i].sig / r)^6 )

        #FIXME: write if statement to choose potential type

        #Morse
        D = atoms[i].eps
        dr = r - atoms[i].sig
        alpha = 10                  #well width. bigger alpha means gentler sloped well. "bond stiffness"
        expar = exp(- alpha * dr)
        pot = D * (1.0 - expar) * (1.0 - expar)
        dudr = 2.0 * D * alpha * expar * (1.0 - expar)
        force = (dudr / r) * dx
        #println("acceleration: ",atoms[i].acc)
        #println("force: ",force)
        #println("")
        atoms[i].acc += force .* dx / atoms[i].mass
      end
    end
  end
end

#print(atoms)
#println("")
function integrate(atoms,dt)
  global traj = zeros(Float64,length(atoms),3)
  for j in 1:length(atoms)
    traj[j,1] = atoms[j].pos[1]
    traj[j,2] = atoms[j].pos[2]
    traj[j,3] = atoms[j].pos[3]
    #println("")
    #println("")
    #print(traj)

    atoms[j].pos = atoms[j].pos + atoms[j].vel * dt + 0.5 * atoms[j].acc * dt^2   #Integrate
    atoms[j].acc_old = atoms[j].acc

    forces(atoms)         #Calculate forces

    atoms[j].vel = atoms[j].vel + 0.5 * (atoms[j].acc + atoms[j].acc_old) * dt  #Update velocities

  end
  return traj
end


#****************************
#Write trajectories
#****************************
function write_output(traj)
  float_natoms = n^3
  natoms = Int.(float_natoms)
  outfile = "traj.xyz"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(traj, 1)
      @printf(f,"%s %lf %lf %lf\n","Ar",traj[i,1],traj[i,2],traj[i,3])
    end
  end
end

#****************************
#Initialize some stuff
#****************************
println("")
println("***************")
println("Running dynamics...")
gaussify(atoms,T)
forces(atoms)

#****************************
#Run the dynamics
#****************************
for q in 1:nsteps
  integrate(atoms,dt)
  println("Steps completed:",q,"/",nsteps)
  write_output(traj)
end

println("")
println("Done!\n")
println("")
