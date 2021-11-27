#To do:
#1.) Make a box
#2.) Check that atoms don't move if force=0
#3.) Check that atoms accelerate away if force = constant
import SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
println("________________________\n\n")
println("Morse Potential MD v0.2\n")

global nsteps = 300       #number of integration points
global cell = 10
global dt = 0.001         #timestep size
global n = 2              #number of atoms to be simulated
T = 3
dx = 3.1                   #spacing between adjacent atoms in cube
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
#alpha = 1                  #well width. bigger alpha means gentler sloped well. "bond stiffness"

#He(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("He",mass,eps,sigma,w,x0,pos,vel,acc,acc_old)
#Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",39.95,128.326802,3.371914,1,x0,pos,vel,acc,acc_old) #define some atom parameters
Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",39.95,50,3,1,x0,pos,vel,acc,acc_old) #define some atom parameters


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
#Take out Center of Mass velocities
#****************************
function com(atoms)
      sumvx= 0.
      sumvy= 0.
      sumvz= 0.
      for i in 1:length(atoms)
        sumvx=atoms[i].vel[1]+sumvx
        sumvy=atoms[i].vel[2]+sumvy
        sumvz=atoms[i].vel[3]+sumvz
      end

      #calculate the center of mass velocity
      sumvx=sumvx/length(atoms)
      sumvy=sumvy/length(atoms)
      sumvz=sumvz/length(atoms)

       #subtract off the center of mass velocity
      for i in 1:length(atoms)
         atoms[i].vel[1]=atoms[i].vel[1]-sumvx
         atoms[i].vel[2]=atoms[i].vel[2]-sumvy
         atoms[i].vel[3]=atoms[i].vel[3]-sumvz
      end
end

#****************************
#Calculate forces (update potentials)
#****************************
function forces(atoms)      #calculate forces
  for i in 1:length(atoms)
    for j in 1:length(atoms)
      if atoms[i] != atoms[j]

        dx = atoms[i].pos[1] - atoms[j].pos[1]
        dy = atoms[i].pos[2] - atoms[j].pos[2]
        dz = atoms[i].pos[3] - atoms[j].pos[3]

        dx = dx - round(Int64,(dx/cell))*cell
        dy = dy - round(Int64,(dy/cell))*cell
        dx = dz - round(Int64,(dz/cell))*cell

        rx = dot(dx,dx)
        ry = dot(dy,dy)
        rz = dot(dz,dz)
        r = sqrt(rx + ry + rz)

        dx = dx / r
        dy = dy / r
        dz = dz / r

        #Lennard-Jones
        #sig_r = atoms[i].sig / r
        #U = 4 * atoms[i].eps * (((sig_r)^12)-((sig_r)^6))
        #force = 24. * atoms[i].eps / (r ^ 2) * ( 2 * (atoms[i].sig / r) ^ 12 - (atoms[i].sig / r)^6 )

        #FIXME: write if statement to choose potential type

        #Morse
        D = atoms[i].eps
        dr = r - atoms[i].sig
        expar = exp(- atoms[i].w * dr)
        U = D * (1.0 - expar) * (1.0 - expar)
        dudr = 2.0 * D * atoms[i].w * expar * (1.0 - expar)   #U = D * (1.0 - expar) * (1.0 - expar)
        force = (dudr / r) * dx

        atoms[i].acc += force / atoms[i].mass             #a = F/m
      end
    end
  end
end

function get_energy(atoms)
  total_vel = 0
  global energy = 0
  for i in 1:length(atoms)
    for j in 1:length(atoms)
      if atoms[i] != atoms[j]
        dx = atoms[i].pos - atoms[j].pos
        r = sqrt(dot(dx,dx))
        dx = dx / r

        #Lennard-Jones
        #sig_r = atoms[i].sig / r
        #U = 4 * atoms[i].eps * (((sig_r)^12)-((sig_r)^6))

        #FIXME: write if statement to choose potential type

        #Morse
        D = atoms[i].eps
        dr = r - atoms[i].sig
        alpha = 1                  #well width. bigger alpha means gentler sloped well. "bond stiffness"
        expar = exp(- alpha * dr)
        U = D * (1.0 - expar) * (1.0 - expar)

        for x in 1:3
          total_vel = sum((atoms[i].vel[x])^2)
        end

        KE = sum(atoms[i].mass * total_vel) / n^3

        energy = KE + U
      end
    end
    return energy
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
com(atoms)
forces(atoms)


#****************************
#Run the dynamics
#****************************
for q in 1:nsteps
  integrate(atoms,dt)
  println("Steps completed:",q,"/",nsteps)
  get_energy(atoms)
  println("Energy: ",energy)
  write_output(traj)
end

println("")
println("Done!\n")
println("")
