#To do:
#1.) Make a box
using SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
using Distributions
using Plots

println("________________________\n\n")
println("Morse Potential MD v0.9\n")

global kb = 1.38064852e-23
global nsteps = 700                #number of integration points
global dt = 0.05                    #timestep size
global n = 5                        #cube root of the number of atoms to be simulated
global dx = 3.0                      #spacing between adjacent atoms in cube
global origin = [0.,0.,0.]
println("Number of steps: ", nsteps)
println("Timestep size: ",dt)
println("")

#****************************
#Make the units correct
#****************************
Tstar = 0.1                 #Temperature
rhostar = 0.6               #density
massstar = 32.              #twisting your mind and smashing your dreams
eps = 128.326802            #well-depth; this is in units of eps/kb, so T will just be Tstar * eps
sig = 3.371914              #place where E = 0, NOT x corresponding to bottom of well
w = 3.                      #arbitrary potential-specific parameter; harmonic freq, morse alpha, etc #alpha is well width. bigger alpha means gentler sloped well. "bond stiffness"

T = Tstar * eps
mass = massstar * (10. / (6.022169 * 1.380662))
rho = rhostar / (sig)^3


#cell = 15.                         #Cell size
cell = (n ^ 3 / rho) ^ (1. / 3.)    #reduced units
global cutoff = cell/2              #how far atoms can see each other


@printf("Reduced Density: %f\nReduced Cell size: %f",rho,cell)


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
#Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",39.95,128.326802,3.371914,1,x0,pos,vel,acc,acc_old) #LJ UFF atom parameters
Ar(x0::Vector,pos::Vector,vel::Vector,acc::Vector,acc_old::Vector) = Element("Ar",mass,eps,sig,1,x0,pos,vel,acc,acc_old) #toy parameters


#****************************
#Make the cubic lattice
#****************************
atoms = []
for i in 1:n
  for j in 1:n
    for k in 1:n
      x0 = [i,j,k] .* dx
      atom = Ar(x0,x0,zeros(3),zeros(3),zeros(3)) #Ar is a function of 5 vectors: x0,pos,vel,acc,acc_old
      push!(atoms,atom)
    end #i
  end   #j
end     #k


#****************************
#Make the velocities gaussian
#****************************
function init_temper(atoms,target_temp)
  for i in 1:length(atoms)
      atoms[i].vel .+= rand(3)
  end

  itemp = 0
  for i in 1:size(atoms,1)
    atoms[i].vel[1] -= 0.5
    atoms[i].vel[2] -= 0.5
    atoms[i].vel[3] -= 0.5
    itemp += atoms[i].mass * dot(atoms[i].vel,atoms[i].vel)
  end

  itemp /= 3 * size(atoms,1) - 3      #Need to comment out the '-3' if there's only 1 atom
  for i in 1:size(atoms,1)
    atoms[i].vel *= sqrt(target_temp/itemp)
  end
end

#****************************
#Rescale velocities for thermostat
#****************************
function temper(atoms,target_temp)
    itemp = 0
  for i in 1:size(atoms,1)
    itemp += atoms[i].mass * dot(atoms[i].vel,atoms[i].vel)
  end

  itemp /= 3 * size(atoms,1) - 3      #Need to comment out the '-3' if there's only 1 atom
  for i in 1:size(atoms,1)
    atoms[i].vel *= sqrt(target_temp/itemp)
  end
end

#****************************
#Take out Center of Mass velocities
#****************************
function com(atoms)
      sumvx = 0.
      sumvy = 0.
      sumvz = 0.

      for i in 1:size(atoms,1)
        sumvx = atoms[i].vel[1] + sumvx
        sumvy = atoms[i].vel[2] + sumvy
        sumvz = atoms[i].vel[3] + sumvz
      end

      #calculate the center of mass velocity
      sumvx = sumvx / size(atoms,1)
      sumvy = sumvy / size(atoms,1)
      sumvz = sumvz / size(atoms,1)

       #subtract off the center of mass velocity
      for i in 1:size(atoms,1)
         atoms[i].vel[1] = atoms[i].vel[1] - sumvx
         atoms[i].vel[2] = atoms[i].vel[2] - sumvy
         atoms[i].vel[3] = atoms[i].vel[3] - sumvz
      end
end


#****************************
#Calculate forces (update potentials)
#****************************
function forces(atoms,cutoff)
  force = [0.,0.,0.]               #initialize forces, might as well make them zero floats
  for i in 1:size(atoms,1)         #zero out accelerations
    atoms[i].acc = zeros(3)
  end

  for i in 1:size(atoms,1)
    for j in 1:size(atoms,1)
      #if i < j
      if atoms[i] != atoms[j]   #should this be i < j?

        d_pos = atoms[i].pos - atoms[j].pos    #[dx,dy,dz], Float64
        dr = sqrt(dot(d_pos,d_pos))     #sqrt(dx^2 + dy^2 + dz^2), Float64

        for d in 1:3              #Periodic Boundaries
          if(atoms[i].pos[d] > (cell/2))
            atoms[i].pos[d] -= cell

          elseif(atoms[i].pos[d] <= (cell/2)) #should be negative?
            atoms[i].pos[d] += cell
          end
        end

        if(dr<cutoff)

          d_pos = d_pos/dr                #[dx/dr, dy/dr , dz/dr]

#         FIXME: write if statement to choose potential type
          #******************************
          #Lennard-Jones
          #******************************
#         dudr = 24. * atoms[i].eps / (r ^ 2) * ( 2 * (atoms[i].sig / r) ^ 12 - (atoms[i].sig / r)^6 )

          #******************************
          #Morse
          #******************************
          D = atoms[i].eps
          expar = exp(- atoms[i].w .* (dr .- atoms[i].sig))       #w = well width. bigger alpha means gentler sloped well. "bond stiffness"
          dudr = (2.0 * D * atoms[i].w) .* expar .* (1.0 - expar)   #U = D * (1.0 - expar) ^ 2

          force[1] = (-dudr) * (d_pos[1]/dr)    #(dU/dr)*(dr/dx)
          force[2] = (-dudr) * (d_pos[2]/dr)
          force[3] = (-dudr) * (d_pos[3]/dr)
          force = force .* (atoms[1].sig/atoms[1].eps)        #force in reduced units

      else(dr > cutoff)
          force = [0.,0.,0.]
        end
      end       #if i!=j

    end         #j
  end           #i
end             #forces()


function get_energy(atoms)
  global energy = 0.0

  for i in 1:size(atoms,1)
    for j in 1:size(atoms,1)
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
        alpha = 2                  #well width. bigger alpha means gentler sloped well. "bond stiffness"
        expar = exp(- alpha * dr)
        U = D * (1.0 - expar) * (1.0 - expar)

        total_vel = 0.0
        summ = 0.0

        for x in 1:3
          total_vel += (atoms[i].vel[x])^2
        end

        KE = sum(atoms[i].mass * total_vel) / n^3

        summ = KE + U
        energy += summ

      end
    end
    return energy
  end
end

#****************************
#Integrate the positions and update the velocities
#****************************
function integrate(atoms,dt)
  for j in 1:length(atoms)
  #for j in 1:size(atoms,1)

    atoms[j].pos = atoms[j].pos .+ atoms[j].vel * dt .+ (0.5 .* atoms[j].acc * (dt^2.))   #Integrate
    atoms[j].acc_old = atoms[j].acc

    forces(atoms,cutoff)         #Calculate forces

    atoms[j].vel = atoms[j].vel .+ (0.5 .* (atoms[j].acc + atoms[j].acc_old) * dt)  #Update velocities
  end

end


#****************************
#Write trajectories
#****************************
function write_pos(atoms)
  natoms = n^3
  outfile = "traj.xyz"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(atoms, 1)
      @printf(f,"%s %lf %lf %lf\n","Ar",atoms[i].pos[1],atoms[i].pos[2],atoms[i].pos[3])
    end
  end
end

function write_vel(atoms)
  natoms = n^3
  outfile = "vel.dat"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(atoms, 1)
      @printf(f,"%s %lf %lf %lf\n","Ar",atoms[i].vel[1],atoms[i].vel[2],atoms[i].vel[3])
    end
  end
end


#****************************
#Initialize some stuff
#****************************
println("")
println("***************")
println("Running dynamics...")

init_temper(atoms,T)
println("Temperature: ",T)
println("")

com(atoms)
forces(atoms,cutoff)        #give a value to old acceleration in integration loop
rm("traj.xyz")


#****************************
#Run the dynamics
#****************************
for q in 1:nsteps
  integrate(atoms,dt)
  temper(atoms,T)
  get_energy(atoms)
  write_pos(atoms)
  #write_vel(atoms)
  println("Steps completed:",q,"/",nsteps,", Energy: ",energy)
  #println("Steps completed:",q,"/",nsteps)
end

#****************************
#Write the histogram
#****************************
function write_vel_hist(atoms)
  natoms = n^3
  outfile = "vel_hist.dat"
  open(outfile, "a") do f
    for i in 1:size(atoms, 1)
      @printf(f,"%lf %lf %lf\n",atoms[i].vel[1],atoms[i].vel[2],atoms[i].vel[3])
    end
  end
end

function write_speed_hist(atoms)
  natoms = n^3
  outfile = "speed_hist.dat"
  open(outfile, "a") do f
    for i in 1:size(atoms, 1)
      @printf(f,"%lf\n%lf\n%lf\n",atoms[i].vel[1],atoms[i].vel[2],atoms[i].vel[3])
    end
  end
end



println("Writing vel hist")
write_vel_hist(atoms)
println("Writing speed hist")
write_speed_hist(atoms)


println("Collecting vel vectors")
list_size = size(atoms,1)
data = []

for n in 1:list_size
  for m in 1:3
    datum = atoms[n].vel #should be the list of velocities from final frame?
    push!(data,datum)
  end
end

println("")
println("Doing flattening")
A = collect(Iterators.flatten(data)) #single vector of velocities from final step

println("Writing histogram")
histogram(A,bins=50)

println("Saving .png")
savefig("/home/mkmostro/classes/SM2/morse_verlet/plot.png")


println("")
println("Done!\n")
println("")
println("Position trajectory is in 'traj.xyz', velocity trajectory is in 'vel.dat'.")
println("")
