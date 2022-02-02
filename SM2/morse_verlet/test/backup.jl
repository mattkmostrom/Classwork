using SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
using Distributions
using Plots
using BenchmarkTools
using Profile

println("________________________\n\n")
println("Morse Potential MD v1.0\n")
println("")

global type = "LJ"            #Potential type, "LJ" or "morse" or "harmonic"
global kb = 1.38064852e-23
global n = 2                        #cube root of the number of atoms to be simulated

global dx = 1.12246204831           #spacing between adjacent atoms in cube

#****************************
#Make the units correct
#****************************
Tstar = 0.8511              #Temperature
rhostar = 0.776             #density
massstar = 32.              #twisting your mind and smashing your dreams
eps = 1.            #well-depth; this is in units of eps/kb, so T will just be Tstar * eps
sig = 1.              #place where E = 0, NOT x corresponding to bottom of well
w = 1.                      #arbitrary potential-specific parameter; harmonic freq, morse alpha, etc #alpha is well width. bigger alpha means gentler sloped well. "bond stiffness"

T = Tstar * eps
#mass = massstar * (10. / (6.022169 * 1.380662))
mass = 1.
rho = rhostar / (sig)^3

time_star = mass * (sig * sig / eps)
time = 100 * time_star
#global dt = 0.005 * time_star
global dt = 0.01
#global nsteps = round(time / dt)    #number of integration points
nsteps = 2
println("Timestep size: ",dt)       #timestep size
println("Number of steps: ", nsteps)

#cell = (n ^ 3 / rho) ^ (1. / 3.)    #Cell size reduced units
cell = 100000.
global cutoff = cell/2              #how far atoms can see each other

@printf("Cutoff Distance: %lf\n",cutoff)
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
  #for j in 1:n
   # for k in 1:n
      #x0 = [i,j,k] .* dx
      x0 = [i,0,0] .* dx
      atom = Ar(x0,x0,zeros(3),zeros(3),zeros(3)) #Ar is a function of 5 vectors: x0,pos,vel,acc,acc_old
      push!(atoms,atom)
    #end #i
  #end   #j
end     #k

natoms = size(atoms,1)

atoms[1].vel = [1.,0.,0.]
atoms[2].vel = [-1.,0.,0.]
#atoms[1].pos = [-1.,0.,0.]
#atoms[2].pos = [0.,0.,0.]


function write_init_pos(atoms)
  natoms = n
  outfile = "init_pos.xyz"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(atoms, 1)
      @printf(f,"%s %lf %lf %lf\n","Ar",atoms[i].x0[1],atoms[i].x0[2],atoms[i].x0[3])
    end
  end
end

touch("init_pos.xyz")
rm("init_pos.xyz")
write_init_pos(atoms)

function write_init_vel(atoms)
  natoms = n
  outfile = "init_vel.xyz"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(atoms, 1)
      @printf(f,"%s %lf %lf %lf\n","Ar",atoms[i].vel[1],atoms[i].vel[2],atoms[i].vel[3])
    end
  end
end

function write_init_mass(atoms)
  natoms = n
  outfile = "init_mass.xyz"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(atoms, 1)
      @printf(f,"%s %lf\n","Ar",atoms[i].mass)
    end
  end
end


touch("init_vel.xyz")
rm("init_vel.xyz")
write_init_vel(atoms)
touch("init_mass.xyz")
rm("init_mass.xyz")
write_init_mass(atoms)

#****************************
#Initialize random velocities
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
#Make random distribution of velocities gaussian
#****************************
function gaussify(atoms)
  
  #https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Basic_form
  for i in 1:size(atoms,1)
    U1 = rand()
    U2 = rand()
    R = sqrt(-2 * log(U1))
    theta = 2*pi*U2
    X1 = R * cos(theta) #both X1 and X2 are similarly gaussian distributed 
    X2 = R * sin(theta) #and either can be used to substitute atoms[i].vel[n]
    atoms[i].vel[1] = X1

    U1 = rand()
    U2 = rand()
    R = sqrt(-2 * log(U1))
    theta = 2*pi*U2
    Y1 = R * cos(theta)
    Y2 = R * sin(theta)
    atoms[i].vel[2] = Y1

    U1 = rand()
    U2 = rand()
    R = sqrt(-2 * log(U1))
    theta = 2*pi*U2
    Z1 = R * cos(theta)
    Z2 = R * sin(theta)
    atoms[i].vel[3] = Z1
  end
end
#gaussify(atoms)

#****************************
#Turn velocity PDF into a speed PDF
#****************************
speeds = []
function speed(atoms)
  for i in 1:size(atoms,1)
    for n in 1:3
      speed = sqrt(atoms[i].vel[n] * atoms[i].vel[n])
      push!(speeds,speed)
    end
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
function anticom(atoms)
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
      if atoms[i] != atoms[j]   #should this be i < j?

        d_pos = atoms[i].pos .- atoms[j].pos    #[dx,dy,dz], Float64

        for d in 1:3              #Periodic Boundaries
          if(atoms[i].pos[d] > (cell/2))
            atoms[i].pos[d] -= cell

          elseif(atoms[i].pos[d] <= (-cell/2)) #should be negative?
            atoms[i].pos[d] += cell
          end
        end

        dr = sqrt(dot(d_pos,d_pos))     #sqrt(dx^2 + dy^2 + dz^2), Float64

        #println("dr: ",dr)

        if(dr<cutoff)

          #d_pos = d_pos/dr                #[dx/dr, dy/dr , dz/dr]
        
          if type == "LJ"
            #******************************
            #Lennard-Jones
            #******************************
            #dudr = 24. * atoms[i].eps / (dr ^ 2) * ( 2 * (atoms[i].sig / dr) ^ 12 - (atoms[i].sig / dr) ^ 6 )
            dudr = 24. * atoms[i].eps * (((atoms[i].sig^6 ) / (dr^7)) - (2. * ((atoms[i].sig^12)/(dr^13))))

          elseif type == "morse"
            #******************************
            #Morse
            #******************************
            D = atoms[i].eps
            expar = exp(- atoms[i].w .* (dr .- atoms[i].sig))       #w = well width. bigger alpha means gentler sloped well. "bond stiffness"
            dudr = (2.0 * D * atoms[i].w) .* expar .* (1.0 - expar)   #U = D * (1.0 - expar) ^ 2
            
          elseif type == "harmonic"
            #U = 0.5 * atoms[i].w * dr * dr
            dudr = atoms[i].w * dr        
          end

          #println("dudr: ",dudr)

          #force = (-dudr .* d_pos[1]) / dr
          force[1] = (-dudr) * (d_pos[1]/dr)    #(dU/dr)*(dr/dx)
          force[2] = (-dudr) * (d_pos[2]/dr)
          force[3] = (-dudr) * (d_pos[3]/dr)
          #force = force .* (atoms[1].sig/atoms[1].eps)        #force in reduced units

          for n in 1:3
            atoms[i].acc[n] = (dr * force[n]) / atoms[i].mass
          end

          #atoms[i].acc = (dr .* force) / atoms[i].mass

      else(dr > cutoff)
          force = [0.,0.,0.]
        end
      end       #if i!=j
    end         #j
  end           #i
  

  println("START FORCES LOOP DUMP ")
  d_pos = atoms[1].pos .- atoms[2].pos
  dr = sqrt(dot(d_pos,d_pos))
  println("")
  println("Position of Atom ",1," x:",atoms[1].pos[1],", y: ",atoms[1].pos[2],", z: ",atoms[1].pos[3])
  println("Position of Atom ",2," x:",atoms[2].pos[1],", y: ",atoms[2].pos[2],", z: ",atoms[2].pos[3])
  println("d_pos: ",d_pos)
  println("Distance between ",1," and ",2,": ",dr)
  println("dr^7: ",dr^7)
  println("dr^13: ",dr^13)
  
  println("")
  println("Force on Atom ",1," x:",force[1],", y: ",force[2],", z: ",force[3])
#=  
  println("")
  println("Acceleration on Atom ",1," x:",atoms[1].acc[1],", y: ",atoms[1].acc[2],", z: ",atoms[1].acc[3])
  println("Acceleration on Atom ",2," x:",atoms[2].acc[1],", y: ",atoms[2].acc[2],", z: ",atoms[2].acc[3])
  println("")
  println("")
  println("END FORCES LOOP DUMP ")
  return force
 =#
end             #forces()


function get_energy(atoms)
  global energy = 0.0
  global KE = 0.0
  global PE = 0.0
  local total_vel = 0.0
  summ = 0.0

  for i in 1:size(atoms,1)
    for j in 1:size(atoms,1)
      if atoms[i] != atoms[j]
        #***********************************************
        #Potential Energy
        #***********************************************
        dx = atoms[i].pos - atoms[j].pos
        dr = sqrt(dot(dx,dx))

        if type == "LJ"
          #Lennard-Jones
          sig_r = atoms[i].sig / dr
          U = 4 * atoms[i].eps * ( (sig_r^12) - (sig_r^6) )

        elseif type == "morse"
          #Morse
          D = atoms[i].eps
          dr = dr - atoms[i].sig
          #alpha = 2                  #well width. bigger alpha means gentler sloped well. "bond stiffness"
          expar = exp(- atoms[i].w * dr)
          U = D * (1.0 - expar) * (1.0 - expar)

        elseif type == "harmonic"
          #Harmonic
          U = 0.5 * atoms[i].w * dr * dr
          
        end #type if

        #prefactor = pi * size(atoms,1) * rhostar * cutoff^(-3)
        #lrc = ((8/9) - (8/3)) * prefactor

        PE = U #+ lrc
        #summ += PE
        total_vel += sqrt(dot(atoms[i].vel,atoms[i].vel))
        println("total velocity: ",total_vel)
      end
    end
    #***********************************************
    #Kinetic Energy
    #***********************************************
    
    
    KE = 0.5 * atoms[i].mass * total_vel ^ 2
    println("PE: ",PE)
    println("KE: ",KE)
    println("")
    energy = KE + PE
  end
end

#****************************
#Integrate the positions and update the velocities
#****************************
function integrate(atoms,dt)
  for j in 1:size(atoms,1)
  #for j in 1:size(atoms,1)
    atoms[j].pos = atoms[j].pos .+ atoms[j].vel * dt .+ (0.5 .* atoms[j].acc * (dt^2.))   #Integrate
    atoms[j].acc_old = atoms[j].acc

    forces(atoms,cutoff)         #Calculate forces

    atoms[j].vel = atoms[j].vel .+ (0.5 .* (atoms[j].acc + atoms[j].acc_old) * dt)  #Update velocities  
    
    write_pos(atoms)
    write_vel(atoms)
    get_energy(atoms)
    
    local natoms = size(atoms,1)
    outfile = "energy.dat"
    open(outfile, "a") do f
      @printf(f,"%lf\n",energy)
    end
    
    outfile = "potential_energy.dat"
    open(outfile, "a") do f
      @printf(f,"%lf\n",PE)
    end

    outfile = "kinetic_energy.dat"
    open(outfile, "a") do f
      @printf(f,"%lf\n",KE)
    end
  #=
    println("START INTEGRATION LOOP DUMP ",j)
    d_pos = atoms[1].pos .- atoms[2].pos
    dr = sqrt(dot(d_pos,d_pos))
    println("")
    println("Position of Atom ",1," x:",atoms[1].pos[1],", y: ",atoms[1].pos[2],", z: ",atoms[1].pos[3])
    println("Position of Atom ",2," x:",atoms[2].pos[1],", y: ",atoms[2].pos[2],", z: ",atoms[2].pos[3])
    println(d_pos)
    println("Distance between ",1," and ",2,": ",dr)
    println("dr^7: ",dr^7)
    println("dr^13: ",dr^13)
    
    println("")
    println("Acceleration on Atom ",1," x:",atoms[1].acc[1],", y: ",atoms[1].acc[2],", z: ",atoms[1].acc[3])
    println("Acceleration on Atom ",2," x:",atoms[2].acc[1],", y: ",atoms[2].acc[2],", z: ",atoms[2].acc[3])
    println("")
    println("")
    println("END INTEGRATION LOOP DUMP ",j,"\n\n\n\n\n\n")
  =#
  
  end
end


#****************************
#Write useful outputs to files
#****************************
function write_pos(atoms)
  natoms = size(atoms,1)
  outfile = "traj.xyz"
  open(outfile, "a") do f
    @printf(f,"%i\n\n",natoms)
    for i in 1:size(atoms, 1)
      @printf(f,"%s %lf %lf %lf\n","Ar",atoms[i].pos[1],atoms[i].pos[2],atoms[i].pos[3])
    end
  end
end

function write_vel(atoms)
  natoms = size(atoms,1)
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


#init_temper(atoms,T)
println("\nTemperature: ",T)
println("")

anticom(atoms)
forces(atoms,cutoff)        #give a value to old acceleration in integration loop

get_energy(atoms)
println("Total inital energy is: ",energy)

println("\nInital Accelerations: ")
for i in 1:size(atoms,1)
  println("Atom ",i,": ",atoms[i].acc)
end

println("")

touch("traj.xyz")
touch("vel.dat")
touch("energy.dat")
touch("avg_energy.dat")
touch("potential_energy.dat")
touch("kinetic_energy.dat")

rm("traj.xyz")
rm("vel.dat")
rm("energy.dat")
rm("avg_energy.dat")
rm("potential_energy.dat")
rm("kinetic_energy.dat")




#****************************
#Run the dynamics
#****************************
energies = []
function run_dynamics(atoms,dt,T,n,nsteps,natoms)
  e_sum = 0.0
  println("#****************************")
  println("STARTING PRODUCTION RUN")
  println("#****************************")
  for q in 1:nsteps 
    
    integrate(atoms,dt)


    #println("\nAccelerations: ")
    #for i in 1:size(atoms,1)
    #  println("Atom ",i,": ",atoms[i].acc)
    #end

    #if q !=0 && q%10==0 && q < (nsteps/2)
    #  temper(atoms,T)
    #end
    
    #get_energy(atoms)
    
    

    println("Steps completed:",q,"/",nsteps,", Energy: ",energy)   
  end

end
run_dynamics(atoms,dt,T,n,nsteps,natoms)

println("\nFinal Accelerations: ")
for i in 1:natoms
  println("Atom ",i,": ",atoms[i].acc)
end

#****************************
#Write the histogram
#****************************
function write_vel_hist(atoms)
  natoms = size(atoms,1)
  outfile = "vel_hist.dat"
  open(outfile, "a") do f
    for i in 1:size(atoms, 1)
      @printf(f,"%lf %lf %lf\n",atoms[i].vel[1],atoms[i].vel[2],atoms[i].vel[3])
    end
  end
end


println("")
println("Writing vel hist")
write_vel_hist(atoms)
println("Writing speed hist")
speed(atoms)


println("Collecting vel vectors")
list_size = size(atoms,1)
data = []

for n in 1:list_size
  for m in 1:3
    datum = atoms[n].vel #should be the list of velocities from final frame?
    push!(data,datum)
  end
end

println("Doing flattening")
A = collect(Iterators.flatten(data)) #single vector of velocities from final step
B = collect(Iterators.flatten(speeds))

println("Writing histogram")
#histogram(A,bins=50)

println("Saving velocity plot")
#savefig("vel_plot.png")

println("Saving speed plot")
#histogram(B,bins=100)
#savefig("speed_plot.png")

println("")
println("Done!\n")
println("")
println("Position trajectory is in 'traj.xyz', velocity trajectory is in 'vel.dat'.")
println("")
