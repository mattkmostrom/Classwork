using SpecialFunctions
using LinearAlgebra
using Printf
using DelimitedFiles
using Distributions
using Plots
using Random

println("________________________\n\n")
println("Morse Potential MD v1.1\n")
println("")

global type = "LJ"           #Potential type: "LJ" or "morse" or "harmonic"
global n = 4                 #cube root of the number of atoms to be simulated

#****************************
#Set some constants and make the units correct
#****************************
Tstar = 1.64              #Temperature
rhostar = 0.9             #density
massstar = 39.9             #twisting your mind and smashing your dreams

eps = 119.8                 #well-depth; in units of eps/kb, so T will just be Tstar * eps
sig = 3.405                 #place where E = 0, NOT x corresponding to bottom of well
w = 1.                      #arbitrary potential-specific parameter: harmonic freq, morse alpha, etc. alpha is well width; bigger alpha means gentler sloped well. "bond stiffness"

T = Tstar * eps
mass = massstar * (10. / (6.022169 / 1.380662)) #KPsa
rho = rhostar / (sig^3)

time_star = mass * (sig * sig / eps)      #this block is just making the 
time = 100 * time_star                    #units what they should be in the paper
global dt = 0.5 * 0.005 * time_star
global nsteps = 3 * round(time / dt)    #number of integration points
global natoms = 32
global cell = (natoms / rho) ^ (1. / 3.)    #Cell size reduced units
global dx = cell/n           #spacing between adjacent atoms in cube

#global dt = 0.005
#nsteps = 100000
println("Timestep size: ",dt)       #timestep size
println("Number of steps: ", nsteps)


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
for i in 1:n/2
  for j in 1:n
    for k in 1:n
      x0 = ([i,j,k] .* dx)  #define some coordinate for the atom in question to sit at
      atom = Ar(x0,x0,zeros(3),zeros(3),zeros(3)) #Ar is a function of 5 vectors: x0,pos,vel,acc,acc_old
      push!(atoms,atom)
    end #i
  end   #j
end     #k


global cutoff = cell/2              #how far atoms can see each other
global cutoff2 = cutoff * cutoff

println("Number of atoms: ",natoms)
@printf("Cutoff Distance: %lf\n",cutoff)
@printf("Reduced Density: %f\nReduced Cell size: %f",rho,cell)
println("\nReduced mass: ",mass)

function write_init_pos(atoms)
  natoms = size(atoms,1)
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


#****************************
#Initialize random velocities
#****************************
function init_temper(atoms,target_temp)
  Random.seed!(1234)
  for i in 1:size(atoms,1)
      atoms[i].vel .+= rand(3)  #slap some random numbers in there
  end

  itemp = 0
  for i in 1:size(atoms,1)
    atoms[i].vel[1] -= 0.5    #change from 0-1 to -0.5-0.5
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
      speed = abs(atoms[i].vel[n]) 
      push!(speeds,speed)
    end
  end
end

#****************************
#Rescale velocities for thermostat
#****************************
function temper(atoms,target_temp)
  itemp = 0.0
  
  for i in 1:size(atoms,1)  #get total kinetic energy
    itemp += atoms[i].mass * dot(atoms[i].vel,atoms[i].vel)
  end
                                      #divide by DOF
  itemp /= 3 * size(atoms,1) - 3      #Need to comment out the '-3' if there's only 1 atom
  
  for i in 1:size(atoms,1)            #rescale temps
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
        sumvx += atoms[i].vel[1]
        sumvy += atoms[i].vel[2]
        sumvz += atoms[i].vel[3]
      end

      #calculate the center of mass velocity
      sumvx /= size(atoms,1)
      sumvy /= size(atoms,1)
      sumvz /= size(atoms,1)

       #subtract off the center of mass velocity
      for i in 1:size(atoms,1)
         atoms[i].vel[1] -= sumvx
         atoms[i].vel[2] -= sumvy
         atoms[i].vel[3] -= sumvz
      end
end


#****************************
#Calculate forces
#****************************
function forces(atoms,cutoff)

  force = [0.,0.,0.]               #initialize forces, might as well make them zero floats
  for i in 1:size(atoms,1)         #zero out accelerations
    atoms[i].acc = zeros(3)
  end

  for i in 1:size(atoms,1)-1  #no double-counting, no self-counting
    for j in i+1:size(atoms,1)
      
        d_pos = atoms[i].pos .- atoms[j].pos    #[dx,dy,dz], Float64

        for n in 1:3          #Calculate forces of nearest image
          if d_pos[n] > cutoff
            d_pos[n] -= cell
          elseif d_pos[n] <= -cutoff
            d_pos[n] += cell
          end
        end      

        dr2 = dot(d_pos,d_pos)
        

      if(dr2<cutoff2)
        dr = sqrt(dr2)     #sqrt(dx^2 + dy^2 + dz^2), Float64
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
          expar = exp(- atoms[i].w * (dr - atoms[i].sig))       #w = well width. bigger alpha means gentler sloped well. "bond stiffness"
          dudr = (2.0 * D * atoms[i].w) * expar * (1.0 - expar)   #U = D * (1.0 - expar) ^ 2

        elseif type == "harmonic"
          dudr = atoms[i].w * dr #U = 0.5 * atoms[i].w * dr * dr
        end

       
        force[1] = (-dudr) * (d_pos[1]/dr)    #(dU/dr)*(dr/dx)
        force[2] = (-dudr) * (d_pos[2]/dr)
        force[3] = (-dudr) * (d_pos[3]/dr)
       
        for n in 1:3
          atoms[i].acc[n] += force[n] / atoms[i].mass
          atoms[j].acc[n] -= force[n] / atoms[j].mass
        end

      else(dr2 >= cutoff2)
        force = [0.,0.,0.]
      end
    end         #j
  end           #i
end             #forces()


function get_energy(atoms)
  global energy = 0.0
  global KE = 0.0
  global PE = 0.0
  local total_vel = 0.0
  summ = 0.0

  for i in 1:size(atoms,1)
    total_vel += dot(atoms[i].vel,atoms[i].vel) #(Vx^2 + Vy^2 + Vz^2) = V^2
  end
  
  KE = 0.5 * atoms[1].mass * total_vel

  for i in 1:size(atoms,1)-1
    for j in i+1:size(atoms,1)
      dx = atoms[i].pos - atoms[j].pos
      
      for n in 1:3        #making sure energy is calculated from 
        if dx[n] > cutoff #nearest particle image
          dx[n] -= cell
        elseif dx[n] <= -cutoff
          dx[n] += cell
        end
      end
      
      dr2 = dot(dx,dx)

      if dr2 <= cutoff2
        dr = sqrt(dr2)
        if type == "LJ"
          #Lennard-Jones
          sig_r = atoms[i].sig / dr
          U = 4. * atoms[i].eps * ( (sig_r^12) - (sig_r^6) )    #Got ecut from UMS, idk
          ecut = 4 * atoms[i].eps * (((atoms[i].sig/cutoff)^12)-((atoms[i].sig/cutoff)^6))

        elseif type == "morse"
          #Morse
          D = atoms[i].eps
          dr = dr - atoms[i].sig
          #alpha = 2                  #well width. bigger alpha means gentler sloped well. "bond stiffness"
          expar = exp(- atoms[i].w * dr)
          U = D * (1.0 - expar) * (1.0 - expar)
          expar = exp(- atoms[i].w * cutoff)
          ecut = D * (1.0 - expar) * (1.0 - expar)

        elseif type == "harmonic"
          #Harmonic
          U = 0.5 * atoms[i].w * dr * dr
          ecut = 0.5 * atoms[i].w * cutoff * cutoff

        end #type if
        
        PE += U #- ecut              

      elseif dr2 > cutoff2
        U = 0.
      end #if cutoff
    end #j
  end #i
  energy = KE + PE
end   #energy()

#****************************
#Integrate the positions and update the velocities
#****************************
function integrate(atoms,dt,q)
  for j in 1:size(atoms,1)          #https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
    atoms[j].pos = atoms[j].pos .+ atoms[j].vel * dt .+ (0.5 .* atoms[j].acc * (dt^2.))   #update positions
    atoms[j].acc_old = atoms[j].acc
  end

  forces(atoms,cutoff)         #Calculate forces

  for j in 1:size(atoms,1)
    atoms[j].vel = atoms[j].vel .+ (0.5 .* (atoms[j].acc + atoms[j].acc_old) * dt)  #Update velocities
  end

  for i in 1:size(atoms,1)                 #Putting particles back in the box
    for d in 1:3                           
      if(atoms[i].pos[d] > (cell/2))
        atoms[i].pos[d] -= cell
      elseif(atoms[i].pos[d] <= (-cell/2))
        atoms[i].pos[d] += cell
      end
    end
  end

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

init_temper(atoms,T)            #Give the atoms some random initial velocities
gaussify(atoms)                 #make those random velocities a gaussian PDF
println("\nTemperature: ",T)
println("")

anticom(atoms)
forces(atoms,cutoff)        #give a value to old acceleration in integration loop


println("")         
touch("traj.xyz")     #so you can start with empty files every time you run md.jl
touch("vel.dat")
touch("energy.dat")
touch("avg_energy.dat")
touch("potential_energy.dat")
touch("kinetic_energy.dat")
touch("itemp.dat")
rm("traj.xyz")
rm("vel.dat")
rm("energy.dat")
rm("avg_energy.dat")
rm("potential_energy.dat")
rm("kinetic_energy.dat")
rm("itemp.dat")

#****************************
#Run the dynamics
#****************************
energies = []
function run_dynamics(atoms,dt,T,n,nsteps,natoms)
  itemp = 0.0
  e_sum = 0.0
  av_E = 0.0
  
  for q in 1:nsteps
    integrate(atoms,dt,q)

    if q !=0 && q%10==0 && q < (nsteps/4) #thermostat the first few steps
      temper(atoms,T)
    end

    if q !=0 && q >= (nsteps/4)           #collecting energy averaging after
      av_E += PE                          #thermostat turns off
    end                                   

    if q == nsteps                        
      av_E = av_E / (nsteps-((nsteps/4)))   #divide only by steps that we used to add them up 
      av_E /= eps
      println("Av. red. energy: ",av_E)

      sigorcut= sig/cutoff
      lrc = 8. /3. * eps * pi * natoms
      lrc *= rhostar * (1. /3. * sigorcut^9 - sigorcut^3)

      #prefactor = pi * size(atoms,1) * rhostar
      #lrc = (8/9)*prefactor*(cutoff/sig)^(-9) - (8/3)*prefactor*(cutoff/sig)^(-3) #from some paper brian linked        
      
      println("lrc: ",lrc)
      println("lrc per particle: ",lrc/natoms)
      lrc /= eps
      println("reduced lrc: ",lrc)
      println("reduced lrc per particle: ",lrc/natoms)
      av_E = av_E + lrc
      println("\n\n")
    end                                   

    for i in 1:size(atoms,1)                                    #this is the same as what's in temper()
      itemp += atoms[i].mass * dot(atoms[i].vel,atoms[i].vel)   #it's just here because julia is difficult 
    end                                                         #when it comes to variable scoping
    itemp /= 3 * size(atoms,1) - 3  #Need to comment out the '-3' if there's only 1 atom

    natoms = size(atoms,1)
    outfile = "itemp.dat"
    open(outfile, "a") do f
      for i in 1:size(atoms, 1)           #records the temperature at each step,
        @printf(f,"%f\n",itemp)           #but step number isn't explicitly written
      end                                 #so you might need to write it out yourself
    end

    println("Steps completed:",q,"/",nsteps,", Total Energy: ",energy,", Temperature: ",itemp)
  end
  println("\nAverage Energy: ",av_E)
  println("Average Energy per particle: ",av_E / natoms)
end
run_dynamics(atoms,dt,T,n,nsteps,natoms)


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

println("Writing velocity histogram")
histogram(A,bins=50)
println("Saving velocity plot")
savefig("vel_plot.png")

println("Writing velocity histogram")
histogram(B,bins=100)
println("Saving speed plot")
savefig("speed_plot.png")

println("")
println("Done!\n")
println("")

@printf("Cutoff Distance: %lf\n",cutoff)
@printf("Reduced Density: %f\nReduced Cell size: %f",rho,cell)
println("\nReduced mass: ",mass)

println("\nTrajectory: vmd traj.xyz \nEnergy: xmgrace *energy.dat")
println("")
