using ParallelTemperingMonteCarlo
using Random

#demonstration of the new verison of the new code   

#-------------------------------------------------------#
#-----------------------MC Params-----------------------#
#-------------------------------------------------------#

Random.seed!(1234)

# number of atoms
n_atoms = 32
pressure = 101325
#pressure = 100000
#AtoBohr = 1.8897259886 #without B, unit of distance is Bohr
AtoBohr = 1.0    #with B, #When the unit of distance is still Angstrom

# temperature grid
ti = 50
tf = 3750
n_traj = 24

temp = TempGrid{n_traj}(ti,tf) 

# MC simulation details

mc_cycles = 10000 #default 20% equilibration cycles on top


mc_sample = 1  #sample every mc_sample MC cycles

#move_atom=AtomMove(n_atoms) #move strategy (here only atom moves, n_atoms per MC cycle)
displ_atom = 0.1 # Angstrom
n_adjust = 100

max_displ_atom = [0.1*sqrt(displ_atom*temp.t_grid[i]) for i in 1:n_traj]

mc_params = MCParams(mc_cycles, n_traj, n_atoms, mc_sample = mc_sample, n_adjust = n_adjust)


#-------------------------------------------------------------#
#----------------------Potential------------------------------#
#-------------------------------------------------------------#

c=[-1.81754233e-01, -2.32682279e+02,  7.49842579e+03, -1.56977364e+04, -6.15601605e+05,  3.50732411e+06]
#c=[-10.5097942564988, 989.725135614556, -101383.865938807, 3918846.12841668, -56234083.4334278, 288738837.441765]
#c=[-123.63510161951,21262.8963716972,-3239750.64086661,189367623.844691,-4304257347.72069,35314085074.72069]
pot = ELJPotentialEven{6}(c)


link="/Users/phoebemeister/Julia_MC/ParallelTemperingMonteCarlo.jl/scripts/look-up_table-he_correct.txt"
potlut=LookuptablePotential(link)
#-------------------------------------------------------------#
#------------------------Move Strategy------------------------#
#-------------------------------------------------------------#
separated_volume=false
ensemble = NPT(n_atoms,pressure*2.2937122783969076e-13/AtoBohr^3,separated_volume)
move_strat = MoveStrategy(ensemble)

#-------------------------------------------------------------#
#-----------------------Starting Config-----------------------#
#-------------------------------------------------------------#
#starting configurations
#icosahedral ground state of Ne13 (from Cambridge cluster database) in Angstrom
#use same for helium
pos_he32 =  [[ -4.3837,       -4.3837,       -4.3837],
  [-2.1918,       -2.1918,       -4.3837],
  [-2.1918,       -4.3837,       -2.1918],
  [-4.3837,       -2.1918,       -2.1918],
  [-4.3837,       -4.3837,        0.0000],
  [-2.1918,       -2.1918,        0.0000],
  [-2.1918,       -4.3837,        2.1918],
  [-4.3837,       -2.1918,        2.1918],
  [-4.3837,        0.0000,       -4.3837],
  [-2.1918,        2.1918,       -4.3837],
  [-2.1918,        0.0000,       -2.1918],
  [-4.3837,        2.1918,       -2.1918],
  [-4.3837,        0.0000,        0.0000],
  [-2.1918,        2.1918,        0.0000],
  [-2.1918,        0.0000,        2.1918],
  [-4.3837,        2.1918,        2.1918],
 [0.0000,       -4.3837,       -4.3837],
 [2.1918,       -2.1918,       -4.3837],
 [2.1918,       -4.3837,       -2.1918],
 [0.0000,       -2.1918,       -2.1918],
 [0.0000,       -4.3837,        0.0000],
 [2.1918,       -2.1918,        0.0000],
 [2.1918,       -4.3837,        2.1918],
 [0.0000,       -2.1918,        2.1918],
 [0.0000,        0.0000,       -4.3837],
 [2.1918,        2.1918,       -4.3837],
 [2.1918,        0.0000,       -2.1918],
 [0.0000,        2.1918,       -2.1918],
 [0.0000,        0.0000,        0.0000],
 [2.1918,        2.1918,       0.0000],
 [2.1918,        0.0000,        2.1918],
 [0.0000,        2.1918,        2.1918]]


pos_he32 = pos_he32 * AtoBohr

#binding sphere
box_length = 8.7674 * AtoBohr
box_height = 15.34597014 * AtoBohr
#CubicBC (only box_length) -> RhombicBC (also box_height)
bc_he32 = RhombicBC(box_length, box_height)   

length(pos_he32) == n_atoms || error("number of atoms and positions not the same - check starting config")

start_config_1 = Config(pos_he32, bc_he32)
start_config_2 = Config(pos_he32, bc_he32)
start_config = [start_config_1,start_config_2]


#----------------------------------------------------------------#
#-------------------------Run Simulation-------------------------#
#---------------------------------------x-------------------------#
mc_states, results = ptmc_run!(mc_params,temp,start_config,potlut,ensemble)

#to check code in REPL
#@profview ptmc_run!(mc_params,temp,start_config,pot,ensemble)
#@benchmark ptmc_run!(mc_params,temp,start_config,pot,ensemble)

## # plot(temp.t_grid,results.heat_cap)
# heatmap(results.ev_histogram[ ## ]) <- vary the array number to move up and down
# temp_result,cp = multihistogram_NPT(ensemble, temp, results, 10^(-6), false)
# plot(temp_result,cp)
