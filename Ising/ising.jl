using Random
using ArgParse
using DataStructures
using DelimitedFiles

# Note that there is already an implementation of Wolf algorithm
# at <https://github.com/cossio/IsingModels.jl>

"""
Parse command-line arguments.
See <https://carlobaldassi.github.io/ArgParse.jl/latest/>
"""
function parse_commandline()
	s = ArgParseSettings()

	@add_arg_table s begin
		"--length", "-l"
			help = "Linear number of spins"
			arg_type = Int
			required = true
		"--beta", "-b"
			help = "Inverse temperature"
			arg_type = Float64
			required = true
		"-J"
			help = "Coupling constant"
			arg_type = Float64
			default = 1.
		"-H"
			help = "Magnetic field"
			arg_type = Float64
			default = 0.
		"--method"
			help = "Method ('markov', 'heatbath' or 'cluster')"
			arg_type = String
			default = "markov"
		"--nit" 
			help = "Number of iterations"
			arg_type = Int
			default = 100
		"--nth" 
			help = "Number of iterations of thermalization"
			arg_type = Int
			default = 0
		"--skip" 
			help = "Compute observables every given number of iterations"
			arg_type = Int
			default = 10
		"--nsim" 
			help = "Number of simulations"
			arg_type = Int
			default = 1
		"--nmoms" 
			help = "Number of moments to compute"
			arg_type = Int
			default = 4
		"--determ"
			help = "Deterministic initial condition (spins up)"
			action = :store_true
		"--output"
			help = "Name of output file"
			arg_type = String
			default = "ising"
		"--spins"
			help = "Export the raw spins at the end of each simulation"
			action = :store_true
		"--verbose", "-v"
			help = "Verbose mode"
			action = :store_true
	end

	return parse_args(s)
end

"State of the system at a given time."
mutable struct State
    spins::Array{Int8, 2}  # Spins
	magnetization::Int  # Magnetization
	energy::Float64  # Energy
end

"Moments of magnetization and energy."
mutable struct Observables
	n::Int  # Number of configs that are summed
	moms_mag::Array{Float64}  # Moments of magnetization
	moms_energy::Array{Float64}  # Moments of energy
end

"Divide moments by number of configs."
function scale!(o::Observables)
	o.moms_mag ./= o.n
	o.moms_energy ./= o.n
	o.n = 1
end

"Add o2 to o1."
function add!(o1::Observables, o2::Observables)
	o1.n += o2.n
	o1.moms_mag .+= o2.moms_mag
	o1.moms_energy .+= o2.moms_energy
end

"Data for simulation. Note that `state` and `obs` are mutable."
struct Simu
	t::Int  # Time
	L::Int  # Length of the system
	beta::Float64  # Inverse temperature
	J::Float64  # Coupling constant
	h::Float64  # Magnetic field
	state::State  # State of the system
	obs::Observables  # Observables
end

"Initialize simulation from command-line arguments."
function Simu(args::Dict)
	L = args["length"]
	s = Simu(0, L, args["beta"], args["J"], args["H"],
			 State(ones(L, L), 0, 0.),
			 Observables(0, zeros(args["nmoms"]), zeros(args["nmoms"])))

	if !args["determ"]
		rand!(s.state.spins, [-1, 1])
	end

	compute_magnetization_and_energy!(s)
	return s
end

"Compute the energy of the system."
function compute_magnetization_and_energy!(s::Simu)
	s.state.magnetization = sum(s.state.spins)
	E = -s.h * s.state.magnetization
	E -= s.J * sum(s.state.spins .* circshift(s.state.spins, (1, 0)))
	E -= s.J * sum(s.state.spins .* circshift(s.state.spins, (0, 1)))
	s.state.energy = E
end

"Return neighbors of a given site (periodic boundary conditions)."
function neighbors(k, L)
	[
	 CartesianIndex(mod1(k[1]-1, L), k[2]),
	 CartesianIndex(mod1(k[1]+1, L), k[2]),
	 CartesianIndex(k[1], mod1(k[2]+1, L)),
	 CartesianIndex(k[1], mod1(k[2]-1, L))
	]
end

"Local Metropolis algorithm (see Krauth, algorithm 5.7)."
function update_markov!(s::Simu)
	k = rand(CartesianIndices(s.state.spins))
	snbr = sum(s.state.spins[j] for j in neighbors(k, s.L))
	field = s.h + s.J * snbr
	dE = 2. * s.state.spins[k] * field
	u = exp(-s.beta * dE)
	if rand() < u
		s.state.spins[k] *= -1
		s.state.energy += dE
	end
end

"Heat bath algorithm (see Krauth, algorithm 5.8)."
function update_heatbath!(s::Simu)
	k = rand(CartesianIndices(s.state.spins))
	snbr = sum(s.state.spins[j] for j in neighbors(k, s.L))
	field = s.h + s.J * snbr
	p = 1. / (1. + exp(-2. * s.beta * field))
	s_old = s.state.spins[k]
	s_new = 2 * (rand() < p) - 1
	s.state.spins[k] = s_new
	s.state.energy += field * (s_old - s_new)
end

"Cluster algorithm (see Krauth, algorithm 5.9)."
function update_cluster!(s::Simu)
	p = -expm1(-2. * s.beta)  # Probability to add site into cluster

	C = Deque{CartesianIndex{2}}()
	P = Deque{CartesianIndex{2}}()  # Or Stack
	visited = zeros(Bool, size(s.state.spins))

	k = rand(CartesianIndices(s.state.spins))
	push!(C, k)
	push!(P, k)
	visited[k] = true

	while !isempty(P)
		i = pop!(P)
		for j in neighbors(i, s.L)
			if (!visited[j] && s.state.spins[i] == s.state.spins[j]
				&& rand() < p)
				push!(C, j)
				push!(P, j)
				visited[j] = true
			end
		end
	end

	# If there is a magnetic field, we need to enforce a Metropolis rule
	dEin = 2 * s.h * length(C) * s.state.spins[k] 
	u = exp(-s.beta * dEin)
	if rand() < u
		for j in C
			s.state.spins[j] *= -1
		end
	end

	# To be optimized
	compute_magnetization_and_energy!(s)
end

"Compute the observables."
function compute_observables!(s::Simu)
	s.obs.n += 1

	m = s.state.magnetization
	x = 1.
	for i in 1:length(s.obs.moms_mag)
		x *= m
		s.obs.moms_mag[i] += x
	end

	E = s.state.energy
	x = 1.
	for i in 1:length(s.obs.moms_energy)
		x *= E
		s.obs.moms_energy[i] += x
	end
end

"Run the simulation."
function run!(s::Simu, method::String, nit::Int, record::Bool=false,
		skip::Int=1)
	if method == "markov"
		for i in 1:nit
			for j in 1:(s.L*s.L)
				update_markov!(s)
			end
			if record && (i % skip == 0)
				compute_observables!(s)
			end
		end
	elseif method == "heatbath"
		for i in 1:nit
			for j in 1:(s.L*s.L)
				update_heatbath!(s)
			end
			if record && (i % skip == 0)
				compute_observables!(s)
			end
		end
	elseif method == "cluster"
		for i in 1:nit
			update_cluster!(s)
			if record && (i % skip == 0)
				compute_observables!(s)
			end
		end
	end
end

"Export observables to file."
function export_observables(fname::String, obs::Observables)
	open(fname; write=true) do f
		hd = ("# " * join(["m^$i" for i in 1:length(obs.moms_mag)], " ")
			  * " " * join(["m^$i" for i in 1:length(obs.moms_energy)], " ")
			  * "\n")
		write(f, hd)
		writedlm(f, [obs.moms_mag; obs.moms_energy], " ")
	end
end

function main()
	args = parse_commandline()

	if args["verbose"]
		println("Arguments:")
		println(args)
	end

	if !(args["method"] in ["markov", "heatbath", "cluster"])
		println("Method should be 'markov', 'heatbath' or 'cluster'.")
		return
	end

	obs = Observables(0, zeros(args["nmoms"]), zeros(args["nmoms"]))

	for n in 1:args["nsim"]
		if args["verbose"]
			print("\r", n)
		end
		simu = Simu(args)
		run!(simu, args["method"], args["nth"], false)
		run!(simu, args["method"], args["nit"], true, args["skip"])
		scale!(simu.obs)
		add!(obs, simu.obs)
		if args["spins"]
			writedlm(args["output"] * "_spins_$n.txt", simu.state.spins)
		end
	end
	
	scale!(obs)
	export_observables(args["output"] * "_obs.txt", obs)
end

main()
