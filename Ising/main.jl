using ArgParse
using DelimitedFiles

include("./Ising.jl")
using .Ising

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
			default = 1
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

function args_to_string(args::Dict)
	"# " * join(["$a=$v" for (a, v) in args], ", ") * "\n"
end

"Export observables to file."
function export_observables(fname::String, obs::Observables, hd::String="")
	open(fname; write=true) do f
		write(f, hd)
		hd2 = ("# " * join(["m^$i" for i in 1:length(obs.moms_mag)], " ")
		  	   * " " * join(["m^$i" for i in 1:length(obs.moms_energy)], " ")
			   * "\n")
		write(f, hd2)
		writedlm(f, [obs.moms_mag; obs.moms_energy], " ")
	end
end

function main()
	args = parse_commandline()
	hd = args_to_string(args)

	if args["verbose"]
		println("Arguments:")
		print(hd)
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
		simu = Simu(args["length"], args["beta"], args["J"], args["H"],
					args["nmoms"], args["determ"])
		run!(simu, args["method"], args["nth"], false)
		run!(simu, args["method"], args["nit"], true, args["skip"])
		scale!(simu.obs)
		add!(obs, simu.obs)
		if args["spins"]
			open(args["output"] * "_spins_$n.txt"; write=true) do f
				write(f, hd)
				writedlm(f, simu.state.spins)
			end
		end
	end
	
	scale!(obs)
	export_observables(args["output"] * "_obs.txt", obs)
end

main()
