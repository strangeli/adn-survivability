using ClusterManagers
using Distributed

if length(ARGS) > 0
	N_workers = parse(Int, ARGS[1]) - 1
else
	N_workers = 1
end

addprocs(SlurmManager(N_workers))

println(nprocs(), " process(es)")

@everywhere println("Testing slurm")
