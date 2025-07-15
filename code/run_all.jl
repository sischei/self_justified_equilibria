#!/usr/bin/env julia
################################################################################
# run_all.jl
#
# Usage:
#   julia run_all.jl           # runs for Ψ = 0, 0.5, 0.9
#   julia run_all.jl 0.5       # runs only for Ψ = 0.5
################################################################################

# 1) Decide which Ψs to run
psi_list = length(ARGS) >= 1 ? [ parse(Float64, a) for a in ARGS ] : [0.0, 0.5, 0.9]

# 2) Loop over each Ψ, call mainhom.jl
for psi in psi_list
    println("\n\n===================== Running mainhom.jl for Ψ = $psi =====================")
    # Dispatch the same Julia executable, passing the argument through
    julia_cmd = Base.julia_cmd()                 # ensures we call this same julia
    cmd = `$(julia_cmd) --threads=$(Threads.nthreads()) mainhom.jl $psi`
    run(cmd)

    println(">> Done Ψ = $psi")
end

println("\nAll simulations complete.")
