using Pkg
Pkg.activate(".", io=devnull)
Pkg.instantiate()

using ArgParse
using JSON


function parse_commandline()
    s = ArgParseSettings()
    add_arg_table!(s,
        ["--verbosity", "-v"],
        Dict(
            :help => "Verbosity level (0 = quiet, 1 = info, 2 = debug)",
            :arg_type => Int,
            :default => 1
        ),
        ["--N", "-N"],
        Dict(
            :help => "Numbers of edges along the diaphragm",
            :arg_type => Int,
            :default => 4
        ),
        ["--nref", "-n"],
        Dict(
            :help => "Number of refinements (for extrapolation study, 3 or 4 (for uncertainty quantification); 1 = no refinement)",
            :arg_type => Int,
            :default => 4
        ),
        ["--mesh", "-m"],
        Dict(
            :help => "Mesh type (uniform or biased)",
            :arg_type => String,
            :default => "biased"
        ),
        ["--support", "-s"],
        Dict(
            :help => "Support type (soft or hard)",
            :arg_type => String,
            :default => "hard"
        ),
        ["--element", "-e"],
        Dict(
            :help => "Element type (q4rs or t3ff)",
            :arg_type => String,
            :default => "q4rs"
        ),
        ["--compute", "-c"],
        Dict(
            :help => "Compute type (deflection, resultants, plots)",
            :arg_type => String,
            :default => ""
        ),
    )
    return parse_args(s)
end

p = parse_commandline()

if p["compute"] == ""
    println("Help is available: julia --project driver.jl --help")
    exit(0)
end

if p["verbosity"] > 0
    @info "Argument values:"
    for (k, v) in p
        @info "  $k = $v"
    end
end

N = p["N"]
nref = p["nref"]
mesh = p["mesh"]
element = p["element"]
support = p["support"]
ns = N .* [2^(j-1) for j in 1:nref]

p["verbosity"] > 0 && @info "Number of elements along the diaphragm: $(ns)"

using ScordelisLoRoofBenchmark

code = "$(support)-$(mesh)-$(element)-$(N)-$(nref)"

if p["compute"] == "deflection"
    @info "Computing deflection..."
    results = ScordelisLoRoofBenchmark.compute_deflection(support=support, element=element, mesh=mesh, ns=ns, verbosity=p["verbosity"])
    storeas = "deflection_extrapolations-$(code).json"
    @info "Deflection extrapolation stored in $storeas"
    JSON.json(storeas, results; pretty=true, allownan=true)
end

if p["compute"] == "resultants"
    @info "Computing resultants..."
    ScordelisLoRoofBenchmark.compute_resultants(support=support, element=element, mesh=mesh, ns=ns, verbosity=p["verbosity"])
    @info "Resultants along edges saved to CSV files."
    @info "Computing extrapolations..."
    basef = "scolo-$(element)-$(support)-$(mesh)"
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(basef = basef, ns = ns, res = "n", verbosity=p["verbosity"])
    storeas = "n_extrapolations-$(code).json"
    JSON.json(storeas, results; pretty=true, allownan=true)
    @info "Membrane force resultant extrapolations saved to $storeas."
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(basef = basef, ns = ns, res = "m", verbosity=p["verbosity"])
    storeas = "m_extrapolations-$(code).json"
    JSON.json(storeas, results; pretty=true, allownan=true)
    @info "Membrane force resultant extrapolations saved to $storeas."
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(basef = basef, ns = ns, res = "q", verbosity=p["verbosity"])
    storeas = "q_extrapolations-$(code).json"
    JSON.json(storeas, results; pretty=true, allownan=true)
    @info "Shear force resultant extrapolations saved to $storeas."
    @info "Use --compute=plots to generate plots."
end

if p["compute"] == "plots"
    @info "Generating plots..."
    for n in ns
        basefn = "scolo-$(element)-$(support)-$(mesh)" * "-$(n)"
        ScordelisLoRoofBenchmark.plot_resultants(basefn=basefn, res="n", verbosity=p["verbosity"])
        ScordelisLoRoofBenchmark.plot_resultants(basefn=basefn, res="m", verbosity=p["verbosity"])
        ScordelisLoRoofBenchmark.plot_resultants(basefn=basefn, res="q", verbosity=p["verbosity"])
    end
    @info "Resultants along edges saved to PDF files."
end

        
