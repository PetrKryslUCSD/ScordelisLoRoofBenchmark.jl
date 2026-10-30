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
        ["--thickness", "-t"],
        Dict(
            :help => "Thickness of the shell",
            :arg_type => Float64,
            :default => 0.25
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

thickness = p["thickness"]
N = p["N"]
nref = p["nref"]
mesh = p["mesh"]
element = p["element"]
support = p["support"]
ns = N .* [2^(j-1) for j in 1:nref]

p["verbosity"] > 0 && @info "Thickness: $(thickness)"
p["verbosity"] > 0 && @info "Number of elements along the diaphragm: $(ns)"

using ScordelisLoRoofBenchmark: compute_deflection, compute_resultants, extrapolate_resultants, plot_resultants

code = "scolo-$(thickness)-$(support)-$(mesh)-$(element)-N=$(N)-nref=$(nref)"
p["verbosity"] > 0 && @info "Code: $(code)"

if p["compute"] == "deflection"
    @info "Computing deflection..."
    results = compute_deflection(code; thickness=thickness, support=support, element=element, mesh=mesh, ns=ns, verbosity=p["verbosity"])
    storeas = "$(code)-deflection_extrapolations.json"
    @info "Deflection extrapolation stored in \n      $storeas"
    JSON.json(storeas, results; pretty=true, allownan=true)
end

if p["compute"] == "resultants"
    @info "Computing resultants..."
    compute_resultants(code; thickness=thickness, support=support, element=element, mesh=mesh, ns=ns, verbosity=p["verbosity"])
    @info "Resultants along edges saved to CSV files."
    @info "Computing extrapolations..."
    results = extrapolate_resultants(code; ns = ns, res = "n", verbosity=p["verbosity"])
    storeas = "$(code)-n_extrapolations.json"
    JSON.json(storeas, results; pretty=true, allownan=true)
    @info "Membrane force resultant extrapolations saved to \n      $storeas."
    results = extrapolate_resultants(code; ns = ns, res = "m", verbosity=p["verbosity"])
    storeas = "$(code)-m_extrapolations.json"
    JSON.json(storeas, results; pretty=true, allownan=true)
    @info "Membrane force resultant extrapolations saved to \n      $storeas."
    results = extrapolate_resultants(code; ns = ns, res = "q", verbosity=p["verbosity"])
    storeas = "$(code)-q_extrapolations.json"
    JSON.json(storeas, results; pretty=true, allownan=true)
    @info "Shear force resultant extrapolations saved to \n      $storeas."
    @info "Use --compute=plots to generate plots."
end

if p["compute"] == "plots"
    @info "Generating plots..."
    for n in ns
        basefn = "$(code)" * "-$(n)"
        plot_resultants(basefn=basefn, res="n", verbosity=p["verbosity"])
        plot_resultants(basefn=basefn, res="m", verbosity=p["verbosity"])
        plot_resultants(basefn=basefn, res="q", verbosity=p["verbosity"])
    end
    @info "Resultants along edges saved to PDF files."
end

        
