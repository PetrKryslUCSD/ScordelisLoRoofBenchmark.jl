using Pkg
Pkg.activate(".", io=devnull)
Pkg.instantiate()

using ArgParse
using JSON
using ScordelisLoRoofBenchmark

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--verbosity"
        help = "Verbosity level (0 = quiet, 1 = info, 2 = debug)"
        arg_type = Int
        default = 1
        "--N"
        help = "Numbers of edges along the diaphragm"
        arg_type = Int
        default = 4
        "--nref"
        help = "Number of refinements (for extrapolation study, 3 or 4 (for uncertainty quantification); 1 = no refinement)"
        arg_type = Int
        default = 4
        "--mesh"
        help = "Mesh type (uniform or biased)"
        arg_type = String
        default = "uniform"
        "--support"
        help = "Support type (soft or hard)"
        arg_type = String
        default = "hard"
        "--element"
        help = "Element type (q4rs or t3ff)"
        arg_type = String
        default = "q4rs"
        "--action"
        help = "Action type (deflection, resultants, extrapolate, plot)"
        arg_type = String
        default = ""
        "--arguments"
        help = "Path to the arguments JSON file"
        arg_type = String
        default = ""
    end
    return parse_args(s)
end

p = parse_commandline()

if p["action"] == ""
    println("Help is available: julia --project driver.jl --help")
    exit(0)
end



if p["arguments"] != ""
    p["verbosity"] > 0 && @info "Loading arguments from $(p["arguments"])"
    action = p["action"]
    p = open(p["arguments"], "r") do file
        JSON.parse(file)
    end
    p["action"] = action 
end

if p["verbosity"] > 0
    @info "Line argument values:"
    for (k, v) in p
        @info "  $k = $v"
    end
end

N = p["N"]
nref = p["nref"]
mesh = p["mesh"]
element = p["element"]
support = p["support"]
ns = [N^j for j in 1:nref]

p["verbosity"] > 0 && @info "Number of elements along the diaphragm: $(ns)"

if p["action"] == "deflection"
    @info "Running the action deflection..."
    results = ScordelisLoRoofBenchmark.deflection(support=support, element=element, mesh=mesh, ns=ns, verbosity=p["verbosity"])
    @info "Deflection results stored in deflection_results.json"
    DataDrop.store_json("deflection_results.json", results)
    open("deflection_results.json", "w") do file
        JSON.print(file, results, 4)
    end
    @info "Arguments saved in arguments.json"
    open("arguments.json", "w") do file
        JSON.print(file, results, 4)
    end
    @info "The next run (with --action=OTHERACTION) may use\n the arguments file as --arguments=arguments.json"
end

if p["action"] == "resultants"
    @info "Running the action resultants..."
    ScordelisLoRoofBenchmark.resultants(support=support, element=element, mesh=mesh, ns=ns, verbosity=p["verbosity"])
    @info "Resultants along edges saved to CSV files."
    @info "Arguments saved in arguments.json"
    DataDrop.store_json("arguments.json", p)
    @info "Use --action=extrapolate to generate extrapolations."
    @info "Use --action=plot to generate plots."
    @info "The next run (with --action=extrapolate or --action=plot) may use\n the arguments file as --arguments=arguments.json"
end

if p["action"] == "extrapolate"
    @info "Running the action extrapolate..."
    basef = "scolo-$(element)-$(support)-$(mesh)"
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(basef = basef, ns = ns, res = "n", verbosity=p["verbosity"])
    JSON.json("n_extrapolations.json", results; pretty=true, allownan=true)
    @info "Membrane force resultant extrapolations saved to n_extrapolations.json."
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(basef = basef, ns = ns, res = "m", verbosity=p["verbosity"])
    JSON.json("m_extrapolations.json", results; pretty=true, allownan=true)
    @info "Membrane force resultant extrapolations saved to m_extrapolations.json."
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(basef = basef, ns = ns, res = "q", verbosity=p["verbosity"])
    JSON.json("q_extrapolations.json", results; pretty=true, allownan=true)
    @info "Shear force resultant extrapolations saved to q_extrapolations.json."
    @info "Use --action=plot to generate plots."
    @info "The next run (with --action=extrapolate or --action=plot) may use\n the arguments file as --arguments=arguments.json"
end

        
