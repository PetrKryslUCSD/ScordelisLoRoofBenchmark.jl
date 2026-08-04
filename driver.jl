using ArgParse
using DataDrop

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--verbose"
        help = "Verbosity level (0 = quiet, 1 = info, 2 = debug)"
        arg_type = Int
        default = 0
        "--N"
        help = "Numbers of edges along the diaphragm"
        arg_type = Int
        default = 2
        "--nref"
        help = "Number of refinements"
        arg_type = Int
        default = 4
        "--mesh"
        help = "Mesh type (uniform, graded, or biased)"
        arg_type = String
        default = "uniform"
        "--support"
        help = "Support type (soft, hard)"
        arg_type = String
        default = "hard"
        "--element"
        help = "Element type (q4rs, etc.)"
        arg_type = String
        default = "q4rs"
        "--action"
        help = "Action type (deflection, resultants)"
        arg_type = String
        default = "deflection"
    end
    return parse_args(s)
end

p = parse_commandline()

if p["verbose"] > 0
    @info "Line argument values:"
    for (k, v) in p
        @info "  $k = $v"
    end
end

using ScordelisLoRoofBenchmark

N = p["N"]
nref = p["nref"]
mesh = p["mesh"]
element = p["element"]
support = p["support"]
ns = [N^j for j in 1:nref]

p["verbose"] > 0 && @info "Number of elements along the diaphragm: $(ns)"

if p["action"] == "deflection"
    @info "Running the action deflection..."
    results = ScordelisLoRoofBenchmark.deflection(support=support, element=element, mesh=mesh, ns=ns, verbose=p["verbose"])
    d = Dict(
        "ns" => results[1],
        "deflections" => results[2],   
        "extrapolation" => results[3]
    )
    @info "Deflection results stored in deflection_results.json"
    DataDrop.store_json("deflection_results.json", d)
end


if p["action"] == "resultants"
    @info "Running the action resultants..."
    ScordelisLoRoofBenchmark.resultants(support=support, element=element, mesh=mesh, ns=ns, verbose=p["verbose"])
    @info "Resultants along edges saved to CSV files."
    @info "Use action=extrapolate to generate extrapolations."
    @info "Use action=plot to generate plots."
end


if p["action"] == "extrapolate"
    @info "Running the action resultants..."
    results = ScordelisLoRoofBenchmark.extrapolate_resultants(support=support, element=element, mesh=mesh, ns=ns, verbose=p["verbose"])
    @info "Resultants along edges saved to CSV files."
    @info "Use action=extrapolate to generate extrapolations."
    @info "Use action=plot to generate plots."
end

        
