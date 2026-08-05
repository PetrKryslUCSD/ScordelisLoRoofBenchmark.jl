"""
The barrel vault (Scordelis-Lo) roof is one of the benchmarks for linear elastic
analysis of shells. 

# Problem description

The physical basis of the problem is a deeply arched roof supported only
by diaphragms at its curved edges (an aircraft hangar), deforming under its own
weight. It is interesting to observe that the geometry is such that the
centerpoint of the roof moves upward under the self-weight (downwardly directed)
load. Perhaps this is one reason why the problem is not straightforward
numerically.
"""
module ScordelisLoRoofBenchmark

using LinearAlgebra
using FinEtools
using FinEtools.AlgoBaseModule: solve_blocked!
using FinEtools.AlgoBaseModule: richextrapol
using FinEtoolsDeforLinear
using FinEtoolsFlexStructures.FESetShellT3Module: FESetShellT3
using FinEtoolsFlexStructures.FESetShellQ4Module: FESetShellQ4
using FinEtoolsFlexStructures.FEMMShellT3FFModule
using FinEtoolsFlexStructures.FEMMShellQ4RSModule
using FinEtoolsFlexStructures.RotUtilModule: initial_Rfield, update_rotation_field!
using FinEtools.MeshExportModule: VTKWrite
using FinEtools.MeshExportModule.VTKWrite: vtkwrite
using FinEtools.MeshExportModule.CSV: savecsv
using RichardsonExtrapolationUQ
using DelimitedFiles
using PGFPlotsX

# Parameters:
const E = 4.32e8;
const nu = 0.0;
const thickness = 0.25; # geometrical dimensions are in feet
const R = 25.0;
const L = 50.0;

cylindrical!(csmatout, XYZ, tangents, feid, qpid) = begin
    r = vec(deepcopy(XYZ));
    r[2] = 0.0;
    r[3] += R # see def of Z below
    csmatout[:, 3] .= vec(r)/norm(vec(r))
    csmatout[:, 2] .= (0.0, 1.0, 0.0) #  this is along the axis
    cross3!(view(csmatout, :, 1), view(csmatout, :, 2), view(csmatout, :, 3))
    return csmatout
end

function _exe(;basef="scolo", element=:q4rs, mesh=:uniform, n=8, support = :soft, deflection_only=false, verbosity=0)
    verbosity > 0 && @info "Element: $element; Support: $support; Mesh: $mesh, $n elements per side"
    bias = 100
    if mesh == :uniform
        fens, fes = Q4block(40/360*2*pi, L/2, n, 2*n);
        tolerance = L/2/n/100
    elseif mesh == :biased
        xs = 40/360*2*pi .- reverse(biasedspace(0.0, 40/360*2*pi, n+1, bias))
        ys = biasedspace(0.0, L/2, 2*n+1, bias)
        fens, fes = Q4blockx(xs, ys);
        tolerance = (minimum(abs.(diff(xs))) + minimum(abs.(diff(ys)))) / 100
    else 
        @error "Unknown mesh"
    end
    if element == :t3ff
        fens, fes = Q4toT3(fens, fes);
    end
    verbosity > 0 && @info "Number of nodes: $(count(fens)); Number of elements: $(count(fes))"
    bfes = meshboundary(fes)
    ela0 = selectelem(fens, bfes; facing=true, direction=Float64[-1, 0])
    ela1 = selectelem(fens, bfes; facing=true, direction=Float64[+1, 0])
    ell0 = selectelem(fens, bfes; facing=true, direction=Float64[0, -1])
    ell1 = selectelem(fens, bfes; facing=true, direction=Float64[0, +1])
    fens.xyz = xyz3(fens)
    for i in 1:count(fens)
        a=fens.xyz[i, 1];
        y=fens.xyz[i, 2];
        fens.xyz[i, :] .= (R*sin(a), y, R*(cos(a)-1))
    end

    mater = MatDeforElastIso(DeforModelRed3D, E, nu)
    ocsys = CSys(3, 3, cylindrical!)

    if Symbol(element) == :q4rs
        sfes = FESetShellQ4()
    elseif Symbol(element) == :t3ff
        sfes = FESetShellT3()
    end
    accepttodelegate(fes, sfes)
    if Symbol(element) == :q4rs
        femm = FEMMShellQ4RSModule.make(IntegDomain(fes, GaussRule(2, 2), thickness), ocsys, mater)
        stiffness = FEMMShellQ4RSModule.stiffness
        associategeometry! = FEMMShellQ4RSModule.associategeometry!
    elseif Symbol(element) == :t3ff
        femm = FEMMShellT3FFModule.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
        stiffness = FEMMShellT3FFModule.stiffness
        associategeometry! = FEMMShellT3FFModule.associategeometry!
    end
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz, 1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz, 1), 6))

    # Apply EBC's
    # rigid diaphragm: SOFT simple support
    l1 = connectednodes(subset(bfes, ell0))
    dof = support == :soft ? [1, 3] : [1, 3, 5]
    for i in dof
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to Y
    l1 = connectednodes(subset(bfes, ell1))
    for i in [2, 4, 6]
        setebc!(dchi, l1, true, i)
    end
    # plane of symmetry perpendicular to X
    l1 = connectednodes(subset(bfes, ela0))
    for i in [1, 5, 6]
        setebc!(dchi, l1, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    massem = SysmatAssemblerFFBlock(nfreedofs(dchi))
    vassem = SysvecAssemblerFBlock(nfreedofs(dchi))

    # Assemble the system matrix
    associategeometry!(femm, geom0)
    # vtkwrite("debug-normals-q4rs-$(n).vtu", fens, fes; vectors = [("normals", deepcopy(femm._normals[:, 1:3]))])
    Kff = stiffness(femm, massem, geom0, u0, Rfield0, dchi);

    # Midpoint of the free edge
    nl = selectnode(fens; box=Float64[sin(40/360*2*pi)*25 sin(40/360*2*pi)*25 L/2 L/2 -Inf Inf], inflate=tolerance)
    lfemm = if Symbol(element) == :q4rs
        FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    elseif Symbol(element) == :t3ff
        FEMMBase(IntegDomain(fes, TriRule(1)))
    end
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    
    deflatB = dchi.values[nl, 3][1]
    verbosity > 0 && @info "Deflection at B: $(round(deflatB, digits = 7))"

    if deflection_only
        return deflatB
    end
    
    verbosity > 0 && @info "Generating resultants and visualizations"

    midsection_nodes = connectednodes(subset(bfes, ell1))
    midsection_nodes_x = sortperm(fens.xyz[midsection_nodes, 1])
    midsection_nodes_ordered = midsection_nodes[midsection_nodes_x]
    midsection_angles = asin.(fens.xyz[midsection_nodes_ordered, 1] ./ R) * 180 / pi / 40

    diaphragm_nodes = connectednodes(subset(bfes, ell0))
    diaphragm_nodes_x = sortperm(fens.xyz[diaphragm_nodes, 1])
    diaphragm_nodes_ordered = diaphragm_nodes[diaphragm_nodes_x]
    diaphragm_angles = asin.(fens.xyz[diaphragm_nodes_ordered, 1] ./ R) * 180 / pi / 40

    peak_nodes = connectednodes(subset(bfes, ela0))
    peak_nodes_y = sortperm(fens.xyz[peak_nodes, 2])
    peak_nodes_ordered = peak_nodes[peak_nodes_y]
    peak_dists = fens.xyz[peak_nodes_ordered, 2] ./ (L/2)

    free_nodes = connectednodes(subset(bfes, ela1))
    free_nodes_y = sortperm(fens.xyz[free_nodes, 2])
    free_nodes_ordered = free_nodes[free_nodes_y]
    free_dists = fens.xyz[free_nodes_ordered, 2] ./ (L/2)

    # Visualization
    verbosity > 0 && @info "Generating visualizations"
    basef = basef * "-$(element)-$(support)-$(mesh)-$(n)"
    # Generate a graphical display of resultants
    pfemm = if Symbol(element) == :q4rs
        # Here we generate an ad hoc machine: based on a single integration point, it improves the post processing of the stresses, especially the shear membrane.
        pfemm = FEMMShellQ4RSModule.make(IntegDomain(fes, GaussRule(2, 1), thickness), ocsys, mater)
        associategeometry!(pfemm, geom0)
    elseif Symbol(element) == :t3ff
        femm
    end
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :moment, nc, outputcsys=ocsys)
        push!(scalars, ("m$nc", fld.values))
        verbosity > 0 && @info "m$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        verbosity > 0 && @info "Saving CSV files for m$nc"
        savecsv("$(basef)-midsection-m$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-m$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-m$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-m$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    verbosity > 0 && @info "Saving vtu files for m resultants"
    vtkwrite("$(basef)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
        push!(scalars, ("n$nc", fld.values))
        verbosity > 0 && @info "n$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        verbosity > 0 && @info "Saving CSV files for n$nc"
        savecsv("$(basef)-midsection-n$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-n$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-n$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-n$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    verbosity > 0 && @info "Saving vtu files for n resultants"
    vtkwrite("$(basef)-n.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        push!(scalars, ("q$nc", fld.values))
        verbosity > 0 && @info "q$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        verbosity > 0 && @info "Saving CSV files for q$nc"
        savecsv("$(basef)-midsection-q$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-q$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-q$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-q$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    verbosity > 0 && @info "Saving vtu files for q resultants"
    vtkwrite("$(basef)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    
    return nothing
end

const COLORS = ["red", "green", "blue", "black", "cyan", "magenta", "yellow", "gray"]
const MARKERS = [
    "diamond",
    "triangle",
    "square",  
    "triangle*",
    "x",
    "square*",
    "diamond*",
]
const MARK_REPEAT = [37, 43, 57, 79, 23, 27, 29] .+ 173

function plot_resultants(basef = "", 
              res = "q")
    ncs = res == "q" ? (1:2) : (1:3)
    for edge in ["midsection", "diaphragm", "peak", "free"]
        objects = []
        for nc in ncs
            f = "$(basef)-$(edge)-$(res)$(nc).csv"
            A = readdlm(f, ',', Float64; skipstart=1)
            @pgf o = PGFPlotsX.Plot(
            {
            color = COLORS[nc],
            mark=MARKERS[nc], mark_size=2.5, mark_repeat=MARK_REPEAT[nc],
            line_width  = 1.0
            },
            Coordinates([v for v in  zip(A[:,1], A[:,2])])
            )
            push!(objects, o)
            # push!(objects, LegendEntry("$(res)$(nc)"))
        end
        @pgf ax = Axis(
            {
                title = "$(edge)",
                xlabel = "Normalized distance",
                ylabel = "Resultant",
                # xmin = -0.01,
                # xmax = 1.01,
                xmode = "linear",
                ymode = "linear",
                yminorgrids = "true",
                grid = "both",
                # legend_style = {
                #     at = Coordinate(1.005, 0.5),
                #     anchor = "west",
                # },
            },
            objects...
        )
        display(ax)
        pgfsave("$(basef)-$(edge)-$(res).pdf", ax)
    end
    return true
end

const CORNERS = [
    (:A, "midsection", :b),
    (:B, "midsection", :e),
    (:C, "diaphragm", :e),
    (:D, "diaphragm", :b),
    ]

function extrapolate_resultants(;basef = "", ns = [16, 32, 64, 128], res = "q", verbosity=0)
    ncs = res == "q" ? (1:2) : (1:3)
    results = Dict()
    for corner in CORNERS
        edge = corner[2]
        components = Dict()
        for nc in ncs
            r = Float64[]
            for n in ns
                f = "$(basef)-$(n)-$(edge)-$(res)$(nc).csv"
                A = readdlm(f, ',', Float64; skipstart=1)
                push!(r, corner[3] == :b ? A[1, 2] : A[end, 2])
            end
            e = nothing
            try    
                e = RichardsonExtrapolationUQ.richextrapol_uq(r, 1.0 ./ ns)
            catch
                e = (success = false, message = "no result")
            end
            if e.success
                verbosity > 0 && @info "$(corner[1]): $(res)$(nc)"
                verbosity > 0 && @info "q_ci = $(e.q_star) +/- $(e.q_star_ci), beta_ci = $(e.beta_star) +/- $(e.beta_star_ci)"
            else
                verbosity > 0 && @warn "$(corner[1]): $(res)$(nc) - no result"
            end
            components[nc] = e
        end
        results[corner[1]] = components
    end
    for edge in ["diaphragm"]#["midsection", "diaphragm", "peak", "free"]
        components = Dict()
        for nc in ncs
            rmx = Float64[]
            rmn = Float64[]
            for n in ns
                f = "$(basef)-$(n)-$(edge)-$(res)$(nc).csv"
                A = readdlm(f, ',', Float64; skipstart=1)
                push!(rmx, maximum(A[:, 2]))
                push!(rmn, minimum(A[:, 2]))
            end
            r = rmx
            e = nothing
            try    
                e = RichardsonExtrapolationUQ.richextrapol_uq(r, 1.0 ./ ns)
            catch
                e = (success = false, message = "no result")
            end
            if e.success
                verbosity > 0 && @info "$(edge): $(res)$(nc) maximum"
                verbosity > 0 && @info "q_ci = $(e.q_star) +/- $(e.q_star_ci), beta_ci = $(e.beta_star) +/- $(e.beta_star_ci)"
            else
                verbosity > 0 && @warn "$(edge): $(res)$(nc) maximum - no result"
            end
            components[nc] = e
            r = rmn
            e = nothing
            try    
                e = RichardsonExtrapolationUQ.richextrapol_uq(r, 1.0 ./ ns)
            catch
                e = (success = false, message = "no result")
            end
            if e.success
                verbosity > 0 && @info "$(edge): $(res)$(nc) minimum"
                verbosity > 0 && @info "q_ci = $(e.q_star) +/- $(e.q_star_ci), beta_ci = $(e.beta_star) +/- $(e.beta_star_ci)"
            else
                verbosity > 0 && @warn "$(edge): $(res)$(nc) minimum - no result"
            end
        end
        results[Symbol(edge)] = components
    end
    return results
end

function resultants(;ns=[128, 256, 512, 1024], mesh=:uniform, element=:q4rs, support=:soft, verbosity=0)
    deflection_only = false
    if Symbol(element) == :q4rs
    elseif Symbol(element) == :t3ff
    else
        throw(ArgumentError("Unsupported element type"))
    end
    for n in ns
        _exe(element=Symbol(element), mesh=Symbol(mesh), n=n, support=Symbol(support), deflection_only=deflection_only, verbosity=verbosity)
    end
    return nothing
end

function deflection(;ns=4 .* [16, 32, 64, 128, ], mesh=:uniform, element=:q4rs, support=:soft, verbosity=0)
    deflection_only = true 
    deflectionsB = Float64[]
    if Symbol(element) == :q4rs
    elseif Symbol(element) == :t3ff
    else
        throw(ArgumentError("Unsupported element type"))
    end
    for n in ns
        v = _exe(element=Symbol(element), mesh=Symbol(mesh), n=n, support=Symbol(support), deflection_only=deflection_only, verbosity=verbosity)
        push!(deflectionsB, v)
    end
    if length(ns) > 3
        extrapolation = RichardsonExtrapolationUQ.richextrapol_uq(deflectionsB, 1.0 ./ ns)
    elseif length(ns) > 2
        extrapolation = richextrapol(deflectionsB, 1.0 ./ ns)
    else
        extrapolation = nothing
    end
    return Dict(
        "ns" => ns,
        "deflections" => deflectionsB,   
        "results_of_extrapolation" => extrapolation
    )
end

nothing

end # module ScordelisLoRoofBenchmark
