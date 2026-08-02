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

function _execute_q4rs(mesh=:uniform, n=8, support = :soft, deflection_only=false)
    formul = FEMMShellQ4RSModule
    @info "Support: $support; Mesh: $mesh, $n elements per side"
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
    @info "Number of nodes: $(count(fens)); Number of elements: $(count(fes))"
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

    sfes = FESetShellQ4()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, GaussRule(2, 2), thickness), ocsys, mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

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
    lfemm = FEMMBase(IntegDomain(fes, GaussRule(2, 2)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)

    result = dchi.values[nl, 3][1]
    @info "Deflection at B: $(round(result, digits = 7))"

    if !deflection_only
        @info "Generating resultants and visualizations"
    else
        return result
    end

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
    basef = "scolo_q4rs-$(support)-$(mesh)-$(n)"
    # Generate a graphical display of resultants
    # Here we generate an ad hoc machine: based on a single integration point, it improves the post processing of the stresses, especially the shear membrane.
    pfemm = formul.make(IntegDomain(fes, GaussRule(2, 1), thickness), ocsys, mater)
    associategeometry!(pfemm, geom0)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :moment, nc, outputcsys=ocsys)
        push!(scalars, ("m$nc", fld.values))
        @info "m$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-m$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-m$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-m$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-m$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    vtkwrite("$(basef)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
        push!(scalars, ("n$nc", fld.values))
        @info "n$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-n$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-n$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-n$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-n$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    vtkwrite("$(basef)-n.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        push!(scalars, ("q$nc", fld.values))
        @info "q$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-q$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-q$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-q$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-q$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    vtkwrite("$(basef)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    
    result
end

function _execute_t3ff(mesh=:uniform, n=8, support = :soft, deflection_only=false)
    formul = FEMMShellT3FFModule
    @info "Mesh: $mesh; $n elements per side"
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
    fens, fes = Q4toT3(fens, fes);
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

    sfes = FESetShellT3()
    accepttodelegate(fes, sfes)
    femm = formul.make(IntegDomain(fes, TriRule(1), thickness), ocsys, mater)
    stiffness = formul.stiffness
    associategeometry! = formul.associategeometry!

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
    lfemm = FEMMBase(IntegDomain(fes, TriRule(1)))
    fi = ForceIntensity(Float64[0, 0, -90, 0, 0, 0]);
    Ff = distribloads(lfemm, vassem, geom0, dchi, fi, 2);

    # Solve
    Uf = Kff \ Ff
    scattersysvec!(dchi, Uf, DOF_KIND_FREE)
    U = gathersysvec(dchi, DOF_KIND_ALL)

    result = dchi.values[nl, 3][1]
    @info "Solution: $(result), $(round(result, digits = 7))"

    if !deflection_only
        @info "Generating resultants and visualizations"
    else
        return result
    end

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
    basef = "scolo_t3ff-$(support)-$(mesh)-$(n)"
    # Generate a graphical display of resultants
    # No point in generating another machine; this is the best we can do.
    pfemm = femm
    associategeometry!(pfemm, geom0)
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :moment, nc, outputcsys=ocsys)
        push!(scalars, ("m$nc", fld.values))
        @info "m$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-m$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-m$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-m$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-m$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    vtkwrite("$(basef)-m.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    scalars = []
    for nc in 1:3
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :membrane, nc, outputcsys=ocsys)
        push!(scalars, ("n$nc", fld.values))
        @info "n$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-n$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-n$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-n$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-n$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    vtkwrite("$(basef)-n.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    # for nc in 1:3
    #     scalars = []
    #     connectivity, points, values = _elementwise_arrays(pfemm, geom0, dchi, :membrane, nc, tolerance)
    #     push!(scalars, ("dn$nc", deepcopy(values)))
    #     vtkwrite("$(basef)-dn$nc.vtu", connectivity, points, VTKWrite.Q4; scalars=scalars)
    # end
    scalars = []
    for nc in 1:2
        fld = fieldfromintegpoints(pfemm, geom0, dchi, :shear, nc, outputcsys=ocsys)
        push!(scalars, ("q$nc", fld.values))
        @info "q$nc Range: $(minimum(fld.values)) to $(maximum(fld.values))"
        savecsv("$(basef)-midsection-q$(nc).csv", a=midsection_angles, v=fld.values[midsection_nodes_ordered])
        savecsv("$(basef)-diaphragm-q$(nc).csv", a=diaphragm_angles, v=fld.values[diaphragm_nodes_ordered])
        savecsv("$(basef)-peak-q$(nc).csv", a=peak_dists, v=fld.values[peak_nodes_ordered])
        savecsv("$(basef)-free-q$(nc).csv", a=free_dists, v=fld.values[free_nodes_ordered])
    end
    vtkwrite("$(basef)-q.vtu", fens, fes; scalars=scalars, vectors=[("u", dchi.values[:, 1:3])])
    
    result
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

function extrapolate_resultants(basef = "", ns = [16, 32, 64, 128], res = "q")
    ncs = res == "q" ? (1:2) : (1:3)
    results = [] 
    for corner in CORNERS
        edge = corner[2]
        for nc in ncs
            r = Float64[]
            for n in ns
                f = "$(basef)-$(n)-$(edge)-$(res)$(nc).csv"
                A = readdlm(f, ',', Float64; skipstart=1)
                push!(r, corner[3] == :b ? A[1, 2] : A[end, 2])
            end
            extrapolation = nothing
            try    
                extrapolation = RichardsonExtrapolationUQ.richextrapol_uq(r, 1.0 ./ ns)# richextrapol(r, 1.0 ./ ns)
            catch
                extrapolation = (data = r,)
            end
            push!(results, (corner = corner[1], nc = nc, extrapolation = extrapolation))
        end
    end
    for edge in ["midsection", "diaphragm", "peak", "free"]
        for nc in ncs
            rmx = Float64[]
            rmn = Float64[]
            for n in ns
                f = "$(basef)-$(n)-$(edge)-$(res)$(nc).csv"
                A = readdlm(f, ',', Float64; skipstart=1)
                push!(rmx, maximum(A[:, 2]))
                push!(rmn, minimum(A[:, 2]))
            end
            extrapolationmx = nothing
            try    
                extrapolationmx = RichardsonExtrapolationUQ.richextrapol_uq(rmx, 1.0 ./ ns)# richextrapol(r, 1.0 ./ ns)
            catch
                extrapolationmx = (data = rmx,)
            end
            extrapolationmn = nothing
            try    
                extrapolationmn = RichardsonExtrapolationUQ.richextrapol_uq(rmn, 1.0 ./ ns)# richextrapol(r, 1.0 ./ ns)
            catch
                extrapolationmn = (data = rmn,)
            end
            push!(results, (edge = Symbol(edge), nc = nc, extrapolationmx = extrapolationmx, extrapolationmn = extrapolationmn))
        end
    end
    return results
end

function resultants(;ns=[128, 256, 512, 1024], mesh=:uniform, element=:q4rs, support=:soft)
    deflection_only = false
    for n in ns
        if Symbol(element) == :q4rs
            v = _execute_q4rs(Symbol(mesh), n, Symbol(support), deflection_only)
        elseif Symbol(element) == :t3ff

            v = _execute_t3ff(Symbol(mesh), n, Symbol(support), deflection_only)
        else
            throw(ArgumentError("Unsupported element type"))
        end
    end
    @info "Resultants along edges saved to CSV files. Use plot_resultants() to generate plots."
    return ns
end

function deflection(;ns=4 .* [16, 32, 64, 128, ], mesh=:uniform, element=:q4rs, support=:soft, deflection_only = true)
    deflectionsB = Float64[]
    for n in ns
        if Symbol(element) == :q4rs
            v = _execute_q4rs(Symbol(mesh), n, Symbol(support), deflection_only)
        elseif Symbol(element) == :t3ff
            v = _execute_t3ff(Symbol(mesh), n, Symbol(support), deflection_only)
        else
            throw(ArgumentError("Unsupported element type"))
        end
        push!(deflectionsB, v)
    end
    if length(ns) > 3
        extrapolation = RichardsonExtrapolationUQ.richextrapol_uq(deflectionsB, 1.0 ./ ns)
    elseif length(ns) > 2
        extrapolation = richextrapol(deflectionsB, 1.0 ./ ns)
    else
        extrapolation = nothing
    end
    return ns, deflectionsB, extrapolation
end

nothing

end # module ScordelisLoRoofBenchmark
