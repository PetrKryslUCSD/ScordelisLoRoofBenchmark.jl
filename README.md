# ScordelisLoRoofBenchmark.jl

The barrel vault (Scordelis-Lo) roof benchmark implementation.

## Running

First of all, one needs an installation of [Julia](julialang.org).

Run the file `driver.jl` like this to obtain help:
```
julia --project driver.jl --help
```

One can control 
- What to compute (deflection, resultants, plots).
- Support over the diaphragm (hard or soft).
- The number of finite element edges along the diaphragm edge.
- Finite element to employ (quadrilateral q4rs or triangle t3ff).
- The mesh type  (biased or biased).
- The number of refinements of the mesh (4 to get the results with 
    uncertainty quantification in the form of a 95% confidence interval).
    The refinements are by bisection.
- The verbosity of the output (0=minimal, 1=more verbose).

Any command will activate and instantiate the environment,
which may lead to a longish pre-compilation period.

To compute just the deflection at the point B, run the following:
```
julia --project=. driver.jl -c deflection
```
This command is executed with the default settings.
The following output will be printed:
```
[ Info: Argument values:
[ Info:   support = hard
[ Info:   element = q4rs
[ Info:   N = 4
[ Info:   mesh = biased
[ Info:   verbosity = 1
[ Info:   compute = deflection
[ Info:   nref = 4
[ Info: Number of elements along the diaphragm: [4, 8, 16, 32]
[ Info: Computing deflection...
[ Info: Element: q4rs; Support: hard; Mesh: biased, 4 elements per side
[ Info: Number of nodes: 45; Number of elements: 32
[ Info: Deflection at B: -0.2454529
[ Info: Element: q4rs; Support: hard; Mesh: biased, 8 elements per side
[ Info: Number of nodes: 153; Number of elements: 128
[ Info: Deflection at B: -0.2744281
[ Info: Element: q4rs; Support: hard; Mesh: biased, 16 elements per side
[ Info: Number of nodes: 561; Number of elements: 512
[ Info: Deflection at B: -0.2958704
[ Info: Element: q4rs; Support: hard; Mesh: biased, 32 elements per side
[ Info: Number of nodes: 2145; Number of elements: 2048
[ Info: Deflection at B: -0.3001441
[ Info: Deflection extrapolation stored in deflection_extrapolations-hard-biased-q4rs-4-4.json
```

To choose something else, use the argument flags. For instance:
```
julia --project=. driver.jl -s hard -e t3ff -N 16 -c deflection
```
will select hard simple support over the diaphragm, the triangle element,
and the initial mesh will have 16 edges along the diaphragm and 2*16 edges axially.
This time, the printouts will be:
```
[ Info: Argument values:
[ Info:   support = hard
[ Info:   element = t3ff
[ Info:   N = 16
[ Info:   mesh = biased
[ Info:   verbosity = 1
[ Info:   compute = deflection
[ Info:   nref = 4
[ Info: Number of elements along the diaphragm: [16, 32, 64, 128]
[ Info: Computing deflection...
[ Info: Element: t3ff; Support: hard; Mesh: biased, 16 elements per side
[ Info: Number of nodes: 561; Number of elements: 1024
[ Info: Deflection at B: -0.2742959
[ Info: Element: t3ff; Support: hard; Mesh: biased, 32 elements per side
[ Info: Number of nodes: 2145; Number of elements: 4096
[ Info: Deflection at B: -0.293887
[ Info: Element: t3ff; Support: hard; Mesh: biased, 64 elements per side
[ Info: Number of nodes: 8385; Number of elements: 16384
[ Info: Deflection at B: -0.2994473
[ Info: Element: t3ff; Support: hard; Mesh: biased, 128 elements per side
[ Info: Number of nodes: 33153; Number of elements: 65536
[ Info: Deflection at B: -0.3009366
[ Info: Deflection extrapolation stored in deflection_extrapolations-hard-biased-t3ff-16-4.json
```

The compute the resultants and their extrapolations, use this command line:
```
julia --project=. driver.jl -s hard -e t3ff -N 16 -c resultants
```
This will result in voluminous output. To suppress that output, add `-v=0`:
```
julia --project=. driver.jl -s hard -e t3ff -N 16 -c resultants -v=0
```
That will result in the more manageable output
```
[ Info: Computing resultants...
[ Info: Resultants along edges saved to CSV files.
[ Info: Computing extrapolations...
[ Info: Membrane force resultant extrapolations saved to n_extrapolations-hard-biased-t3ff-16-4.json.
[ Info: Membrane force resultant extrapolations saved to m_extrapolations-hard-biased-t3ff-16-4.json.
[ Info: Shear force resultant extrapolations saved to q_extrapolations-hard-biased-t3ff-16-4.json.
[ Info: Use --compute=plots to generate plots.
```
Inspect the json files to discover how the resultants converged at the four corners
and the diaphragm edge.

For instance, these are the first few lines of `q_extrapolations-hard-biased-t3ff-16-4.json`:
```
{
  "A": {
    "1": {
      "success": true,
      "q_star": 9.74836811791732,
      "q_star_ci": 35.631803373774126,
      "beta_star": 0.755411607792921,
      "beta_star_ci": 0.28217444426345195,
      "edat": [
        {
```
The data here is for the corner A, component of the shear force number 1.
`q_star` is the estimate of the true solution, `q_star_ci` is the estimate 
of the confidence interval halfwidth such that the solution is
`q_star` +/- `q_star_ci`. `edat` is the complete data from which the extrapolation 
and the estimation of the confidence interval was performed.

Plots of the resultants along the edges can be produced using the command line:
```
julia --project=. driver.jl -s hard -e t3ff -N 16 -c plots 
```
Note that this command can only be issued after CSV data files
have been generated by the command with `-c resultants`
was executed, not before.

*Warning*: A LaTeX installation is required, so an error may occur if
no LaTeX distribution could be detected.

