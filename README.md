# ScordelisLoRoofBenchmark.jl

The barrel vault (Scordelis-Lo) roof benchmark implementation.

To compute just the deflection at the point B, run the following:
```
ScordelisLoRoofBenchmark.deflection(support="hard", element=:q4rs, mesh="biased")
```

Run a sequence of simulations to compute the resultants.
The results are saved in CSV files, for each mesh resolution and for each edge one.
```
ScordelisLoRoofBenchmark.resultants(support="hard", element=:q4rs, mesh="biased")
```


The CSV files may then be used for the plotting of the resultants along the four edges.
```
ScordelisLoRoofBenchmark.plot_resultants("scolo_q4rs-hard-biased-1024", "n")
ScordelisLoRoofBenchmark.plot_resultants("scolo_q4rs-hard-biased-1024", "m")
ScordelisLoRoofBenchmark.plot_resultants("scolo_q4rs-hard-biased-1024", "q")
```


Finally, the CSV files may then be used for the extraction of the quantities at the corners
(and the intermediate points along the diaphragm) to plot

```
julia> ScordelisLoRoofBenchmark.extrapolate_resultants("scolo_q4rs-hard-biased", ns, "n")
```