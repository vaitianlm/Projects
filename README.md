# Ising_model

## Ising_model.cpp
Uses Metropolis algorithm to draw energy and magnetisation samples from the 2D-Ising model.

Build command:
```
g++ ising_model.cpp -fopenmp src/utils.cpp -I include -o ising_model.exe -fopenmp
```
Run command:
```
./ising_model.exe <min temperature> <max temperature> <#temp-steps> <lattice length> <initial config [o/u]> <store samples [true/false]> <MCMC cycles> <burn-in cycles>
```
Outputs a file named "L\<lattice length\>_N\<MCMC_cycles\>_\<#temp_steps\>values.txt", containing temperature, mean epsilon, mean m, heat capacity per spin, magnetic susceptility per spin, epsilon squared and m squared for every temperature, in that order. If store samples is set to true and min temperature is equal to max temperature, it will also output a file containing MCMC-cycles done, epsilon, running mean epsilon, heat capacity until current cycle and magnetic susceptibility until current cycle, for every MCMC-cycle done, in that order.

## compare.py
Reads file containing MCMC-cycles done, epsilon, running mean epsilon, heat capacity until current cycle and magnetic susceptibility until current cycle, for every MCMC-cycle done, and plots the wanted quantity and it's analytical solution.

## plot_burn_in.py
Reads four files containing MCMC-cycles done, epsilon, running mean epsilon, heat capacity until current cycle and magnetic susceptibility until current cycle, for every MCMC-cycle done, and plots epsilon and mean epsilon for the four different files in the same plot.

## distribution.py
Plots histogram approximating pdf of epsilon of the $20 \times 20$ Ising-lattice.

## plot_values.py
Reads files containing temperature, mean epsilon, mean m, heat capacity per spin, magnetic susceptility per spin, epsilon squared and m squared for different temperatures, and plots one of the quantities against temperature. Also performs secound order fitting of the data, and plots the fit.

## critical_point.py
Estimates critical temperature at $L \rightarrow \infty$, and makes a plot of it.

