# Quantum Diffraction

## set_params.cpp
Creates parameter file. Edit the file to set a filename and choose parameters.

Build command:
```
g++ set_params.cpp -larmadillo -o set_params.exe
```
Run command:
```
./set_params.exe
```

## main.cpp
Runs simulation of Gaussian wave packet traveling through slits. Needs a parameter file from set_params.cpp to run.

Build command:
```
g++ main.cpp src/Quantum_box.cpp -I include -larmadillo -o main.exe
```
Run command:
```
./main.exe params/<input_filename.txt> files/<output_filename.bin> <track deviation [true/false]>
```

Make sure that there is a folder named "files", or don't include it in the input filepath for the parameter file. The output file will be a binary file, make sure to input the .bin suffix. The output file will contain an Armadillo cx_mat containing the wavefuntion $u_{ij}^n$ at all timesteps. If the third command line argument is "true", then a file named "deviation.bin" will be created, containing a cx_vector filled with the deviation of the probability function from 1, i.e. $|1-u^2|$.

## animation.py
Animates the simulated wavefunction. Needs one cx_cube binary file.

## detector.py
Creates plot of three probability functions. Needs three cx_cube binary files.

## plot_deviation.py
Creates plot of the deviation file.

## plot_slice.py
Plots wavefunction at three different timesteps. Choose wether to plot magnitude, real part or imaginary part. Needs one cx_cube binary file.
