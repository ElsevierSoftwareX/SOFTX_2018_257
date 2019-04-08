#### Requirements

The only software requirement is a C++ compiler.
The default compiler, set in `code/Makefile`, is GNU g++.

The Eigen library is included in the source code distribution, under `code/Eigen`.

#### Automatic build and run:

Go to the `code` folder and execute the `run.sh` script. 
This will build and run the two executables (`tendonmech` and `tendonmech_test`). 
The results are printed at the standard output and also stored in the corresponding text files. 
A detail description is provided in the following two sections.

#### Build:

Please edit `code/Makefile`, if necessary.
You can specify the C++ compiler and flags and enable report of timings.

From the `code` folder of the package, type `make` to compile and build the code.

This will create two executables (`tendonmech` and `tendonmech_test`) in the same folder.
The two executables differ only in the problem size. The first one (`tendonmech`) is used for production runs while the second one is created only for testing purposes.
More specifically, `tendonmech_test` is compiled with `NSLICES=10` and `NRADIANS=11`.
For the production version (`tendonmech`), the code uses the default values for both of these parameters.

#### Run:

To run the executables, just type `./tendonmech` or `./tendonmech_test`, from the `code` folder. 

The executables read the default configuration file (`InputHelixData.txt`), available in the `code` folder of the package.

The contents of `InputHelixData.txt` are as follows:
```
Angle = 60
Center =10
FiberContent = 40
Ef = 2000
Em = 0.25
```

Users can specify a configuration file as runtime argument to the `tendonmech` executable, for instance `./tendonmech myconfig.txt`

#### File Description:

- `code/*`: Source files, Makefile and License file of the TendonMech program
- `code/InputHelixData.txt`: Default configuration file
- `code/sample_test_results.txt`: reference output results of `tendonmech_test` (printer to standard output)
- `data/*`: Data (text) files with position of fiber positions
- `results/*`: Generated text files with the standard output messages of the `tendomech_test` and `tendonmech` executables
