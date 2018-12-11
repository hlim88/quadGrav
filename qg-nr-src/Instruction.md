## Numerical Realtivity with BO-AMR dirver

This instruction contains to build NR infrastructure for our project. 
We import AMR infrastructure from HAD which is written in Fortran 90

### Dependencies

This code has several dependencise. It has been tested with several different 
architecture include Linux, Mac OSX, and different supercomputer.

For example, if your machines have below module, you are ready to use

```
  1) gdb/8.1                        4) mkl/11.2.0                     7) gsl/1.16
  2) openmpi/3.0_gcc-6              5) compiler_intel/15.0.3          8) fftw/3.3.4
  3) defaultenv/4(default)          6) mpi/mpich-3.1.4_intel-15.0.3
```

Here we use intel compliers with `MPICH` for parallel but it could be replaced by
compatiable `OpenMPI`

Please consult your computing resource people to check whether these dependencies
are available or not.

I also refer Matt's homepage for [HAD installation](http://relativity.phys.lsu.edu/postdocs/matt/software/had-docs/had/docs/index.html). This is quite old documentation but procedure
should be similar

### Build the code

Currently, we have only one project code `hyperQuadG` in `qg-had/src` directory. 
In `src` directory, you can find two addtional directories : `amr` and `hyper`
We hardly recommend that **DO NOT** modify this directory unless you understand
what you are doing.

To complie this code, you need to export paths for lib/complier/flags.
To do that easily, we create `ENV` file in main directory 
`quadGrav/qg-nr-src/qg-had`

Please check example `ENV` file [here](https://github.com/hlim88/quadGrav/blob/master/qg-nr-src/qg-had/.ENV_m7).
You can generate your own enviroment file based on your system.

Note that you also need export project folder name and bin file name. In this work,
we have only one so you can simply add

```
export EQS="hyperQuadG"
export EXECNAME="hquadg"
```

Also, you can add additional project code like above

Once you set `ENV` file correctly, then simply type
```
./.ENV_m7;
make
```

If you didn't see any error during compling (it takes time for first time),
you are done.

### Run several test

Once you successfully build the code, you can find bin executable `hquadg`
in `bin` directory. Parameter files are still tuning but you can try one of
par file in `src/hyperQuadD/par`

## Contact

If you have any question, please contact Hyun Lim (hylim1988@gmail.com) 
