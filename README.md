# AutoDock-Aurora: AutoDock molecular docking for the NEC SX-Aurora TSUBASA

Ported from OpenCL code of AutoDock. Specifically:

* Host code and Solis-Wets local search from [ocladock-fpga (lga-sdx182)](https://git.esa.informatik.tu-darmstadt.de/docking/ocladock-fpga/-/tree/lga-sdx182)
* ADADELTA local search from [AutoDock-GPU (v1.1)](https://github.com/ccsb-scripps/AutoDock-GPU/tree/v1.1/device)

## Instructions

### Cloning

```bash
git clone --single-branch --branch sx-aurora --recurse-submodules --shallow-submodules https://github.com/esa-tu-darmstadt/AutoDock-Aurora.git

cd Autodock-Aurora

make inputs
```

### Compiling for "Full Debug" mode

```
make CONFIG=FDEBUG
```

* Host binary is compiled with symbols enabled for **gdb**
* Device binary is compiled with no vectorization (`-O0 -g`).

To compile only the device code:

```
make -C device CONFIG=FDEBUG kernel_ga
```

### Compiling for FTRACE

PROGINF is always enabled as it costs no overhead.

For FTRACE pass `TRACE=YES` (uppercase "YES"!) as a make variable.

```
make PDB=1yv3 NRUN=16 TRACE=YES eval
```

### Bit-reproducible results

Setting the make variable `REPRO=YES` will initialize the random seeds
to fixed numbers instead of using the time. This must be used in conjunction
with disabling OpenMP on the SX-Aurora Vector Engine (eg. `OMP=NO`, which
is default). Disabling OpenMP is mandatory, on the VE the cores (i.e. different
LGA runs) share the same random generator and the order in which the cores call
it is undetermined.

When running with OMP disabled, please also set `VE_OMP_NUM_THREADS=1`, for example:

```
env VE_OMP_NUM_THREADS=1 make PDB=1yv3 NRUN=16 REPRO=YES OMP=NO eval
```

### Run VE kernel in debugger

This option is useful for debugging the kernel or finding the approximate or exact
location of exceptions on the VE side. This option and FTRACE output (TRACE=YES)
are mutually exclusive.

```
# Normally compiled VE kernel
make PDB=1yv3 NRUN=8 DEBUGVE=YES OMP=YES eval

# Full debug compile of VE kernel (-g -O0)
make PDB=1yv3 NRUN=8 DEBUGVE=YES OMP=YES CONFIG=FDEBUG eval

# Strict in order execution of VE instructions
env VE_ADVANCEOFF=YES make PDB=1yv3 NRUN=8 DEBUGVE=YES OMP=YES eval
```

When the VE kernel starts, you will get a prompt from the VE gdb. Start by typing
`run`. You may set breakpoints at this point.


### Compiling and evaluating

```
make PDB=1yv3 NRUN=16 eval
```

`make eval` automatically configures all input arguments and set of parameters to launch a docking job.

By default `CONFIG`=RELEASE.

* `PDB` is a molecular input ID 
* `NRUN` is the number of molecular docking runs

### Choosing local search method

```
make TESTLS=ad POPSIZE=2048 PDB=1mzc NRUN=50 TRACE=YES eval
```

`TESTLS` can be either `sw` (Solis-Wets) or `ad` (ADADELTA).

### Compiling with OpenMP + running validation using a single molecule

```
make PDB=1yv3 NRUN=16 OMP=YES eval
```

### Compiling with OpenMP + running validation using five molecules

It runs above command bit for five different molecules, PDB = {`1yv3`, `1ywr`, `1mzc`, `1jyq`, `3er5`}.

```
./val.sh
```


## Running on multiple VEs

Multi-VE runs can be used to increase the performance by spreading the work across multiple VEs.

Specify the VE IDs in the environment variable `VE_NODE_IDS` as a comma separated list. For example,
running on four VEs and using up all cores on them requires:

```
export VE_NODE_IDS=0,1,2,3
```

The environment variable can also be used to run autodock on a particular VE, by default VE#8 is chosen.
For example:
```
export VE_NODE_IDS=6
```


### Multiple processes on one VE

When VEs are using uniform memory access mode (non-NUMA), the default number of OpenMP threads
a process can use is equal to the number of cores, eg. 8 cores on VE20B. By limiting the number
of OpenMP threads one can run multiple processes on one VE, for example:

```
# Running 2 OpenMP threads per process
export VE_OMP_NUM_THREADS=2

# Running 4 processes on VE #1
export VE_NODE_IDS=1,1,1,1
```


### NUMA mode

Vector Engines in NUMA mode by default only allow half of the cores to be in one OpenMP process.
A VE10B or VE20B a maximum of 4 OpenMP threads per process. In order to benefit of the NUMA
effects of reduced memory network/ports conflicts, one can use the multi-process mode.
Simply run two processes on each VE:

```
export VE_NODE_IDS=0,0
```
