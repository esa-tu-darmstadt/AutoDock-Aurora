# AutoDock-Aurora-2: AutoDock running on NEC SX-Aurora

Ported from OpenCL code of []().

## Original instructions

[Original README of AutoDock-GPU](./README-ORIGINAL.md).

## SX-Aurora specific instructions

### Cloning

**FIXME!**

```bash
git clone --single-branch --branch sx-aurora --recurse-submodules --shallow-submodules https://gitlab.com/postdoc_tud/molecular-docking/autodock-aurora/autodock-aurora-2.git

cd autodock-aurora

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


### Run VE Kernel in debugger

This option is useful for debugging the kernel or finding the approximate or exact
location of exceptions on the VE side. This option and FTRACE output (TRACE=YES)
are mutually exclusive.

```
# normally compiled VE kernel
make PDB=1yv3 NRUN=8 DEBUGVE=YES OMP=YES eval

# full debug compile of VE kernel (-g -O0)
make PDB=1yv3 NRUN=8 DEBUGVE=YES OMP=YES CONFIG=FDEBUG eval

# strict in order execution of VE instructions
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

### Compiling with OpenMP + running validation using a single molecule

```
make PDB=1yv3 NRUN=16 OMP=YES eval
```

## Compiling with OpenMP + running validation using five molecules

It runs above command bit for five different molecules, PDB = {`1yv3`, `1ywr`, `1mzc`, `1jyq`, `3er5`}.

```
./val.sh
```


