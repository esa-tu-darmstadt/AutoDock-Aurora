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

### Compiling for PROGINF and FTRACE

PROGINF is always enabled as it costs no overhead.

For FTRACE pass `TRACE=YES` (uppercase "YES"!) as a make variable.
```
make PDB=1yv3 NRUN=16 TRACE=YES eval
```

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


