# Set SDAccel up

```zsh
source /opt/cad/xilinx/sdx/SDx/2017.4/settings64.sh
```

# Documentation

A detailed explanation can be found in:

_"SDAccel Environment Tutorial: Introduction (UG1021)"_

Chapter 2 (_Flow Overview_) -> Lab 2 (_Introduction to the SDAccel Makefile_)

# Important commands

* See basic rules under [Makefile](./ofdock_taskpar_xl/Makefile)
* See user compilation flags under [boards.mk](./common_xilinx/utility/boards.mk)
* See OpenCL ICD options under [opencl.mk](./common_xilinx/libs/opencl/opencl.mk)

## SW and HW Emulation

```zsh
make swemu
```

```zsh
make hwemu
```

## FPGA building

```zsh
make hw
```

## Where to find reports

A Makefile flow is used. A system estimate can be generated this way:

```zsh
make kerrpt-sw
```

In order to generate additional reports, a [`sdaccel_ini`](./ofdock_taskpar_xl/sdaccel.ini) was created.
