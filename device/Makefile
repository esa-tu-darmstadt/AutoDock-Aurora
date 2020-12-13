# AutoDock-Aurora Makefile

# ------------------------------------------------------
# Compiler paths
AURORA_INC_PATH=/opt/nec/ve/veos/include
AURORA_LIB_PATH=/opt/nec/ve/veos/lib64

VE_COMPILER=/opt/nec/ve/bin/ncc

VE_EXEC=/opt/nec/ve/bin/ve_exec
FTRACE=/opt/nec/ve/bin/ftrace

# ------------------------------------------------------
# Project directories
COMMON_DIR=../common
BIN_DIR=../bin

# Device sources
K_OBJS = calc_pc.o energy_ia.o energy_ie.o lga.o libkernel_ga.o perform_ls.o rand_gen.o \
 energy_ia_ser.o energy_ie_ser.o calc_pc_ser.o
K_GA_NAME = libkernel_ga

KRNL_GA_LIB = $(BIN_DIR)/$(K_GA_NAME).so

# ------------------------------------------------------
# Configuration
# FDEBUG (full) : enables debugging on both host + device
# LDEBUG (light): enables debugging on host
# RELEASE
CONFIG=RELEASE

KERNEL_TRACE=-ftrace -DENABLE_TRACE
KERNEL_DEBUG_TRACEBACK=-traceback=verbose
KERNEL_DIAGNOSTIC=-report-all
KERNEL_DEBUG_MSG=-Wall -Wcomment -Werror -Wunknown-pragma
KERNEL_DEBUG_UNUSED=-Wunused-variable -Wunused-parameter -Wunused-but-set-parameter -Wunused-but-set-variable -Wunused-value

KERNEL_DEBUG_BASIC=-g $(KERNEL_DEBUG_MSG) $(KERNEL_DEBUG_UNUSED) $(KERNEL_DIAGNOSTIC) -finline-functions -mvector-packed -mvector-low-precise-divide-function
KERNEL_DEBUG_INFO_AND_NOOPT=-g -O0
KERNEL_DEBUG_COMPLETE=$(KERNEL_DEBUG_BASIC) $(KERNEL_DEBUG_INFO_AND_NOOPT)
KERNEL_DEBUG_TRACE=$(KERNEL_TRACE) $(KERNEL_DEBUG_TRACEBACK)

ifeq ($(CONFIG),FDEBUG)
#	OPT_KRNL =$(KERNEL_DEBUG_COMPLETE) -DPRINT_ALL_KRNL -DPRINT_ALL_PC -DPRINT_ALL_IA  -DPRINT_ALL_IE -DPRINT_ALL_LS
	OPT_KRNL =$(KERNEL_DEBUG_COMPLETE)
else ifeq ($(CONFIG),LDEBUG)
	OPT_KRNL =$(KERNEL_DEBUG_BASIC)
else ifeq ($(CONFIG),RELEASE)
	OPT_KRNL =-O4 $(KERNEL_DEBUG_BASIC)
# Do not use "make CONFIG=PROFILE eval"
# Compiling this way is intended only 
# to work with "make profile"
else ifeq ($(CONFIG),PROFILE)
	OPT_KRNL =-O4 $(KERNEL_DEBUG_BASIC) $(KERNEL_DEBUG_TRACE)
else 
	OPT_KRNL =
endif

OMP=NO
FLAGOMP = 

ifeq ($(OMP),YES)
	OPT_KRNL +=-fopenmp
	FLAGOMP = -fopenmp
endif

ifeq ($(TRACE),YES)
	OPT_KRNL += -ftrace -DENABLE_TRACE
	BUILDREQS = build_veorun_ftrace
else
	BUILDREQS =
endif

# ------------------------------------------------------
# Reproduce results (remove randomness)
REPRO=NO

ifeq ($(REPRO),YES)
	REP =-DREPRO
else
	REP =
endif

# ------------------------------------------------------
# Priting out its git version hash
GIT_VERSION := $(shell git describe --abbrev=40 --dirty --always --tags)
CFLAGS+=-DVERSION=\"$(GIT_VERSION)\"

# ------------------------------------------------------
# Note that the definition of "newline" contains two blank lines.
# Using $(newline) will expand into only one newline.
define newline


endef

# ------------------------------------------------------
# Compiling device code

all: kernel_ga $(BUILDREQS)

ADDLIBS = -lm
ifeq ($(OMP),YES)
ADDLIBS += /opt/nec/ve/nlc/2.1.0/lib/libasl_openmp.a
else
ADDLIBS += /opt/nec/ve/nlc/2.1.0/lib/libasl_sequential.a
endif

INCLUDES = -I$(COMMON_DIR) -I. -I/opt/nec/ve/nlc/2.1.0/include

kernel_ga: clean $(KRNL_GA_LIB)

$(KRNL_GA_LIB): $(K_OBJS)
	$(VE_COMPILER) -fpic -shared -o $@ $(OPT_KRNL) $(K_OBJS) $(ADDLIBS)
	@echo $(newline)

#lga.o: lga.c
#	$(VE_COMPILER) -fpic $(INCLUDES) $(KERNEL_DEBUG_BASIC) $(KERNEL_DEBUG_INFO_AND_NOOPT) -o $@ -c $<

#energy_ie.o: energy_ie.c
#	$(VE_COMPILER) -fpic $(INCLUDES) $(KERNEL_DEBUG_BASIC) $(KERNEL_DEBUG_INFO_AND_NOOPT) -o $@ -c $<

#energy_ia.o: energy_ia.c
#	$(VE_COMPILER) -fpic $(INCLUDES) $(KERNEL_DEBUG_BASIC) -g -O2 -o $@ -c $<

# TODO: better dependencies

%.o: %.c
	$(VE_COMPILER) -fpic $(INCLUDES) $(OPT_KRNL) -o $@ -c $<

build_veorun_ftrace:
	touch _dummy.c
	/opt/nec/ve/bin/ncc -o _dummy.o -c _dummy.c
	/opt/nec/ve/bin/mk_veorun_static -ftrace $(FLAGOMP) -o $(BIN_DIR)/veorun_ftrace _dummy.o
	rm -f _dummy.o _dummy.c


clean:
	rm -f $(KRNL_GA_LIB) *.o *.L $(BIN_DIR)/veorun_ftrace
