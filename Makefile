# AutoDock-Aurora Makefile

# ------------------------------------------------------
# Compiler paths
AURORA_INC_PATH=/opt/nec/ve/veos/include
AURORA_LIB_PATH=/opt/nec/ve/veos/lib64

VH_COMPILER=g++
VE_COMPILER=/opt/nec/ve/bin/ncc

VE_EXEC=/opt/nec/ve/bin/ve_exec
FTRACE=/opt/nec/ve/bin/ftrace

# ------------------------------------------------------
# Project directories
COMMON_DIR=./common
WRAPPERVEO_DIR=./wrapveo
HOST_INC_DIR=./host/inc
HOST_SRC_DIR=./host/src
KRNL_DIR=./device
KCMN_DIR=$(COMMON_DIR)
BIN_DIR=./bin

# Host sources
WRAPPERVEO_SRC=$(wildcard $(WRAPPERVEO_DIR)/*.cpp)
HOST_SRC=$(wildcard $(HOST_SRC_DIR)/*.cpp)
SRC=$(WRAPPERVEO_SRC) $(HOST_SRC)

IFLAGS=-I$(COMMON_DIR) -I$(HOST_INC_DIR) -I$(WRAPPERVEO_DIR) -I$(AURORA_INC_PATH)
LFLAGS=-L$(AURORA_LIB_PATH)
CFLAGS=$(IFLAGS) $(LFLAGS)

# Device sources
K0_NAME=libhello
K1_NAME=libkernel1
K2_NAME=libkernel2
K3_NAME=libkernel3
K4_NAME=libkernel4
K7_NAME=libkernel7

KRNL0_SRC=$(KRNL_DIR)/$(K0_NAME).c
KRNL0_LIB=$(BIN_DIR)/$(K0_NAME).so

KRNL1_SRC=$(KRNL_DIR)/$(K1_NAME).c
KRNL1_LIB=$(BIN_DIR)/$(K1_NAME).so

KRNL2_SRC=$(KRNL_DIR)/$(K2_NAME).c
KRNL2_LIB=$(BIN_DIR)/$(K2_NAME).so

KRNL3_SRC=$(KRNL_DIR)/$(K3_NAME).c
KRNL3_LIB=$(BIN_DIR)/$(K3_NAME).so

KRNL4_SRC=$(KRNL_DIR)/$(K4_NAME).c
KRNL4_LIB=$(BIN_DIR)/$(K4_NAME).so

KRNL7_SRC=$(KRNL_DIR)/$(K7_NAME).c
KRNL7_LIB=$(BIN_DIR)/$(K7_NAME).so

K_NAMES=-DK0=$(K0_NAME) \
        -DK1=$(K1_NAME) \
				-DK2=$(K2_NAME) \
				-DK3=$(K3_NAME) \
				-DK4=$(K4_NAME) \
				-DK5=$(K5_NAME) \
				-DK6=$(K6_NAME) \
				-DK7=$(K7_NAME)

# Kernel flags
KFLAGS=-DKRNL_LIB_DIRECTORY=$(BIN_DIR) \
			 -DKRNL_SRC_DIRECTORY=$(KRNL_DIR) \
			 -DKCMN_SRC_DIRECTORY=$(KCMN_DIR) \
			 $(K_NAMES)

TARGET := autodock-aurora
BIN := $(wildcard $(TARGET)*)

# ------------------------------------------------------
# Configuration
# FDEBUG (full) : enables debugging on both host + device
# LDEBUG (light): enables debugging on host
# RELEASE
CONFIG=RELEASE

HOST_DEBUG_BASIC=-O0 -g3 -Wall
HOST_DEBUG_MEDIUM=$(HOST_DEBUG_BASIC) -DPRINT_ALL_VEO_API
HOST_DEBUG_ALL=$(HOST_DEBUG_MEDIUM) -DDOCK_DEBUG

KERNEL_TRACE=-ftrace -DENABLE_TRACE
KERNEL_DEBUG_TRACEBACK=-traceback=verbose
KERNEL_DIAGNOSTIC=-report-diagnostics -report-format
KERNEL_DEBUG_MSG=-Wall -Wcomment -Werror -Wunknown-pragma
KERNEL_DEBUG_UNUSED=-Wunused-variable -Wunused-parameter -Wunused-but-set-parameter -Wunused-but-set-variable -Wunused-value

KERNEL_DEBUG_BASIC=$(KERNEL_DEBUG_MSG) $(KERNEL_DEBUG_UNUSED) -finline-functions -mvector-low-precise-divide-function
KERNEL_DEBUG_COMPLETE=$(KERNEL_DEBUG_BASIC) $(KERNEL_DIAGNOSTIC)
KERNEL_DEBUG_TRACE=$(KERNEL_TRACE) $(KERNEL_DEBUG_TRACEBACK) $(KERNEL_DIAGNOSTIC) 

ifeq ($(CONFIG),FDEBUG)
	OPT_HOST =$(HOST_DEBUG_ALL)
#	OPT_KRNL =$(KERNEL_DEBUG_BASIC) $(KERNEL_DEBUG_TRACE) -DPRINT_ALL_KRNL -DDEBUG_ENERGY_KERNEL -DDEBUG_ENERGY_KERNEL1 -lm
	OPT_KRNL =$(KERNEL_DEBUG_COMPLETE) -lm
else ifeq ($(CONFIG),LDEBUG)
	OPT_HOST =$(HOST_DEBUG_MEDIUM)
	OPT_KRNL =$(KERNEL_DEBUG_BASIC) -lm
else ifeq ($(CONFIG),RELEASE)
	OPT_HOST =-O3
	OPT_KRNL =$(KERNEL_DEBUG_BASIC) -lm
# Do not use "make CONFIG=PROFILE eval"
# Compiling this way is intended only 
# to work with "make profile"
else ifeq ($(CONFIG),PROFILE)
	OPT_HOST =-O3
	OPT_KRNL =$(KERNEL_DEBUG_BASIC) $(KERNEL_DEBUG_TRACE)  -lm
else 
	OPT_HOST =-O3
	OPT_KRNL =-lm
endif

OMP=NO

ifeq ($(OMP),YES)
	OPT_KRNL +=-fopenmp
else
	OPT_KRNL +=
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
# Compiling host and device codes

all: veodock kernel0 kernel1 kernel2 kernel3 kernel4 kernel7

# Notice: kernel0 is compiled without -shared
# otherwise PROGINFO (via "make profile")
# does not work!
kernel0: $(KRNL0_SRC)
	$(VE_COMPILER) -fpic  -I$(COMMON_DIR) -I$(KRNL_DIR) -o $(KRNL0_LIB) $(KRNL0_SRC) $(OPT_KRNL)
	@echo $(newline)

kernel1: $(KRNL1_SRC)
	$(VE_COMPILER) -fpic -shared -I$(COMMON_DIR) -o $(KRNL1_LIB) $(KRNL1_SRC) $(OPT_KRNL)
	@echo $(newline)

kernel2: $(KRNL2_SRC)
	$(VE_COMPILER) -fpic -shared -I$(COMMON_DIR) -o $(KRNL2_LIB) $(KRNL2_SRC) $(OPT_KRNL)
	@echo $(newline)

kernel3: $(KRNL3_SRC)
	$(VE_COMPILER) -fpic -shared -I$(COMMON_DIR) -o $(KRNL3_LIB) $(KRNL3_SRC) $(OPT_KRNL)
	@echo $(newline)

kernel4: $(KRNL4_SRC)
	$(VE_COMPILER) -fpic -shared -I$(COMMON_DIR) -o $(KRNL4_LIB) $(KRNL4_SRC) $(OPT_KRNL)
	@echo $(newline)

kernel7: $(KRNL7_SRC)
	$(VE_COMPILER) -fpic -shared -I$(COMMON_DIR) -o $(KRNL7_LIB) $(KRNL7_SRC) $(OPT_KRNL)
	@echo $(newline)

export VE_PROGINF := DETAIL

displayproginf:
	@echo " "
	@echo "-------------------"
	@echo "PROGINF"
	@echo "-------------------"
	@echo "VE_PROGINF is set to $$VE_PROGINF"
	@echo " "
	${VE_EXEC} $(KRNL0_LIB)

displayftrace:
	@echo " "
	${FTRACE} -f ./ftrace.out

# Requires make CONFIG=PROFILE
profile: displayproginf displayftrace

# ------------------------------------------------------
veodock: $(SRC)
	make clean
	@echo $(newline)
	$(VH_COMPILER) \
	-o$(BIN_DIR)/$(TARGET) \
	$(CFLAGS) \
	$(SRC) \
	$(OPT_HOST) $(REP) $(KFLAGS) \
	-Wl,-rpath=$(AURORA_LIB_PATH) -lveo
	@echo $(newline)

# ------------------------------------------------------
# Notice that gzip -k option is not supported in
# the workstation hosting Aurora.
# Thus, the command with such option is commented out

EVAL_INPUTS_DIR=./AD-GPU_set_of_42/data
inputs:
	git submodule update --init --recursive
	for dir in $(EVAL_INPUTS_DIR)/* ; do (cd $$dir && gunzip protein.*.map.gz); done
	#for dir in $(EVAL_INPUTS_DIR)/* ; do (cd $$dir && gunzip -k protein.*.map.gz); done

PDB      := 3ce3
NEV      := 2048000
NRUN     := 100
NGEN     := 99999
POPSIZE  := 150
GFPOP    := 1
TESTNAME := test
TESTLS   := sw

test: all
	$(BIN_DIR)/$(TARGET) \
	-ffile ./input/$(PDB)/derived/$(PDB)_protein.maps.fld \
	-lfile ./input/$(PDB)/derived/$(PDB)_ligand.pdbqt \
	-nev ${NEV} \
	-nrun $(NRUN) \
	-ngen $(NGEN) \
	-psize $(POPSIZE) \
	-resnam $(TESTNAME) \
	-gfpop $(GFPOP) \
	-lsmet $(TESTLS)

eval: all
	$(BIN_DIR)/$(TARGET) \
	-lsmet sw \
	-lsit 300 -lsrat 100.000000 -smooth 0.500 \
	-nev ${NEV} -ngen $(NGEN) -nrun $(NRUN) -psize $(POPSIZE) \
	-lfile $(EVAL_INPUTS_DIR)/$(PDB)/rand-2.pdbqt \
	-xraylfile $(EVAL_INPUTS_DIR)/$(PDB)/flex-xray.pdbqt \
	-ffile $(EVAL_INPUTS_DIR)/$(PDB)/protein.maps.fld \
	-resnam $(PDB)-"`date +"%Y-%m-%d-%H:%M"`"

clean:
	rm -f $(BIN_DIR)/* *.L initpop.txt ftrace.out