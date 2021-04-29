# AutoDock-Aurora Makefile

# ------------------------------------------------------
# Compiler paths
VH_COMPILER=g++

AURORA_INC_PATH=/opt/nec/ve/veos/include
AURORA_LIB_PATH=/opt/nec/ve/veos/lib64

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

HOST_FLAGS=-g -Wall -Wcomment #-Werror
IFLAGS=-I$(COMMON_DIR) -I$(HOST_INC_DIR) -I$(WRAPPERVEO_DIR) -I$(AURORA_INC_PATH)
LFLAGS=-L$(AURORA_LIB_PATH)
CFLAGS=$(HOST_FLAGS) $(IFLAGS) $(LFLAGS)

# Device sources
K0_NAME=libhello
K_GA_NAME=libkernel_ga
K_LS_NAME=libkernel_ls

KRNL_GA_LIB=$(BIN_DIR)/$(K_GA_NAME).so

K_NAMES=-DK0=$(K0_NAME) -DK_GA=$(K_GA_NAME)

# Kernel flags
KFLAGS=-DKRNL_LIB_DIRECTORY=$(BIN_DIR) \
			 -DKRNL_SRC_DIRECTORY=$(KRNL_DIR) \
			 -DKCMN_SRC_DIRECTORY=$(KCMN_DIR) \
			 $(K_NAMES)

TARGET := autodock-aurora-2
BIN := $(wildcard $(TARGET)*)

# ------------------------------------------------------
# Configuration
# FDEBUG (full) : enables debugging on both host + device
# LDEBUG (light): enables debugging on host
# RELEASE
CONFIG=RELEASE

HOST_DEBUG_BASIC=-O0 -g3 -Wall
HOST_DEBUG_MEDIUM=$(HOST_DEBUG_BASIC) -DPRINT_ALL_VEO_API -DDOCK_DEBUG
HOST_DEBUG_ALL=$(HOST_DEBUG_MEDIUM)

ifeq ($(CONFIG),FDEBUG)
	OPT_HOST =$(HOST_DEBUG_ALL)
else ifeq ($(CONFIG),LDEBUG)
	OPT_HOST =$(HOST_DEBUG_MEDIUM)
else ifeq ($(CONFIG),RELEASE)
	OPT_HOST =-O3
# Do not use "make CONFIG=PROFILE eval"
# Compiling this way is intended only 
# to work with "make profile"
else ifeq ($(CONFIG),PROFILE)
	OPT_HOST =-O3
else 
	OPT_HOST =-O3
endif
PADDING ?= 0
OPT_HOST += -DPADDING=$(PADDING)

OMP=NO

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
CFLAGS+=-DVERSION=\"$(GIT_VERSION)\" #-pg

# ------------------------------------------------------
# Note that the definition of "newline" contains two blank lines.
# Using $(newline) will expand into only one newline.
define newline


endef

# ---------------------------------
# Setting VE specific env variables
# ---------------------------------
export VE_PROGINF := DETAIL
ENVSET = env VE_PROGINF=DETAIL
ifeq ($(DEBUGVE),YES)
ENVSET += VEORUN_BIN="/opt/nec/ve/bin/gdb /opt/nec/ve/veos/libexec/aveorun"
endif
ifeq ($(TRACE),YES)
ENVSET += VEORUN_BIN=$(BIN_DIR)/veorun_ftrace
POSTRUN = $(FTRACE) -f ftrace.out
else
POSTRUN =
endif


# ------------------------------------------------------
# Compiling host and device codes

all: veodock kernel_ga

kernel_ga:
	make -C device OMP=$(OMP) TRACE=$(TRACE)

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
	$(VH_COMPILER) -std=c++11 \
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

PDB      := 3ce3
NEV      := 2048000
NRUN     := 100
NGEN     := 99999
POPSIZE  := 2048
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
	$(ENVSET) $(BIN_DIR)/$(TARGET) \
	-lsmet sw \
	-lsit 300 -lsrat 100.000000 -smooth 0.500 \
	-nev ${NEV} -ngen $(NGEN) -nrun $(NRUN) -psize $(POPSIZE) \
	-lfile $(EVAL_INPUTS_DIR)/$(PDB)/rand-2.pdbqt \
	-xraylfile $(EVAL_INPUTS_DIR)/$(PDB)/flex-xray.pdbqt \
	-ffile $(EVAL_INPUTS_DIR)/$(PDB)/protein.maps.fld \
	-resnam $(PDB)-"`date +"%Y-%m-%d-%H:%M"`"
	$(POSTRUN)

clean:
	make -C device clean
	rm -f $(BIN_DIR)/* initpop.txt ftrace.out
