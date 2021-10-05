#!/bin/bash
# Remove license preamble from all AutoDock-Aurora source and header files
# if license preamble is present

# Reference
# sed '/PATTERNSTART/,/PATTERNEND/{0,/PATTERNEND/d}' file

# license preamble
LICENSE_PREAMBLE="./PREAMBLE_LICENSE"

# common-header files
COMMON_HEADER_DIR="./common"
COMMON_HEADERS="$COMMON_HEADER_DIR/*.h"

# kernel-source files
KRNL_SOURCE_DIR="./device"
KRNL_HEADERS="$KRNL_SOURCE_DIR/*.h"
KRNL_SOURCES="$KRNL_SOURCE_DIR/*.c"

# host-header files
# excluding "correct_grad_axisangle.h" as its license is taken from AutoDock-GPU
HOST_HEADER_DIR="./host/inc"
HOST_HEADERS="$HOST_HEADER_DIR/calcenergy.h $HOST_HEADER_DIR/ext_headers.h $HOST_HEADER_DIR/getparameters.h $HOST_HEADER_DIR/miscellaneous.h $HOST_HEADER_DIR/packbuff.hpp $HOST_HEADER_DIR/performdocking.h $HOST_HEADER_DIR/processgrid.h $HOST_HEADER_DIR/processligand.h $HOST_HEADER_DIR/processresult.h"

# host-source files
HOST_SOURCE_DIR="./host/src"
HOST_SOURCES="$HOST_SOURCE_DIR/*.cpp"

# full list of source files
AUTODOCKAURORA_SOURCE="$COMMON_HEADERS $KRNL_HEADERS $KRNL_SOURCES $HOST_HEADERS $HOST_SOURCES"

# Print variables
#echo $COMMON_HEADER_DIR/*.h
#echo $KRNL_HEADER_DIR/*.h
#echo $KRNL_SOURCE_DIR/*.c
#echo $HOST_HEADER_DIR/*.h*
#echo $HOST_SOURCE_DIR/*.cpp
#echo $AUTODOCKAURORA_SOURCE

# Remove license-preamble
# Excluding sources that do not have it
for f in $AUTODOCKAURORA_SOURCE; do
	if (grep -q "Copyright (C)" $f); then
		echo "Removing existing license-preamble from $f"
		sed '/\/\*/,/\*\//{0,/\*\//d}' "$f" > "$f.old"
		mv "$f.old" "$f"
		echo "Done!"
	else
		echo "License-preamble was not found in $f"
		echo "No license-preamble is removed."
	fi
	echo " "
done













































































































































































































































































HOST_SOURCE_DIR="./host/src"
HOST_SOURCES="$HOST_SOURCE_DIR/*.cpp"

# full list of source files
#AUTODOCKGPU_SOURCE="$KRNL_HEADER_DIR/*.h $KRNL_SOURCE_DIR/*.cl $HOST_HEADER_DIR/*.h $HOST_SOURCE_DIR/*.cpp"

AUTODOCKGPU_SOURCE="$KRNL_HEADERS $KRNL_SOURCES $HOST_HEADERS $HOST_SOURCES"

# Print variables
#echo $KRNL_HEADER_DIR/*.h
#echo $KRNL_SOURCE_DIR/*.cl
#echo $HOST_HEADER_DIR/*.h
#echo $HOST_SOURCE_DIR/*.cpp
#echo $AUTODOCKGPU_SOURCE

# Add license-preamble
# Excluding sources that already have it, and
# excluding the automatically-generated ./host/inc/stringify.h
for f in $AUTODOCKGPU_SOURCE; do
	if [ "$f" != "$HOST_HEADER_DIR/stringify.h" ]; then

		if (grep -q "Copyright (C)" $f); then
            		echo "Removing existing license-preamble from $f"
		        awk '/^\/\*/{c++} c!=1; /^ \*\//{c++}' ${f} > tmp.txt
        		cp tmp.txt $f
        		echo "Done!"	
		else
			echo "License-preamble was not found in $f"
			echo "No license-preamble is removed."
		fi
		echo " "
	fi
done

