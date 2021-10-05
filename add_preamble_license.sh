#!/bin/bash
# Copy license preamble to all AutoDock-Aurora source and header files
# if license preamble is not present

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
HOST_SOURCES="$HOST_SOURCE_DIR/calcenergy.cpp $HOST_SOURCE_DIR/getparameters.cpp $HOST_SOURCE_DIR/main.cpp $HOST_SOURCE_DIR/miscellaneous.cpp $HOST_SOURCE_DIR/performdocking.cpp $HOST_SOURCE_DIR/processgrid.cpp $HOST_SOURCE_DIR/processligand.cpp $HOST_SOURCE_DIR/processresult.cpp"

# full list of source files
AUTODOCKAURORA_SOURCE="$COMMON_HEADERS $KRNL_HEADERS $KRNL_SOURCES $HOST_HEADERS $HOST_SOURCES"

# Add license-preamble
# Excluding sources that already have it
for f in $AUTODOCKAURORA_SOURCE; do
	if (grep -q "Copyright (C)" $f); then
		echo "License-preamble was found in $f"
		echo "No license-preamble is added."
	else
		echo "Adding LGPL license-preamble to $f ..."
		cat $LICENSE_PREAMBLE "$f" > "$f.new"
		mv "$f.new" "$f"
		echo "Done!"
	fi
	echo " "
done

