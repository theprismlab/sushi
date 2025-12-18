#!/bin/bash

echo Starting replicate collapse...

if [ -z "$BUILD_DIR" ]
then
	echo BUILD_DIR not specified
    exit -1
fi

if [ -z "$LFC" ]
then
	echo LFC parameter empty
    exit -1
fi


#Enforces abs paths
if [[ "$LFC" = /* ]]
then
	LFC=$(ls $LFC)
else
	LFC=$BUILD_DIR/$LFC
fi


echo Build dir is: $BUILD_DIR
echo LFC is: $LFC

echo Rscript collapse_replicates/collapse_replicates.R -c $LFC	\
--out $BUILD_DIR \
--sig_cols $SIG_COLS \
--cell_line_cols $CELL_LINE_COLS \
--collapsed_l2fc_file ${COLLAPSED_LFC}


Rscript collapse_replicates/collapse_replicates.R -c $LFC	\
--out $BUILD_DIR \
--sig_cols $SIG_COLS \
--cell_line_cols $CELL_LINE_COLS \
--collapsed_l2fc_file ${COLLAPSED_LFC}