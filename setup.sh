#!/bin/bash

DUNEDIR=/hosts/morgoth/raid5/dune/

./autogen.sh --gnu \
  --with-dunecommon=$HOME/dune-common \
  --with-dunegrid=$HOME/dune-grid \
	--with-problem-dim=$1 \
	--with-world-dim=$2 \
  --with-grape=$DUNEDIR/modules/grape \
  --x-includes=/usr/X11R6/include \
  --x-libraries=/usr/X11R6/lib \
  --with-alberta=$DUNEDIR/modules/alberta \
  --with-alugrid=$DUNEDIR/modules/alugrid \
  --with-amiramesh=$DUNEDIR/modules/amiramesh \
  --with-ug=$DUNEDIR/modules/ug \
  --enable-blas \
#  --nodebug \
#  --optim
