#!/bin/bash

./autogen.sh --gnu --with-dune=$HOME/work/Dune/src/ \
	--with-problem-dim=$1 \
	--with-world-dim=$2 \
  --with-grape=$HOME/DuneAdds/modules/GrapeNew \
  --with-alberta=$HOME/DuneAdds/modules/alberta-1.2 \
  --enable-blas \
#  --with-alugrid=$HOME/DuneAdds/modules/alu3dgrid \
#  --nodebug \
#  --optim
#  --with-alberta=$HOME/DuneAdds/modules/alberta-1.2 \
#	--with-blas-lib=/hosts/sauron/graid/robertk/Traube/liblinuxGL_dune/ \
#./autogen.sh --intel --with-dune=/hosts/sauron/graid/robertk/Dune/robertk/ \
#./autogen.sh --intel --with-dune=/tmp/robertk/Schrott  \
  
