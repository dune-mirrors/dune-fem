#!/usr/bin/env bash

if [ $(uname) -q "Linux" ] ;
  xhost +si:localuser:$USER
  docker run -it --rm -v $PWD:/host -v dunepy:/dunepy \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
    -e userId=$(id -u) -e groupId=$(id -g) \
    registry.dune-project.org/dune-fem/dune-fempy-base:latest
elif [ $(uname) -q "Darwin" ] ;
  echo "on MAC: for X forwarding remember to run 'ipconfig getifaddr en0)' in separate terminal"
  xhost +si:localuser:$USER
  docker run -it --rm -v $PWD:/host -v dunepy:/dunepy \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
    -e DISPLAY=$(ipconfig getifaddr en0):0 --net=host \
    -e userId=$(id -u) -e groupId=$(id -g) \
    registry.dune-project.org/dune-fem/dune-fempy-base:latest
else
  echo "uname=$(uname) not recognized - specific settings for 'Linux' and 'Darwin' available"
  echo "starting docker container without X forwarding."
  docker run -it --rm -v $PWD:/host -v dunepy:/dunepy \
    -e userId=$(id -u) -e groupId=$(id -g) \
    registry.dune-project.org/dune-fem/dune-fempy-base:latest
fi
