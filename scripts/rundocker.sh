#!/usr/bin/env bash

dockerName=registry.dune-project.org/dune-fem/dune-fem-dev:latest

if [ $(uname) = "Linux" ]
then
  xhost +si:localuser:$USER
  docker run -it --rm -v $PWD:/host -v dunepy:/dunepy \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro --device /dev/dri \
    -e userId=$(id -u) -e groupId=$(id -g) $dockerName
  xhost -si:localuser:$USER
elif [ $(uname) = "Darwin" ]
then
  echo "on MAC: for X forwarding remember to run 'ipconfig getifaddr en0)' in separate terminal"
  xhost +si:localuser:$USER
  docker run -it --rm -v $PWD:/host -v dunepy:/dunepy \
    -v /tmp/.X11-unix:/tmp/.X11-unix:ro \
    -e DISPLAY=$(ipconfig getifaddr en0):0 --net=host \
    -e userId=$(id -u) -e groupId=$(id -g) $dockerName
  xhost -si:localuser:$USER
else
  echo "uname=$(uname) not recognized - specific settings for 'Linux' and 'Darwin' available"
  echo "starting docker container without X forwarding."
  docker run -it --rm -v $PWD:/host -v dunepy:/dunepy \
    -e userId=$(id -u) -e groupId=$(id -g) $dockerName
fi
