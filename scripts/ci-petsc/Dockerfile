FROM duneci/base:9
MAINTAINER "Martin Nolte"

USER root

RUN apt-get update && apt-get dist-upgrade --yes --no-install-recommends
RUN apt-get install --yes --no-install-recommends \
  petsc-dev \
  libeigen3-dev

USER duneci
