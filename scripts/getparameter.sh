#!/bin/bash




srcgrep "fem[.]" | grep "Parameter[ ]*::" | sed -e 's/^.*"fem[.]/fem./g' -e 's/".*$//g' | sort | uniq

