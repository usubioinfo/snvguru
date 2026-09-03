#!/bin/bash
set -e
rm -rf _build/html 
sphinx-build . ./_build/html -a 

for name in hcapsulatum mtuberculosis influenzaA influenzaA3; do
    if [ -d "$HOME/pipeline/$name" ]; then
        ln -snf "$HOME/pipeline/$name" "./_build/html/$name"
    elif [ -d "$HOME/pipeline/test/$name" ]; then
        ln -snf "$HOME/pipeline/test/$name" "./_build/html/$name"
    fi
done
