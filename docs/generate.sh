#!/bin/bash
rm -r _build/html 
sphinx-build . ./_build/html -a 
ln -s ~/pipeline/hcapsulatum/ ./_build/html/hcapsulatum
ln -s ~/pipeline/mtuberculosis/ ./_build/html/mtuberculosis
ln -s ~/pipeline/influenzaA/ ./_build/html/influenzaA
python -m http.server 9888
