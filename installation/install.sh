#!/bin/bash
# Install CliXO
git clone https://github.com/fanzheng10/CliXO-1.0.git
cd CliXO-1.0
git checkout 13ab67b2e9e3f3077bca20542412cf6df81113fc
cd ..

# Install DDOT
git clone https://github.com/michaelkyu/ddot.git
cd ddot
make
pip install .