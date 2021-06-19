#!/bin/bash
git clone https://github.com/fanzheng10/CliXO-1.0.git
cd CliXO-1.0
# git checkout 13ab67b2e9e3f3077bca20542412cf6df81113fc
if [ -x clixo ]; then
  echo "CliXO successfully installed!"
else
  echo "Failed to install CliXO. Please open an issue on https://github.com/fanzheng10/CliXO-1.0"
fi