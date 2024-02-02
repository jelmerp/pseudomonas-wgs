#!/bin/bash

# Install Pseudofinder
mkdir software && cd software || exit
git clone https://github.com/filip-husnik/pseudofinder.git
cd pseudofinder/modules || exit
micromamba create -y -p /fs/ess/PAS0471/jelmer/conda/pseudofinder -f software/pseudofinder/modules/environment.yml

echo '#!/bin/sh'" \
export PATH=\"$(pwd):"'$PATH'\"" \
export ctl=\"$(pwd)/codeml-2.ctl\"" >> ${CONDA_PREFIX}/etc/conda/activate.d/env_vars.sh

cp ../pseudofinder/pseudofinder.py /fs/ess/PAS0471/jelmer/conda/pseudofinder/bin/
cd ..
cp -r modules /fs/ess/PAS0471/jelmer/conda/pseudofinder/bin/
