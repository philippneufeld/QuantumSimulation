#!/usr/bin/env bash

# rsync -aP --exclude 'doc' --exclude 'build' --exclude '.git' --exclude 'venv' . calcc:/mnt/localstorage/fast/$USER/QNOSE_QuantumSimulation
rsync -aP --exclude 'doc' --exclude 'build' --exclude '.git' --exclude 'venv' . ~/remote_home/calc/QNOSE_QuantumSimulation
# rsync -aP --exclude 'doc' --exclude 'build' --exclude '.git' --exclude 'venv' . lupus:git/QuantumSimulation
