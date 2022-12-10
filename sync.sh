#!/usr/bin/env bash

# rsync -aP --exclude 'build' --exclude '.git' --exclude 'venv' . calcc:/mnt/localstorage/fast/$USER/QNOSE_QuantumSimulation
rsync -aP --exclude 'build' --exclude '.git' --exclude 'venv' . lupus:git/QuantumSimulation
