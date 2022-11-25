#!/usr/bin/env bash

rsync -aP --exclude 'build' . calcc:/mnt/localstorage/fast/$USER/QNOSE_QuantumSimulation
