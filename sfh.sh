#!/bin/bash

cd ~/aa_code/cocohype/cigre-test-grid-in-julia/cigre_model

ssh strenge@cluster.pik-potsdam.de

cd dgunit

sbatch submit.sh
