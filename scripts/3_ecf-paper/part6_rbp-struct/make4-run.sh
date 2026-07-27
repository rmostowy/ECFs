#!/bin/bash

OUTBASE='/Users/rmostowy/MGG Dropbox/Rafal Mostowy/Projects/Data/ECFs/v4/paper-ecf/part6_rbp-struct'

# Run in the monomer mode
./make3-pymol.sh -r 3000 -o "${OUTBASE}/pymol-out-mono" --chain-mode hide --default-chain A

# Run in the trimer mode
./make3-pymol.sh -r 3000 -o "${OUTBASE}/pymol-out-trimer"
