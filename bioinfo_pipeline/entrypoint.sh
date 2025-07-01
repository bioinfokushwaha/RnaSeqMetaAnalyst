#!/bin/bash
# entrypoint.sh

# Activate conda env and run your command
exec conda run --no-capture-output -n 16_ppl_bioinfo-env "$@"
