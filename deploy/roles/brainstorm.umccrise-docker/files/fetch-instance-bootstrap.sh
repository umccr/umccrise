#!/bin/bash
set -euxo pipefail

STACK="umccrise"

SCRIPT="bootstrap-instance.sh"

wget https://raw.githubusercontent.com/umccr/workflows/master/$STACK/$SCRIPT -O "/$SCRIPT" && chmod +x "/$SCRIPT" && /bin/bash "/$SCRIPT"
