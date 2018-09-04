#!/bin/bash
set -euxo pipefail

export STACK="umccrise"

wget https://raw.githubusercontent.com/umccr/workflows/master/$STACK/bootstrap-instance.sh && chmod +x bootstrap-instance.sh && /bin/sh bootstrap-instance.sh
