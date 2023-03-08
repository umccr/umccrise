PCGR_VERSION="1.4.0"
PCGR_REPO="https://raw.githubusercontent.com/sigven/pcgr/v${PCGR_VERSION}/conda/env/lock"
PREF=$1 # prefix to install both conda envs

if [[ "$OSTYPE" == "darwin"* ]]; then
    PLATFORM="osx"
else
    PLATFORM="linux"
fi

set -x
mamba create --file ${PCGR_REPO}/pcgr-${PLATFORM}-64.lock --prefix ${PREF}/umccrise_pcgr
mamba create --file ${PCGR_REPO}/pcgrr-${PLATFORM}-64.lock --prefix ${PREF}/umccrise_pcgrr
set +x
