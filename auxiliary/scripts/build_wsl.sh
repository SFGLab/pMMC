#!/bin/bash
set -e
# cd to repo root regardless of where the script is invoked from
# (this script lives at auxiliary/scripts/ — two levels below root)
cd "$(dirname "$0")/../.."
rm -rf build_wsl
mkdir build_wsl
cd build_wsl
cmake .. -DCUDA_ARCH=89 2>&1
make -j4 2>&1
echo "BUILD_DONE"
ls -la pMMC-*
