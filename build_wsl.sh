#!/bin/bash
set -e
cd /mnt/d/git/pMMC
rm -rf build_wsl
mkdir build_wsl
cd build_wsl
cmake .. -DCUDA_ARCH=89 2>&1
make -j4 2>&1
echo "BUILD_DONE"
ls -la pMMC-*
