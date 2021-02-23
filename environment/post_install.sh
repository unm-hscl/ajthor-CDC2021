#!/usr/bin/env bash
set -e

### CVX
#
curl -sL  --create-dirs \
  http://web.cvxr.com/cvx/cvx-a64.tar.gz -o CVX/cvx-a64.tar.gz
tar -xzf CVX/cvx-a64.tar.gz -C CVX/
echo 'Downloaded CVX'

### CORA
#
git clone -n https://github.com/TUMcps/CORA.git CORA
echo 'Downloaded CORA'

### MPT3
#
MIRROR=http://people.ee.ethz.ch/~mpt/tbx/pool/
curl -sL --create-dirs \
  ${MIRROR}mpt3/3_2_1/mpt3-3_2_1.tgz -o MPT3/mpt3.tgz \
  ${MIRROR}mpt3doc/3_0_4/mpt3doc-3_0_4.tgz -o MPT3/mpt3doc.tgz \
  ${MIRROR}cddmex/1_0_1/cddmex_1_0_1_GLNXA64.tgz -o MPT3/cddmex.tgz \
  ${MIRROR}fourier/1_0/fourier_1_0_glnxa64.zip -o MPT3/fourier.zip \
  ${MIRROR}glpkmex/1_0/glpkmex_1_0_glnxa64.zip -o MPT3/glpkmex.zip \
  ${MIRROR}hysdel2/2_0_6/GLNXA64/hysdel.zip -o MPT3/hysdel2.zip \
  ${MIRROR}lcp/1_0_3/lcp_1_0_3_GLNXA64.tgz -o MPT3/lcp.tgz \
  ${MIRROR}sedumi/1_3/sedumi_1_3_glnxa64.zip  -o MPT3/sedumi.zip \
  ${MIRROR}espresso/1_0/espresso_linux.tgz -o MPT3/espresso.tgz

tar -xzf MPT3/mpt3.tgz -C MPT3/
tar -xzf MPT3/mpt3doc.tgz -C MPT3/
tar -xzf MPT3/cddmex.tgz -C MPT3/
unzip -qq MPT3/fourier.zip -d MPT3/
unzip -qq MPT3/glpkmex.zip -d MPT3/
unzip -qq MPT3/hysdel2.zip -d MPT3/
tar -xzf MPT3/lcp.tgz -C MPT3/
unzip -qq MPT3/sedumi.zip -d MPT3/
tar -xzf MPT3/espresso.tgz -C MPT3/
echo 'Downloaded MPT3'

### YALMIP
#
git clone -n https://github.com/yalmip/YALMIP.git YALMIP
echo 'Downloaded YALMIP'

### SReachTools
#
git clone -n https://github.com/unm-hscl/SReachTools.git SReachTools
echo 'Downloaded SReachTools'

### Setup MATLAB environment for SReachTools
#
matlab -nodisplay -nosoftwareopengl -r "setup_dependencies;"
