#!/usr/bin/env bash
src="graph-openmp"
out="$HOME/Logs/$src.log"
ulimit -s unlimited
printf "" > "$out"

# Download program
rm -rf $src
git clone https://github.com/puzzlef/$src
cd $src

# Run
g++ -std=c++17 -O3 -fopenmp main.cxx
stdbuf --output=L ./a.out ~/Data/GAP-road.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-twitter.mtx    2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-web.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-kron.mtx       2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/GAP-urand.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-Orkut.mtx      2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/com-Friendster.mtx 2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/arabic-2005.mtx    2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/uk-2005.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/webbase-2001.mtx   2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/it-2004.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/sk-2005.mtx        2>&1 | tee -a "$out"
stdbuf --output=L ./a.out ~/Data/kmer_V1r.mtx       2>&1 | tee -a "$out"
