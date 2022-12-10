# https://www.kaggle.com/wolfram77/puzzlef-graph-openmp
import os
from IPython.display import FileLink
src="graph-openmp"
inp="/kaggle/input/graphs"
out="{}.txt".format(src)
!printf "" > "$out"
display(FileLink(out))
!ulimit -s unlimited && echo ""

# Download program
!rm -rf $src
!git clone https://github.com/puzzlef/$src
!echo ""

# Run
!g++ -std=c++17 -O3 -fopenmp main.cxx
!stdbuf --output=L ./a.out $inp/GAP-road.mtx       2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/GAP-twitter.mtx    2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/GAP-web.mtx        2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/GAP-kron.mtx       2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/GAP-urand.mtx      2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/com-Orkut.mtx      2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/com-Friendster.mtx 2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/arabic-2005.mtx    2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/uk-2005.mtx        2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/webbase-2001.mtx   2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/it-2004.mtx        2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/sk-2005.mtx        2>&1 | tee -a "$out"
!stdbuf --output=L ./a.out $inp/kmer_V1r.mtx       2>&1 | tee -a "$out"
