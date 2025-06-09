source_dir=/data/home/sydong/work/reproducing/forte
#rm -rf $source_dir/CMakeCache.txt $source_dir/CMakeFiles
psi4_plugin=`psi4 --plugin-compile`
$psi4_plugin -S$source_dir \
	-B$source_dir \
	-DCMAKE_C_COMPILER=icx \
	-DCMAKE_CXX_COMPILER=icpx \
	-Dambit_DIR=$HOME/opt/ambit/share/cmake/ambit \
	-DHDF5_ROOT=$HOME/opt/hdf5 \
	-DCMAKE_BUILD_TYPE=Debug \

cmake --build  $source_dir  -j16
