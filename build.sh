#. ~/llvm/clang9.sh
mkdir -p build
cd build

cmake ..                                \
    -DCMAKE_BUILD_TYPE=Release          \
    -GNinja

ninja
