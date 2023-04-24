mkdir build
cd build
cmake -DPython3_EXECUTABLE="$PYTHON" -DCMAKE_BUILD_TYPE=Release $SRC_DIR/python
cmake --build . --config Release --parallel $CPU_COUNT
cmake --install .
cd ..
rm -r build
