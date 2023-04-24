mkdir build
cd build
cmake -DPython3_EXECUTABLE="%PYTHON%" -DCMAKE_BUILD_TYPE=Release %SRC_DIR%\python || exit /b
cmake --build . --config Release --parallel %CPU_COUNT% || exit /b
cmake --install . || exit /b
cd ..
rmdir /S /Q build

