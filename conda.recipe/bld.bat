mkdir build
cd build
cmake -DPython3_EXECUTABLE="%PYTHON%" -DCMAKE_BUILD_TYPE=Release %RECIPE_DIR%\..\python
cmake --build . --config Release --parallel %CPU_COUNT%
cmake --install .
cd ..
rmdir /S /Q build

if errorlevel 1 exit 1