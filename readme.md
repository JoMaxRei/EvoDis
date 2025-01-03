# EvoDis

1. make build dir with `mkdir build`
2. cd into build `cd build`
3. create files with `cmake ..` (on windows: `cmake -DCMAKE_PREFIX_PATH="C:\path\to\vcpkg\packages\gsl_x64-windows" ..`)
    * also add `"C:/path/to/vcpkg/packages/**"` to cpp settings of .vscode
4. build project with `cmake --build .`
5. run project with `./EvoDis`

## On Windows
1. Clone vcpkg `git clone https://github.com/microsoft/vcpkg.git`
2. cd into directory and install cd vcpkg `./bootstrap-vcpkg.bat`
3. install gsl `.\vcpkg.exe install gsl`

# Git

1. Add changes with `git add .`
2. Package changes with `git commit`
3. Push changes with `git push`

# Logging

See https://github.com/abumq/easyloggingpp