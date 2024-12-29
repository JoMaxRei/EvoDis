# EvoDis

1. make build dir with `mkdir build`
2. cd into build `cd build`
3. create files with `cmake ..` (on windows: `cmake -DCMAKE_PREFIX_PATH="C:/msys64/mingw64" ..`)
    * also add `"C:/msys64/mingw64/**"` to cpp settings of .vscode
4. build project with `cmake --build .`
5. run project with `./EvoDis`

## On Windows
1. Download and install MSYS2.
2. Update MSYS2 packages: `pacman -Syu`
3. Install GSL using MSYS2: `pacman -S mingw-w64-x86_64-gsl`
This installs GSL for the 64-bit MinGW toolchain.

# Git

1. Add changes with `git add .`
2. Package changes with `git commit`
3. Push changes with `git push`