
## Compilation Instructions

### Windows (64-bit)

#### Tool setup
To compile on Windows, you need `CMake` (>= 3.12) and the `MinGW-w64` toolchain—specifically `g++` (a version that supports C++20), `gcc`, and `make`. If you have all these, you can skip to the `Compile program` section. If not, the simplest way to get all of these is via [WinLibs](https://winlibs.com):

1. Navigate to [https://winlibs.com](https://winlibs.com)
2. Scroll to the "Downloads" section
3. Download GCC 14.1.0 (with POSIX threads); note that this supports C++20:
    - Choose the UCRT runtime
    - Choose Win64 (can choose the version that includes LLVM/Clang/LLD/LLDB; it doesn't matter)
    - There's no difference between the 7-zip or zip archive format
4. Extract the file somewhere
5. (Optional) Move the `mingw64` folder that's inside the extracted folder somewhere appropriate (e.g. `C:\mingw64`)
6. Add the `bin` folder inside the `mingw64` folder to the `PATH` environment variable.

To verify that everything is working, run the following commands in command prompt, which should all run correctly:

```
cmake --version
g++ --version
gcc --version
mingw32-make --version
```

#### Compile program
Once the above tools are installed, you can compile this program on Windows as follows. In the root folder directory, run the following commands in command prompt:

```cmd
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

This should create an executable called `cs179-project.exe` inside the `build` folder. Note that if you have `make` instead of `mingw32-make`, you can just substitute that instead in the last command above.