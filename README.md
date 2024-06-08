# CS 179 Final Project: Markus Lendermann & Dominic Phung

This program is a GUI-based fluid dynamics simulator (as such, it cannot be used entirely from the command-line but instead requires a display). OpenGL 4.1 is required to run this program.

<img src="assets/demo.png" width="1000"/>

## CUDA Acceleration
Comments detailing how this program will be parallelized to run on a GPU with CUDA can be found in `physics.cpp`.

## Pre-built Binaries

This folder contains the source code to compile and run the program. However, we include pre-built binaries for Windows (x64) and MacOS in the `bin` folder for convenience. These binaries have been tested on a variety of devices without any failure; should they for whatever reason not work, compiling from the source directly and running should work instead.

To run the MacOS binary, you first have to make it executable via `chmod +x ./cs179-project-mac`. Then, follow these instructions to run it (since the app isn't registered with Apple): [https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac](https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac).

## Required Toolchain

This project requires the following tools:

- `CMake` (>= 3.12)
- A C compiler
- A C++ compiler that supports C++17
- A build system (e.g., `make`)

Platform-specific instructions for setting up the required toolchain are provided below. 

## Platform-Specific Toolchain Setup & Compilation Instructions

### Windows (64-bit)

#### Tool setup
The simplest way to set up the toolchain on Windows is by using the `MinGW-w64` toolchain. Follow these steps:

1. Navigate to [https://winlibs.com](https://winlibs.com)
2. Scroll to the "Downloads" section
3. Download GCC 14.1.0 (with POSIX threads), which supports C++20:
    - Choose the UCRT runtime
    - Choose Win64 (either version with or without LLVM/Clang/LLD/LLDB)
    - Choose either the 7-zip or zip archive format
4. Extract the downloaded file
5. (Optional) Move the `mingw64` folder inside the extracted file to a preferred location (e.g. `C:\mingw64`)
6. Add the `bin` folder inside the `mingw64` directory to the `PATH` environment variable

To verify the setup, run the following commands in the Command Prompt:
```
cmake --version
g++ --version
gcc --version
mingw32-make --version
```

If these commands fail, ensure the `PATH` environment variable includes the `mingw64\bin` directory and restart your terminal or VSCode to apply the changes.

#### Compile the program
To compile the program on Windows, run the following commands in the root directory using Command Prompt:

```
mkdir build
cd build
cmake -G "MinGW Makefiles" ..
mingw32-make
```

This will generate an executable called `cs179-project.exe` inside the `build` folder.


### MacOS

#### Tool Setup

The simplest way to get all the necessary tools is to install Homebrew, since that installs the Xcode command line tools. The installation steps for Homebrew are on [https://brew.sh/](https://brew.sh/).

Then, `cmake` can be installed via the command `brew install cmake`. At this point, all the necessary tools should be installed.

#### Compile the program

To compile the program on MacOS, run the following commands in the root directory:

```
mkdir build
cd build
cmake -G "Unix Makefiles" ..
make
```

This will generate an executable called `cs179-project` inside the `build` folder.