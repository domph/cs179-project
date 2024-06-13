# CS 179 Final Project: Markus Lendermann & Dominic Phung

This program is a GUI-based fluid dynamics simulator that provides the option of GPU acceleration if a compatible NVIDIA GPU is detected. OpenGL 4.1 is required to run this program.

<img src="assets/demo.png" width="1000"/>

## Pre-built Binaries

This folder contains the source code to compile and run the program. However, we include pre-built binaries for Windows (x64) and MacOS in the `bin` folder for convenience. These binaries have been tested on a variety of devices without any failure; should they for whatever reason not work, compiling from the source directly and running should work instead.

To run the MacOS binary, you first have to make it executable via `chmod +x ./cs179-project-mac`. Then, follow these instructions to run it (since the app isn't registered with Apple): [https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac](https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac).

## Required Toolchain

This project requires the following tools:

- `CMake` (>= 3.20)
- A C compiler**
- A C++ compiler that supports C++17**
- A build system (e.g., `make`)
- NVIDIA CUDA Toolkit (version depends on the platform)

Platform-specific instructions for setting up the required toolchain are provided below. 

** Only compilers supported by CUDA work (these are detailed below for each platform)

## Platform-Specific Toolchain Setup & Compilation Instructions

### Windows (64-bit)

#### Tool setup

On Windows, `MSVC` is the only suppported compiler. The simplest way to set up the toolchain is by installing Visual Studio (Community) 2022 via [https://visualstudio.microsoft.com/vs/](https://visualstudio.microsoft.com/vs/). In the installer page, make sure `Desktop development with C++` is selected under the `Workloads` section. This will install all the necessary tools for compiling this project. If you already have Visual Studio installed and are unsure if you have the right configuration, you can open up `Visual Studio Installer` and click `Modify` under the `Visual Studio (Community) 2022` section, and unselect (if it was already selected) and reselect the `Desktop development with C++` option.

Visual Studio does not update the `PATH` environment variable globally. Thus, by default, commands cannot be run in the regular `Command prompt`; instead, commands need to be run from `x64 Native Tools Command Prompt for VS 2022`. This specialized command prompt updates the `PATH` environment for the session appropriately. It can be found by the Search feature of Windows or by navigating to the `Start Menu > Visual Studio 2022` folder.

Finally, NVIDIA CUDA Toolkit 12.5 needs to be installed. The installed can be found here: [https://developer.nvidia.com/cuda-downloads](https://developer.nvidia.com/cuda-downloads).

To verify that the setup is working, run the following commands in the `x64 Native Tools Command Prompt`:

```
where cmake
where nmake
where cl
nvcc --version
```

The first 3 commands should all print out paths that start with `C:\Program Files\Microsoft Visual Studio\2022\Community\`, and the final command should print a release number of 12.5.

#### Compile the program
To compile the program on Windows, run the following commands in the root directory of the project using `x64 Native Tools Command Prompt for VS 2022`:

```
mkdir build
cd build
cmake -G "NMake Makefiles" ..
nmake
```

This will generate an executable called `cs179-project.exe` inside the `build` folder.


### MacOS -- WIP

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