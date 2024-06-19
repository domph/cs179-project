# CS 179 Final Project: Markus Lendermann & Dominic Phung

This program is a GUI-based fluid dynamics simulator that allows you to visualize, subject to fluid physics, a collection of particles in an invisible container. The GUI has a control panel where you can customize the simulation, like adding more particles (at specific locations and with specific speeds) or changing the particle size. On devices where a compatible NVIDIA GPU is detected, the program has a "Use GPU" toggle to switch between running the simulation on the CPU vs. GPU. The camera of the scene is modifiable, and the controls to do so are listed in the control panel.

OpenGL 4.1 is required to run this program. 

Demo video:
<video src="assets/cs179-demo.mp4" width="1000" height="679" controls></video>

## Pre-built Binaries

This folder contains the source code to compile and run the program. However, we include pre-built binaries for Windows (x64) and macOS in the `bin` folder for convenience. These binaries have been tested on a variety of devices without any failure; should they for whatever reason not work, compiling from the source directly and running should work instead.

To run the macOS binary, you first have to make it executable via `chmod +x ./cs179-mac`. Then, follow these instructions to run it (since the app isn't registered with Apple): [https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac](https://support.apple.com/guide/mac-help/open-a-mac-app-from-an-unidentified-developer-mh40616/mac). 

For macOS, GPU acceleration is not supported since NVIDIA has long since deprecated the NVIDIA CUDA Toolkit for macOS (and indeed, no modern macOS device has an NVIDIA card).

## Required Toolchain

This project requires the following tools:

- `CMake` (>= 3.20)
- A C compiler**
- A C++ compiler that supports C++17**
- A build system (e.g., `make`)
- For non-macOS: NVIDIA CUDA Toolkit (version depends on the platform)
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


### macOS

#### Tool Setup

The simplest way to get all the necessary tools is to install Homebrew, since that installs the Xcode command line tools. The installation steps for Homebrew are on [https://brew.sh/](https://brew.sh/).

Then, `cmake` can be installed via the command `brew install cmake`. At this point, all the necessary tools should be installed.

#### Compile the program

To compile the program on macOS, run the following commands in the root directory:

```
mkdir build
cd build
cmake -G "Unix Makefiles" ..
make
```

This will generate an executable called `cs179-project` inside the `build` folder.

## Results / Performance Analysis

On devices where GPU mode is available, we can observe the improvements of the GPU version via the frame rate (FPS) increase. 

Here is one data point (with the system running under medium load):
- CPU: Intel i7-10875H
- GPU: NVIDIA GeForce RTX 2060 w/ Max-Q Design
    - Particles: 1500
        - CPU FPS: ~30
        - GPU FPS: ~300
        - Speedup: ~10x
    - Particles: 20000
        - CPU FPS: ~2
        - GPU FPS: ~60
        - Speedup: ~30x

There are a couple things we could improve.

There is space for further micro-optimizations in our physics CUDA kernels. For example, in
`cudaApplyCollisionResponse()`, we could pre-calculate the value of `SHAKE(shake_t)` to avoid GPU
threads from having to do repeated work.

We took particular care in designing the memory layout of the various physics data arrays in GPU memory to ensure that we avoid memory bank conflicts. However, there is still some warp divergence occuring due to the nature of the calculations being performed. Furthermore, we are not guaranteeing cache-alignment of the arrays (as the number of particles is not a multiple of the cache-line size); this is something that can potentially be ensured to further increase performance.

Generally, though, we are ensuring the minimal possible amount of memory operations takes place when we perform physics calculations on the GPU - only the particle position and velocity states need to be continuously sent back to the host at the end of each frame calculation, while all other states "live" in GPU memory. The only time the host and device are fully re-synchronized is when the user toggles between using the GPU or the CPU for rendering.