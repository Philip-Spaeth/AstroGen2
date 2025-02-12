# Getting Started with Astro Genesis 2.0

This guide will help you set up and compile Astro Genesis 2.0 using CMake and GCC on Linux.

## Prerequisites

Before you begin, ensure you have the following software installed on your system:

- CMake (version 3.10 or higher)
- GCC (GNU Compiler Collection)

1. **Install CMake and GCC**:
```sh
sudo apt update
sudo apt install cmake gcc g++
```

## Cloning the Repository

First, clone the repository to your local machine:
   ```sh
   git clone https://github.com/Philip-Spaeth/AstroGen2
```


Navigate to the project directory:
   ```sh
    cd AstroGen2/simulation
```

Build the project using CMake:
   ```sh
   mkdir build
   cd build
   cmake ..
   make
```

## Running the Compiled Program

After successfully building the project, you can run the program.
```
./AstroGen2
```


If you encounter any issues during the build process, refer to the Troubleshooting section.
