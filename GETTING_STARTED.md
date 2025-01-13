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

## Dependencies

Astro Genesis 2.0 relies on several scientific libraries to perform computations and handle data. Ensure the following dependencies are installed:

- GNU Scientific Library (GSL)
- HDF5 (Hierarchical Data Format 5)

1. **Install GNU Scientific Library (GSL)**
   
```sh
sudo apt install libgsl-dev
```

2. **Install HDF5**
HDF5 is used to efficiently store and manage large amounts of scientific data.

```sh
sudo apt install libhdf5-dev
```

## Cloning the Repository

First, clone the repository to your local machine:
   ```sh
   git clone https://github.com/yourusername/astro-genesis-2.0.git
```


Navigate to the project directory:
   ```sh
    cd astro-genesis-2.0
    cd simulation
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
./Astro_Genesis
```


If you encounter any issues during the build process, refer to the Troubleshooting section.
