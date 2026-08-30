# Documentation of Assistance for Project nonFunc

*Author: Dennis Frank*

This repository is part of the project HYBRID ESTIMATION OF NON-FUNCTIONAL PROPERTIES OF SOFTWARE (TIME, POWER AND ENERGY) ON EMBEDDED COMPUTER ARCHITECTURES FOR IMAGE AND VIDEO PROCESSING or in short form also called nonFunc. The assistance for the project shown in this repository is about the installation and use of Chipyard and the simulation of differnt video and image processing algorithms. In this README you can find steps and information about the topics of this repository. **The information can also be found in the documentation but latest changes will be added here first.**

---

## Installation of Chipyard

> ⚠️ **The documentation is still WORK IN PROGRESS, last updated on: 05.06.2026**

The installation of Chipyard is described in the [official documentation](https://chipyard.readthedocs.io/en/latest/Chipyard-Basics/Initial-Repo-Setup.html). Some additional steps were added with the help of ChatGPT. These steps were performed on a system with **Ubuntu 24.04.3 LTS** — a different configuration could require additional steps!

### Step 1: Dependencies

**Scala:**
```bash
curl -O -L https://github.com/chipsalliance/chisel/releases/latest/download/chisel-example.scala

curl -sSLf https://scala-cli.virtuslab.org/get | sh
```

**Temurin:**
```bash
sudo apt install -y wget gpg apt-transport-https

echo "deb https://packages.adoptium.net/artifactory/deb \
$(awk -F= '/^VERSION_CODENAME/{print$2}' /etc/os-release) main" | \
sudo tee /etc/apt/sources.list.d/adoptium.list

wget -qO - https://packages.adoptium.net/artifactory/api/gpg/key/public | \
gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/adoptium.gpg > /dev/null

apt update

apt install temurin-17-jdk
```

**Ubuntu Packages:**
```bash
sudo apt-get install autoconf automake autotools-dev curl libmpc-dev libmpfr-dev libgmp-dev \
  libusb-1.0-0-dev gawk build-essential bison flex texinfo gperf libtool patchutils bc \
  zlib1g-dev device-tree-compiler pkg-config libexpat-dev libfl-dev

sudo apt install libguestfs-tools
```

**Miniconda:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

bash Miniconda3-latest-Linux-x86_64.sh

source ~/miniconda3/etc/profile.d/conda.sh
```

After this, clone the repository and run the installation. **Caution: a lot of disk space is required!**

```bash
git clone https://github.com/ucb-bar/chipyard.git

cd chipyard

./scripts/init-submodules-no-riscv-tools.sh

source ./env.sh

./build-setup.sh riscv-tools
```

---

## Simulation

### RTL Simulation

RTL simulation is done with **Verilator**. This simulates the actual hardware architecture of Rocket Chip — every instruction is executed just like the hardware, working in cycles, making it **cycle-accurate**. The main disadvantage is that simulation takes a lot of time, especially for more complex code or larger software.

First, generate the simulation executable for the Rocket Chip (other configurations are available). Run the following command from the chipyard directory — the executable will be placed in `chipyard/sims/verilator`:

```bash
make -C sims/verilator verilog CONFIG=RocketConfig
```

To simulate your own C/C++ code, compile it for RISC-V to produce an `.elf` file:

```bash
riscv64-unknown-elf-gcc -o <name>.elf <name>.c
```

Then navigate to `sims/verilator` and run the simulation:

```bash
./simulator-chipyard.harness-RocketConfig pk \
    "the-whole-path-to-your-.elf-file-starting-at-home-directory"
```

The `pk` argument stands for **Proxy Kernel** — a tiny OS provided by the RISC-V group that handles basic tasks like managing system calls.

Chipyard also includes ready-to-run example programs for Verilator. See section 2.1.5 of the documentation:
https://chipyard.readthedocs.io/en/stable/Simulation/Software-RTL-Simulation.html

When running the built-in example code, `pk` can be omitted.

> **Note:** It may be necessary to run `cmake .` before `make` inside the `tests` directory of Chipyard.

---

### Functional Simulation

Functional simulation can be done with **Spike** (or QEMU). Unlike RTL simulation, only the RISC-V instruction architecture is simulated — instructions are executed according to the RISC-V specification, verifying that the program is logically and architecturally correct.

```bash
spike pk "the-whole-path-to-your-.elf-file-starting-at-home-directory"
```

If you are already in the directory containing the `.elf` file, you can omit the full path. When running example code from the chipyard `tests` directory with Spike, omit `pk`.

---

### Use of Qemu

An idea was to use Qemu as another option for simulation and find out if it provides any benefits over Spike or Verilator. In general Qemu is a machine emulator and virtualizer. The most common use case is System Emulation where you can emulate  an entire machine with CPU, Memory and more. In this way it is also possible to emulate a RISC-V chip. For our project the main benefit of Qemu is to test the UART communication which is also used for the implementation on the FPGA. This was tested with the Edge Detection Code (can be found in the qemu directory). Moreover, Qemu can be used track executed instructions or measure the amount of instructions but this was not tested yet.

---

## Algorithms in Simulation with Spike

### Fast Fourier Transformation (FFT)

The FFT was implemented using [FFTW](https://www.fftw.org/). Manual installation for RISC-V cross-compilation:

```bash
wget http://www.fftw.org/fftw-3.3.10.tar.gz

tar xzf fftw-3.3.10.tar.gz

cd fftw-3.3.10

export CC=riscv64-unknown-elf-gcc
export AR=riscv64-unknown-elf-ar
export RANLIB=riscv64-unknown-elf-ranlib
export HOST=riscv64-unknown-elf

./configure --host=$HOST --enable-static --disable-shared --prefix=$(pwd)/install

make -j4

make install
```

**Required files** (all in the same directory): `image.py`, `fft_output.txt`, `image_data.h`, `fftw.c`, `gen_freq_pic.c`, and an image (e.g. `test.png`).

**Workflow:**

1. Convert an image to an array of values (you will be prompted for the image filename including extension, e.g. `.png`):
   ```bash
   python3 image.py
   ```
   The output is written to `image_data.h`. It will be overwritten each time you run the script with a new image.

2. Compile the FFT code:
   ```bash
   riscv64-unknown-elf-gcc -static \
       -I./fftw-3.3.10/install/include \
       -L./fftw-3.3.10/install/lib \
       -o fftw.elf fftw.c -lfftw3 -lm
   ```
   > Make sure the conda environment is activated and `source env.sh` has been run in the chipyard directory.

3. Run the simulation:
   ```bash
   spike pk fftw.elf
   ```
   Output is written to `fft_output.txt`.

4. Visualize the frequency output (compile once, then run):
   ```bash
   gen_freq_pic.c -o gen_freq_pic -lfftw3 -lm

   ./gen_freq_pic
   ```

---

### Edge Detection

**Required files** (all in the same directory): `image.py`, `edge.c`, `image_data.h`, `stb_image.h`, `stb_image_write.h`, and an image (e.g. `test.png`).

Download the STB headers:
```bash
wget https://raw.githubusercontent.com/nothings/stb/master/stb_image.h
wget https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h
```

**Workflow:**

1. Convert an image to data (the image **must** use `.jpg` format for the simulation):
   ```bash
   python3 image.py
   ```
   Output is written to `image_data.h` and overwritten on each run.

2. Compile the edge detection code:
   ```bash
   riscv64-unknown-elf-gcc -static -o edge.elf edge.c -lm
   ```
   > Make sure the conda environment is activated and `source env.sh` has been run in the chipyard directory.

3. Run the simulation:
   ```bash
   spike pk edge.elf
   ```
   The result image is automatically saved as `edges.png`.

---

### Optical Flow

- Measures the movement of objects in a scene; approximates motion between two frames.
- Necessary assumptions:
  1. Brightness of an image point remains constant over time.
  2. Displacement and time step are very small.
- Possible use of OpenCV.
- **Lucas-Kanade Method:** Assumes the motion field (optical flow) is constant within a small neighborhood around each pixel.

---

## Algorithms in Simulation with Verilator

### Edge Detection

#### Configuration

To run edge detection in Verilator, the Rocket Chip configuration needs to be modified. A custom configuration was added to `chipyard/generators/chipyard/src/main/scala/config/MemorySystemConfigs.scala`:

```scala
class OwnMemoryRocketConfig extends Config( //ADDED BY ME
  new freechips.rocketchip.subsystem.WithoutTLMonitors ++
  new freechips.rocketchip.subsystem.WithExtMemSize((1<<30) * 1L) ++ // 1GB DRAM
  new freechips.rocketchip.rocket.WithNHugeCores(4) ++
  new freechips.rocketchip.subsystem.WithNMemoryChannels(4) ++
  new chipyard.config.AbstractConfig)
```

This configuration adds external memory with multiple channels and uses 4 cores (vs. the default single core). The first line (`WithoutTLMonitors`) helps speed up the simulation, similar to `FastRTLSimRocketConfig` in `RocketConfigs.scala`.

#### Main Code

The modified code is located in the `verilator_edge_dec` directory and named `edgeBare.c`. Unlike the Spike version, this code avoids system calls like `malloc`, `calloc`, or `free` — since Verilator simulates bare-metal conditions, arrays are used instead. The edge detection logic itself is unchanged.

#### Workflow

All code is located in `chipyard/tests`. Add the following lines to `CMakeLists.txt` under the `#Build` section:

```cmake
add_executable(edgeBare edgeBare.c)
target_link_libraries(edgeBare m)
```

1. Prepare the image data (only `.jpg` files are supported; you will be prompted for the filename, width, and height):
   ```bash
   python3 image.py
   ```
   Output is written to `image_data.h` as a 2D array.

2. Compile the test files:
   ```bash
   cmake .

   make
   ```

3. Run the Verilator simulation from `chipyard/sims/verilator`:
   ```bash
   make run-binary CONFIG=OwnMemoryRocketConfig BINARY=../../tests/edgeBare.riscv \
       LOADMEM=1 TIMEOUT_CYCLES=100000000000
   ```

4. To save the output to a file, append:
   ```bash
   > ../../tests/output.txt
   ```

5. Generate an image from the output:
   ```bash
   gcc generate_image.c -o generate_image

   ./generate_image
   ```

#### Further Code and Annotations

The script `generating.py` can generate simulated image data, saved to `sim_data.h`. Output can currently only be printed to the console — saving to a file is planned for future work.

To generate simulated image data:
```bash
python3 generating.py
```

To visualize the generated data as an image:
```bash
gcc getData.c -o getData

./getData
```

## Software Profiling

### perf

perf is a performance analyzing tool in Linux, which is built upon the kernel infrastructure. In general it is possible to cover hardware level (CPU/ Performance Monitoring Unit) features and also software features (software counters, tracepoints). Since perf is a linux tool, it cannot directly be used on a bare-metal platform but it is possible to use perf on a host system with linux on it. This means concrete numbers like cpu cycles cannot be applied to the FPGA. Nonetheless, perf can be used to analyse the code itself like the amount of instructions, branches or memory accesses.
The command "perf stat" runs the given command and gathers performance counter statistics.
The command executed with code edgeBare.c which was compiled with gcc (more details in the README in the directory software_profiling):

```
perf stat -e instructions,branches,branch-misses ./edge
[INFO] Starting edge detection...
[INFO] Done.

 Performance counter stats for './edge':

         7.690.671      instructions                                                          
           603.126      branches                                                              
            10.158      branch-misses                                                         

       0,002538004 seconds time elapsed

       0,001299000 seconds user
       0,001299000 seconds sys
``` 


Another command is "perf record" which runs the given command and gathers a performance counter profile from it into "perf.data". The information can be displayed with "perf report" afterwards. 

## Datasets
With the help of perf (in this case "perf stat") it was possible to gather data about the algorithms we used before in Spike and/or in Verilator. Therefore we created a script called "create_dataset.py" in the directory software_profiling/datasets/edge and software_profiling/datasets/edge. With this script you can automatically compile the code with differnt images and run it with perf to gather data about instructions, cycles, branches, cache-references and cache-misses (this could be modified, see perf documantation for this). The code is executed multiple times with the same image and in the end the average of the valuesfrom  the differnt runs get calculated. You can find the data about the single runs in edge_dec_data.csv and the calculated average in edge_dec_avg, the same for fftw. It is important to mention that the codes get compiled with gcc and not with riscv64-unknown-elf-gcc, only in this way it is possible to use perf. Every kind of output in the code which was used before is disabled in this case.
You can execute the script with 
```
python3 create_dataset.py
```
If you want to generate the header files for the images add "-header" as flag. You can find images in software_profiling/images and header files get saved in software_profiling/header. 

