# Installing Maillard

Welcome to Maillard! This library relies on a powerful stack of deep learning (`torch`, `mace-torch`) and quantum chemistry (`pyscf`, `xtb`, `crest`) tools.
Because these tools rely on highly optimized C++ and Fortran code, installation steps vary depending on your operating system. This guide is designed to work from scratch, even if you are new to command-line tools. Please read the section for your specific operating system carefully.

## 1. Prerequisites & Terminal Setup

### 🪟 Windows Users
Native Windows is not supported due to complex chemical library requirements. You must install Windows Subsystem for Linux (WSL2):
Open PowerShell as Administrator and run: 
```
wsl --install
```
Restart your computer. Open the newly installed "Ubuntu" app from your Start Menu and use this terminal for all following steps.

### 🐧 Linux & 🪟 Windows (WSL2) Users

You will run Maillard natively. You need to install Miniforge (a lightweight package manager):

Open your terminal (or Ubuntu app) and run these commands line-by-line:
```
wget [https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh](https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh)
bash Miniforge3-Linux-x86_64.sh -b -p $HOME/miniforge3
source $HOME/miniforge3/bin/activate
conda init
```


Close your terminal window completely and open a new one to apply the changes.

### 🍎 macOS Users (Apple Silicon: M1/M2/M3)

Do not attempt to run this natively on macOS ARM64. Instead, we use Docker to run a lightweight Linux environment. 

**Quick Start**: from the repository root, use the helper script:

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh bootstrap
./scripts/docker_maillard.sh shell
```

**Setup Detail**:
1. Install OrbStack (Recommended) or Docker Desktop. 
2. `./scripts/docker_maillard.sh up` creates or starts the Linux container named `maillard_validation`.
3. `./scripts/docker_maillard.sh bootstrap` creates the `maillard` conda environment, installs CPU-safe dependencies, applies the `xtbiff` patch, and patches the known PyTorch loader issues in `e3nn` and `mace`.
4. `./scripts/docker_maillard.sh shell` opens an interactive shell in `/workspace` with the `maillard` environment activated.
5. All following commands (Section 2-5) may be run inside that shell, or through `./scripts/docker_maillard.sh run ...`.

## 2. Downloading Maillard & Creating the Environment
First, download the Maillard repository to your computer and enter the directory:
```
git clone [https://github.com/PabloAMC/Maillard.git](https://github.com/PabloAMC/Maillard.git)
cd Maillard
```


Next, we will build the software environment. If we do this blindly, the system will try to download 8GB of unnecessary Nvidia GPU drivers. To prevent this, we install the lightweight CPU version of PyTorch first.Run these commands line-by-line (the last step may take 5-10 minutes to download everything):
```
# 1. Create a blank environment named "maillard"
conda create -n maillard python=3.12 -y

# 2. Activate the environment safely
source $(conda info --base)/etc/profile.d/conda.sh
conda activate maillard

# 3. Pre-install Conda's jaxlib (fixes AVX errors on Mac) and CPU PyTorch (saves 8GB of space)
conda install -y -c conda-forge jax jaxlib
pip install torch --index-url [https://download.pytorch.org/whl/cpu](https://download.pytorch.org/whl/cpu)

# 4. Install the rest of the required chemistry dependencies
conda env update --file environment.yml
```

## 3. The xtbiff Patch (Required)

Due to a syntax error in the source code of the crest chemistry package, the internal docking algorithm will crash during Quantum Cluster Growth (QCG) calculations. We must manually download the standalone xtbiff tool to bypass this bug. Run this block of commands to automatically download, extract, and patch your environment:
```
# Ensure we have the tools needed to download and unzip files
conda install -y -c conda-forge wget xz

# Download and extract the standalone binary
wget [https://github.com/grimme-lab/xtbiff/releases/download/v1.1/xtbiff.tar.xz](https://github.com/grimme-lab/xtbiff/releases/download/v1.1/xtbiff.tar.xz)
tar -xf xtbiff.tar.xz

# Move it to your active environment's binary folder and make it executable
mv xtbiff $CONDA_PREFIX/bin/xtbiff
chmod +x $CONDA_PREFIX/bin/xtbiff

# Clean up the downloaded archive
rm xtbiff.tar.xz
```

## 4. Patching Python Dependencies (Required)

Due to a recent security update in PyTorch, the neural network libraries e3nn and mace will cause your code to crash with an UnpicklingError when they try to load their models. We must apply a small patch to their code to tell PyTorch they are safe to run.

Run these two commands to inject the weights_only=False safety flag into their source code:
```
# Patch the e3nn library
sed -i "s/torch.load(os.path.join(os.path.dirname(__file__), 'constants.pt'))/torch.load(os.path.join(os.path.dirname(__file__), 'constants.pt'), weights_only=False)/g" $CONDA_PREFIX/lib/python3.12/site-packages/e3nn/o3/_wigner.py

# Patch the mace library
sed -i "s/torch.load(f=model_path, map_location=device)/torch.load(f=model_path, map_location=device, weights_only=False)/g" $CONDA_PREFIX/lib/python3.12/site-packages/mace/calculators/mace.py
```


## 5. Verifying the Installation
Let's ensure everything is working by running a test calculation, followed by the library's automated test suite.

A. Running a Quantum Cluster Growth (QCG) Example
We will create a dedicated folder, generate dummy molecules, and run a cluster simulation.
```
# Create a directory for our test and enter it
mkdir -p example_qcg
cd example_qcg

# Create a dummy solute file (test_crest.xyz)
cat <<EOF > test_crest.xyz
3
water
O 0.0 0.0 0.0
H 0.0 0.0 0.96
H 0.0 0.92 -0.25
EOF

# Create the solvent file (water.xyz)
cat <<EOF > water.xyz
3
water
O 0.0 0.0 0.0
H 0.0 0.0 0.9584
H 0.0 0.9240 -0.2529
EOF

# Run the cluster growth calculation using the patched binary
crest test_crest.xyz -qcg water.xyz -nsolv 1 --xtbiff
```

If successful, the terminal will print out optimization steps and conclude with CREST terminated normally. You will see several new output files in this folder containing your optimized molecular coordinates.

B. Running the Library Tests

Finally, navigate back to the main Maillard directory and run the automated test suite to verify the Python library is communicating with the chemistry tools correctly:
```
# Move back to the main Maillard folder
cd ..

# Run the test suite
python -m pytest tests/
```

If the tests pass smoothly, your environment is perfectly configured. Welcome to Maillard!

## 6. Returning to Your Work

If you close your terminal or restart your computer, you do not need to reinstall anything. Follow the steps below for your operating system to jump back in.

### 🍎 macOS Users (Docker/OrbStack)

```bash
./scripts/docker_maillard.sh up
./scripts/docker_maillard.sh shell
```

### 🐧 Linux & 🪟 Windows (WSL2) Users

Open your terminal (or Ubuntu app).

Navigate to your project folder: cd path/to/Maillard

Reactivate the environment:
```
conda activate maillard
```

If the tests pass smoothly, your environment is perfectly configured. Welcome to Maillard!
