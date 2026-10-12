# IEcontig document

Welcome to the documentation of IEContig.
## What's it IEContig
> IEContig is a C++ programmer, with powful running time, aimming to **Identify** and **Extract** conserved gene region from big-data, de-nove Contigs and Consensus file driving sam/bam file.
## what can it do for your projects?
Well question!   
IEContig was initially designed for a basic component in our projects, eg. Conserved sequences Collection of desultorily set for Phylogenetic tree, praparing coding sequence with rectified frame for Codon Usage analysing, CAI, So on.
## How to use it
```
.\IEContig.exe -i <Input file> -t [NT] -I <CD file>  -o <Out file> -y [NT|AA]  -s 0
```
The reference file (`-I`) contains amino acid sequences. Use `-t NT` for
six-frame translation of nucleotide input, with `-y NT` or `-y AA` output.
Use `-t AA -y AA` to extract a matching region directly from protein input.
`-t AA -y NT` is rejected because amino acids do not uniquely determine DNA.
Uppercase and lowercase sequences are accepted. Empty files and records are
reported as errors; nucleotide records shorter than one codon yield no match.

Run the CLI regression tests against a compiled executable:
```
python3 tests/test_cli.py /absolute/path/to/IEContig
```

## Docker

Build from the repository root:
```bash
docker build -t iecontig:local .
docker run --rm iecontig:local --help
```
The multi-stage build uses the same CMake configuration as a native build,
automatically downloads SeqAn3 3.3.0, compiles with GCC/C++20, and runs
the CLI regression tests. The final image contains the executable and runtime
libraries, with no compiler or source dependencies. It runs as a non-root user.

Put your nucleotide input and amino acid reference files in a `data` directory.
On Linux, use your host user ID so the container can write results into it:
```bash
mkdir -p data
docker run --rm --user "$(id -u):$(id -g)" \
  -v "$PWD/data:/data" iecontig:local \
  -i /data/input.fasta -t NT -I /data/reference.fasta \
  -o /data/result.fasta -y AA -s 0
```
The output is saved as `data/result.fasta` on the host. Use `-y NT` for
nucleotide output, or `-t AA -y AA` for amino acid input. On Docker Desktop
for Windows/macOS, omit `--user` if host UID mapping is not applicable, and
use the appropriate absolute host path for the volume mount.

Building needs Internet access to the Docker registry, Debian package mirrors,
and GitHub. Processing local FASTA files does not require a server or database.
If your environment uses an HTTP/HTTPS proxy, forward its existing settings
with `docker build --build-arg HTTP_PROXY --build-arg HTTPS_PROXY --build-arg NO_PROXY -t iecontig:local .`.
For Git access through a proxy with a custom certificate authority, also add
`--secret id=build_ca,src=/path/to/trusted-ca-bundle.pem`. The CA bundle is
mounted only during dependency download; TLS verification stays enabled.

## Properity & Support
* Mutithreads, defalut cores > 15
* bitvector, intensive memory storage
* Out put nucleotide or amino acid fasta file, with rectified frame
## Build with CMake

Requirements: CMake >= 3.24, Make (or Ninja), GCC >= 10 with C++20/OpenMP support, Git,
and zlib development headers. Python 3 is needed only for regression tests.
For example, on Debian/Ubuntu:
```bash
sudo apt-get update
sudo apt-get install cmake make g++ git zlib1g-dev python3
```

Clone and build without installing SeqAn3 manually:
```bash
git clone https://github.com/ShawnSHAWN-in/Contig_Identify_and_Extract.git
cd Contig_Identify_and_Extract
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
./build/IEContig --help
```
CMake FetchContent downloads SeqAn3 3.3.0 and its pinned SDSL/cereal
submodules automatically into the build directory. No user-specific include
or library paths are embedded in the project. The SeqAn3 commit is fixed for
reproducible builds and compatibility with the argument parser used here.
The first configuration needs GitHub access; subsequent builds reuse the
downloaded sources. A source ZIP works too, provided Git is installed.

Enable and run the regression suite:
```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build --parallel 2
ctest --test-dir build --output-on-failure -V
```
CTest runs the Python CLI suite (17 test methods, including six-frame checks).
Tests are off by default so Python is not required just to compile the program.

If SeqAn3 3.3.0 and its dependencies are already installed, opt into package
discovery rather than downloading sources:
```bash
cmake -S . -B build-system -DIECONTIG_USE_SYSTEM_SEQAN3=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-system --parallel 2
```
For installations outside standard locations, use CMake's normal
`CMAKE_PREFIX_PATH` or `seqan3_DIR` options. After a successful initial download,
`-DFETCHCONTENT_FULLY_DISCONNECTED=ON` permits configuring again offline with
the same build directory. Optional installation uses
`cmake --install build --prefix <install-directory>`.
## future commitment
