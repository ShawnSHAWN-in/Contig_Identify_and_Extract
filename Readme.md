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
## Properity & Support
* Mutithreads, defalut cores > 15
* bitvector, intensive memory storage
* Out put nucleotide or amino acid fasta file, with rectified frame
## Compile Limits
> Cmake >= 3.4
> gcc >= 10
## Installtion
```
cmake .
# One thread to compile. 
make .
# Now, IEContig you have getted.
./IEContig
```
## future commitment