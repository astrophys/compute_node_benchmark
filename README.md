# Compute Node Benchmarking Tests
## Author : Ali Snedden
## License: MIT
## Purpose
We are seeking "Requests For Proposals" from various vendors of high performanch computing systems.
In order to make a selection from the vendors, we are providing small benchmarking test to be run on a typical This code builds a Singularity image 

## Installation
### Dependencies
1. Singularity (tested with Singularity 2.5.1)
2. A `/tmp` which is used as a bind point for the Singularity image. Most Linux machines should have this. This is used to place data in.
3. A machine that you have root permission on. 

### To Install
Run : 
1. `sudo singularity build test.simg Singularity`
2. `sudo chown user:group test.simg`


## Intended Machine to Test
This is intended to test machine with : 
1. `x86_64` architecture, either AMD or Intel CPUs
2. Up to 20 cores are tested. More cores are acceptable, they just aren't tested in this code

These tests are not indended to test :
1. GPUs
2. FPGAs
3. The Filesystem.

## Running
Run at the command line:
`singularity exec --nv -H /home/group/user test.simg`

Please return to us :
1. `driver.log`
2. The output to `stdout`
