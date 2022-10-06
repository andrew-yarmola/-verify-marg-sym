# Symmetric Margulis constant verification

This repository contains the verification code for finding the symmetric Margulis constant
of orientable hyperbolic 3-manifolds. 

Data will be uploaded soon.

# Installation

There are no dependencies of this code outside of the C++ standard library.
There are several binaries that need to be made from the `src` directory:

```
make rootcat
make verify
make tests
```
One should run `test_float` from the `bin` directory before any verification code to make sure that your system correctly catchest any overflow and underflow.
To run the verification, from the `bin` direcotry execute:
```
./rootcat ../data/verify | ./verify
```
Note, this process is long and take 10-15 hours to finish.

For ease of use, script are provided in the scripts directory that perform the above calls alongside some pretty printing.

# Code

## Arithmetic

The arithmetic used in this repository is a modified version of the one published in "Homotopy hyperbolic 3-manifolds are hyperbolic."
In particular, we change the dimensionality and add complex conjugation.
This produces an arithemtic with 8 real parameters.
Note, in "Homotopy hyperbolic 3-manifolds are hyperbolic," the authors use `CWEB` to generate the source code.
For the new code in this repository, we chose not to use `CWEB` as this adds an additonal step and dependency.
Discussion of the details of the code are in the appendix of the main paper.

## Parameter space code

We use the arithmetic above to verify a binary tree, the leaf nodes of which represent (sub)boxes in the parameter space and corresponding eliminaton criteria. See our paper for deals.

### rootcat

This program reads the binary tree as stored in a collection of text `.out` files, or compressed `.out.tar.gz` versions, and prints the tree out to `stdin` in depth-first format.
This data format was chosen because the dataset is quite large.
Each internal node is represented by an `X` and each leaf node includes an elinimation criterion that will be checked by `verify`.

The code for `rootcat` is entirely contained in `rootcat.c`.
It is written in `C` and should catch any and all reading errors if the data is corrupt or the tree is incomplete.

### verify

This program verifies Proposition TBD.
It expects to recieve the complete data tree in depth-first order (e.g. the output of `rootcat`) via `stdin`.
At a leaf node, the program will take the binary coodrinate of the node and use it to construct the parameters corresponding to that box.
See the relevant code in `box.h` and `box.c`.
Validity of this code is discussed in the paper.

The program then rigorously validates the the condition encoded in the leaf node elimination criterion holds over the entire box.
Conditions 0 and 1 are boundary conditions that check if the box lies entirely outside of the valid parameter space on interest.
Condition 2 is an emedded tube condition derived from the work of Rob Meyerhoff.
Lettered conditions are of the form `t(words)` where `t` encodes the type of condition and `words` is a word or list of words.
There are seven types of lettered contitions, see the paper for details.
The code for verifying these conditions can be found in: `elimination.h` and `elimination.c`.

### tests

These programs checks whether your system correctly report that roundoff error has occured.
Note, the roundoff checking function has been updated from the version in "Homotopy hyperbolic 3-manifolds are hyperbolic" to only include testing for x64 machines. 
