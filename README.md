# AST9240 -- Cosmology II

The repo contains code(s) created based on the Fortran templates provided for AST9240 -- Cosmology II course from University of Oslo. It consists of 4 projects in total:

1. _The background evolution of the Universe_
2. _The recombination history of the Universe_
3. _The evolution of structures in the Universe_
4. _The CMB power spectrum_

Reports for all 4 parts located in `reports` folder within this repo.

The assignments and some materials can be found inside `materials` folder within this repo.

The total list of materials for the course can be found on official course page here:
https://www.uio.no/studier/emner/matnat/astro/AST5220/v19/pensumliste/index.html

## Prerequisites

To install the code you need to have the following libraries preinstalled first:

1. CFitsIO;
2. HEALPix;
3. BLAS/LAPACK;

You will also need to have a Fortran compiler. The primary compiler used was from Intel, but GNU GFortran should also work, given the correct compiler flags.

Once you installed all prerequisites, modify `Makefile` to incorporate the changes.

## Compile

To compile the project simply type:
```
$ make
```
in the root. This will produce executable called `cmbspec`.

Note: I have modified the original template `Makefile` to split `src` codes from compiled `mod` and `obj` files. This means that you will need to create new directory inside root called `obj` for `make` to work. 

## Run

To run the code type:
```
$ ./cmbspec
```
It will start all necessary processing. The resulting files will be stored inside `data` folder within root.

To plot the results for each subproject, use `ipynb` contained within this repo.

Note: you need to create `data` directory before you started processing otherwise it will not work.
