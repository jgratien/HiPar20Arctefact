# HiPar20Arctefact

This project is related to the experiments realized in the submitted paper :
"Introducing multi-level parallelism, at coarse, fine and instruction level to enhance the performance of
iterative solvers for large sparse linear systems on Multi and Many core~architecture"

There are two ditectories:
- data : contains several linear systems extracted from a reservoir simulation of the SPE10 comparative study case;
- src : contains tools to generate laplacian linear system from a unit cartesian grid of Nx*Ny*Nz blocks or to import systems from the files in data.

Extracted systems files :
- Spe10System-T63L1-UnNorm.xml and Spe10System-T63L1-UnNorm.h5
- Spe10x8System-T36L1-UnNorm.xml and Spe10x8System-T36L1-UnNorm.h5

HDF files (*.h5) can be download at \url{https://drive.google.com/drive/folders/1k-BcLST0ZoGuItKp3MeSRHPygFuNGy2g?usp=sharing}
