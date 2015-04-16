converterIO
===========

File Conversion Tool for PHASTA Files. 

1) Description:

This library provides two IO converters for PHASTA files, from the original Posix format to the SyncIO format and vice versa.
The original Posix format is characterized by 1 geombc and 1 restart file per mesh part, whereas SyncIO relies on MPI-IO and can handle several parts per file.

2) Compilation

Both IO converters are compiled by default along with the phasta executable using cmake.
Note that the phastaIO library is also required to compile this library
Both conversion executables are named respectively converterO2N (Posix->SyncIO) and converterN2O (SyncIO->Posix).

3) Usage

a) converterO2N

Let us assume we have already generated a typical #-procs_case directory that contains phasta files under the Posix format (one single geombc and restart file per part), where # is the number of parts in the partitioned mesh.
Let us also assume that we want @ SyncIO files, where "@" is an integer.
The first requirement is that mod(#,@)=0 so that we end up #/@ parts per SyncIO file.

In order to converter the #-procs_case directory, we simply need to run the converterO2N executable in the following way:

  mpirun -np NP converterO2N

which will read the PHASTA files in the #-procs_case directory and generate a new directory named #-procs_case-SyncIO-@ which contains @ new geombc and @ new restart files under the SyncIO format.
The number of processes NP can be relatively flexible but must follow the following rules:
i) @ <= NP <= #
ii) mod(NP,@)=0
iii) nod(#,NP)=0

For instance, if we have 16 parts and we want 2 SyncIO files (8 parts per SyncIO files), then NP can be equal to 2, 4, 8 and 16.
Note that PHASTA always expects a directory named #-procs_case so that the output directory #-procs_case-SyncIO-@ has to be renamed accordingly.

ConverterO2N also requires an input file named IO.O2N.input. This input files includes
- the number of part in the mesh (i.e. #), 
- the time step number of the restart files, 
- the number of output SyncIO files (i.e. @),
- the number of fields (integer and double) in both the geombc and restart files to convert, and
- the list of all the fields in both the geombc and restart files to convert.

It can be a tidious task to list all the fields that are essential for the PHASTA simulation. 
No worry, there is a bash script which can generate automatically this input file, which is provided with this converterIO code.
The script's name is create_IO_O2N_input.sh. Run this script a first time without any argument to learn about its usage. 
It basically requires the time step number of the restart files, the number of parts #, and the number of output SyncIO files.
We can always add more fields if needed, especially for the restart files (ybar, dwal, error, etc). 
Note that the double fields for either the geombc or restart files must always preceed the integer field in the list and the fields counter needs also to be increased in the header of IO.O2N.input.

Beware! This script relies on grep, awk, etc and needs to be applied to ALL the phasta files contained in the #-procs_case directory.
Indeed, for mixed topology meshes, there is no garantee that a single geombc file will contain all the topologies present in the mesh.
Therefore, it is strongly advised to run this script on a limited part count mesh (up to a few thousands is ok) and not on multiple million part meshes (don't do that!).
Instead, we can always recycle the same input file from a small partitioning to a large one as long as the topology of the mesh has not changed.

b) converterN2O

Let us assume again that we have a #-procs_case directory that contains @ phasta files under the SyncIO format this time.

In order to deconverter the #-procs_case directory from the SyncIO format to the Posix format, we simply need to run the converterN2O executable in the following way:

  mpirun -np # converterN2O

which will read the PHASTA files in the #-procs_case directory and generate a new directory named #-procs_case-1PPP which contains # geombc and # restart files under the Posix format.

Although there is a serial reader for SyncIO file that is used in Paraview for instance, converterN2O can currently read the phasta files under the SyncIO format only collectively.
This means that the number of processes in this case must be equal to the total number of parts, like for PHASTA.
This could be improved in the future for more flexibility like for converterO2N.

Similarly to create_IO_O2N_input.sh, there is also a bash script named create_IO_N2O_input.sh which can generate the input file for converterN2O.
This script is simpler since only the fields in the restart file really need to be deconverted (blocks and double).
Indeed, you should already have the original geombc files under the Posix format somewhere so it is useless (and time consuming) to generate them again.
The usage of this script is the same as above. 
By default, it will add to the input file all the blocks present in the restart files and assume they are all double.
Afterwards, we can manually adapt this script for any specific need.


Questions, comments or suggestions are always welcome.

Feel free to contact:

  Michel Rasquin
  michel.rasquin@colorado.edu























