# miniAeroSeq

A sequential version of Mantevo's mini-app miniAero. See:
* https://github.com/Mantevo/miniAero
* https://mantevo.github.io/

## Issues

- (Possibly Stack Related) memory issue which causes an error in the `zero_cell_flux` functor - can be verified with valgrind/gdb.
- Code within `CELL_FLUX` preprocessor tags is unconverted - the makefile is configured to run without this flag.
- General memory leaks



