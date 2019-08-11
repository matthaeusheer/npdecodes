# NPDECODES

This a git repository for the flipped-classroom course "Numerical Methods for Partial Differential Equations" taught at ETH Zurich by Prof. Ralf Hiptmair, Seminar for Applied Mathematics.

This repository is intended for the C++ codes associated with homework problems and
examples shown in class.

## Subdirectories

- `homeworks`: contains many different directories, all with "telltale" names, each
  containing a particular homework problem.
- `lecturecodes`: contains non-LehrFEM++-based C++ codes for (partial) inclusion in the
  lecture material


## Installation of MathGL2

If not present on your machine, MathGL2 will be installed through the provided `CMakeLists.txt`.
Some dependencies, however, are not automatically installed. You will have to install them yourself, if 
you don't already have them, see the following paragraph.

### `Couldn't find ZLIB/PNG/OPENGL`

A frequent error message produced during the installation of MathGL is of the form
```
-- Could NOT find PNG (missing: PNG_LIBRARY PNG_PNG_INCLUDE_DIR) 
CMake Error at CMakeLists.txt:341 (message):
  Couldn't find PNG library.
```
or `ZLIB` or `OPENGL` instead of `PNG`.
You'll need the development versions of these packages.
After you have them installed, remove your previous build and run `cmake` anew.

For Fedora (tested):
* `OPENGL`: `sudo dnf install freeglut-devel`
* `PNG`: `sudo dnf install libpng-devel`
* `ZLIB`: `sudo dnf install zlib-devel`

For Ubuntu (not tested):
* `OPENGL`: `sudo apt install freeglut3-dev`
* `PNG`: `sudo apt install libpng-dev` (or `libpng++-dev`)
* `ZLIB`: `sudo apt install zlib1g-dev`

For Mac (tested):
Make sure Xcode 10.1 or higher is installed and run 
`xcode-select --install`.
* `PNG`: `brew install libpng`
* (`ZLIB`: `brew install zlib-devel`, probably not necessary)