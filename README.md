Add these to "modules/":

- genetic_algorithm.h
- genetic_algorithm.cpp
- params.cpp

You can run from command line or an IDE but I found that the easiest way to run was from bash, so I updated the Makefile. Simply replace the Makefile with the one uploaded into this repository and type into bash:

`make gen`

`./gen.exe`

The Makefile has been optimized for running on macOS ARM64, and I had to work around the way the nlohmann::json library was set up originally. 
