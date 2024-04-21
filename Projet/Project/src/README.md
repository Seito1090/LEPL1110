# Project LEPL1110, Groupe 55, déformation des ailes d'un avion

## Installation
To make the project work on your machine follow these steps:
1. Make sure you are in the directory where CMakeLists.txt is located / at the root of the project.
2. Simply run this command : 'cmake -B build/ ' to make the build and all the necessary files to build the project.
3. Go to the build directory and run 'make' to build the project.
4. You can now run the project with the 'myFem' executable.

OR open the root of this project in vscode, and click on the "Run" button on the bottom left corner of the screen.

Warning: To be sure that the project finds your data files, make sure to run the executable from the root of the project, in this case it's just ./build/myFem.


## Data
The code will read the data from the 'data' directory located in the src directory of the project as this was the structure asked. Just make sure to put there 2 files with the following names:
- 'mesh.txt' : the mesh file 
- 'problem.txt' : the data file

In any case, the preProcessor should create these files for you.
