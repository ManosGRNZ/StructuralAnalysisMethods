# BiaxialResponseSpectrum

Python code that takes the acceleration ground motion records as an input and calculates the direct integration acceleration response spectrum for each direction (say Fault-Normal and Fault-Parallel or simply X and Y). At the same time it calculates the maximum acceleration for each T value (i.e. each single degree of freedom) and displays a "resultant" acceleration response spectrum.

It is obvious that the direction of the maximum, resultant acceleration does not have a single direction. Each value corresponds to a direction which the internal calculation can provide but is not included in the output. Mostly because I cannot think at the moment how to display this in a meaningful (visual) way, other than releasing a list of a hundred angles.  

It is important to note that the "resultant" acceleration response spectrum does not comply with any seismic code. This is a personal research project done "for the fun of it" rather than to be useful in any professional engineering setting. However, it does provide an insight about what the response of a structure may be under a given ground motion, especially when a 3D time history analysis is performed. 

This is a free (GNU/GPL v.3) project and everyone can use, enjoy and contribute under the terms of the license. 

Practical Notes: 

Mosts tests have been performed under Linux. 
It has been tested on Windows and it works.
Adjust the path variable for your computer. The provided file includes my paths for Linux and Windows. 
The default file runs many ground motions in parallel. You may wish to adjust the number of ground motions based on the parallel processing of your computer. Weak laptops still run the code as is but do not benefit much if all 11 ground motions run in parallel. It may be better to run 3 or 4 at a time. Strong CPUs will run fine. 
The default file runs each ground motion and displays a separate matplotlib window with the results. You can comment-out some sections and enable the last section to obtain the average acceleration spectrums of all ground motions, instead for each ground motion separately. This method does not use parallel processing and may take a while to run. 
