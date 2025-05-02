The two scripts will together extract and fit the cross-sections of agricultural drainage ditches. There are two seperate scripts that follow each other because in MATLAB it is more difficult to work with LAS files.

Step 1: Python file
This scripts needs three inputs: 

- An excell file that contains all the 2D coordinates of the endpoints of the different cross-sections you want to extract. The excel should look like this:

  id      X      Y
  1      20.42   4.37
  1      20.43   4.36
  2      18.63   5.91
  2      18.63   5.96
  .      .       .
  .      .       .
  .      .       .
  n      16.57  4.83
  n      16.53  4.81

  - A LAS file that preferably is already clipped to the areas just around the ditches in order to speed up the processing times (especially if working with point clouds with high point densities)
  - The buffer width which has a standard value of 0.35
 
There are a couple places in the code where the data is stored already to then use again for the next step. Do not forget to change the path names to your pc and give them logical names.
The final output needed in the MATLAB file is a text file with the 2D coordinates of the points classified to each cross-section you put in the excell file augmented with metadata if you would want to use these for something else
 
Step 2: MATLAB file

This script needs 1 input:

- The final text file generated in the Python script

Optional changes to parameters are the fitting parameters like knots, percentile, sections,...

The output of this script will be figures of all the cross-sections with the fitted point clouds and a table "Transecttable" which describes the course of each profile in 10000 steps.
