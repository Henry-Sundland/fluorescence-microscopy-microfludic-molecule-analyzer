# fluorescence microscopy microfludic molecule analyzer
 For a molecule on frame of a tif video, tracks the molecule via its center of mass, calculates mean squared displacement, allows user to fit to find diffusion, calculates molecular weight based off of inputted bp/intensity calibration, calculates gyration tensor and calculates radius of gyration. Saves all of this info to a cell array

 Code loops through tif file videos (of one molcule) in the folder that these three codes are placed into. You'll have to paste the path to the folder containing the tif videos. You'll have to enter a calibration (i.e. an average intensity set to a known base pair to create a conversion of bp/intensity) in order to convert to molecular weight based off of intensity estimate. You will also have to input the pixel length in units of nanometers. The end result will be in fundamental units of micrometers, for calculations.


 All three codes need to be in the same folder. The main code is the "microfluidic_flourescence....bla bla" code. Just open that, input the values asked and run that script.


 Any questions, then please message me on here or on my email henryss634@gmail.com . This code is what I wrote to analyze my research data in grad school, and I'm open to revising it, making it more general and easy to use or explaining it in more detail. I want to also make a python version of this code to put on here.

 Code goes like this:
    - When all three are placed into folder, and you run the main code, it asks you to: enter intensity molecular weight conversion and pixel length in units of nm 
    - Then it asks you if you would like to watch the videos of your molecules moving about on each frame
    - code eliminates noise by use of a gaussian blur method and noise subtraction method
    - calculates center of mass based off of intensity of molecule
    - calculates gyration tensor based off of center of mass location
    - calulates radius of gyration of molecule from eigen values of the 2d gyration tensor (units of micro meters)
    - calculates mean squared displacment via particle tracking the center of mass on each frame of tif video
    - displays mean squared displacement in x dimension and the y dimension, then asks user which curve they would like to fit to to find diffusion value
    - the asks which range user would like to fit a line to (the most linear section...make sure to have a good amount of linear points)
    - calculates diffusion value of this fit, as well as the uncertainty
    - saves all of this stuff into a cell array. Look at the comments at the end of the code to see which cell array index refers to which data
    - saves cell array to folder that the code is placed in
