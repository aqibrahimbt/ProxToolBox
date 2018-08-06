Last Update by Alexander Dornheim: 11.02.2018

ProxPython is the Python version of the ProxToolbox.

How to get started?

The easiest way to use ProxPython is to run some of the included demos. There are different families of problems and various demos.
To run a demo open a command line. Change your current folder to ProxPython/TestSuite by

cd ProxPython/TestSuite

assuming ProxPython is in your current folder.

Then just type

python3 demo_name

to run the demo called "demo_name".

For example

python3 BLOCK_RAAR_parallel_ART_demo.py

will run a Computed Tomography demo using the RAAR algortihm and BLOCK wise projections.


Working demos (found in the TestSuite folder) are:

I. Phase Retrieval:
    1. JWST data set: James Webb Space Telescope simulation data set
        b. Phase_JWST_AltP_demo.py (alternating projections)
        c. Phase_JWST_RAAR_demo.py  (relaxed averaged alternating reflections proposed by Luke, 2005)
        
        Not implemented:
        a. Phase_JWST_RCAAR_demo.py (cyclic RAAR algorithm proposed by Borwein and Tam)       RCAAR not yet in ProxPython
        d. Phase_JWST_DRAP_demo.py  (Douglas-Rachford-Alternating Projectons hybrid proposed by Nguyen, 2017)     DRAP not yet in ProxPython

    2. Siemens data set:  far field Siemens star tests
		a.  Phase_Siemens_amplitude_demo.py (phase object)
        b.  Phase_Siemens_nonneg_demo.py (nonnegative object)
		c.  Phase_Siemens_real_demo.py (real object - possibly negative valued!?)

    3. Goettingen data set:
        a. tasse_ap_demo (far field optical diffraction - experimental data)
		b. tasse_raar_demo (far field optical diffraction - experimental data)
		c. Near_field_synthetic_demo (near field cell with noise - simulated)        
        d. Phase_near_field_Siemens_experimental_demo(near feild x-ray of Siemens star phase object with structured beam)
		e. Near_field_living_worm_demo (near field x-ray image of a living worm-  not sure parameters are correct) 

II. Computed Tomography:
    1.  BLOCK_RAAR_parallel_ART_demo.py (parallelized RAAR, Luke 2005, for CT)
    2.  BLOCK_RAAR_sequential_ART_demo.py (RAAR, Luke 2005, for CT)
    3.  ART_alternating_proj_demo.py (cyclic projections)
    4.  ART_averaged_proj_demo.py (Cimmino algorithm)
	5.  BLOCK_GRAAL_parallel_ART_demo.py (Golden Ratio Algorithm of Malitskyi, 2017)
    7.  BLOCK_AP_parallel_ART_demo.py (Block parallelized Cimmino)

    Not implemented:
    6.  BLOCK_ART_averaged_proj_in (Block Cimmino algorithm)    seems to do the same as demo 7 in ProxToolbox

V. Sudoku:
	1. Sudoku_demo.py (RAAR for Sudoku following Elser's 2007 model)



Structure

What is the basic structure of ProxPython?

I. proxtoolbox folder

This is the heart of the proxtoolbox. Most of the code is found here. There are four subfolders:

1. Algortihms
    This folder contains the Algortihm class and various instances of this class. The proxtoolbox allows the use of different algortihms. At moment the following algortihms are working:

    a) AP: alternating projections
    b) AP_expert: alternating projection
    c) GRAAL: Golden Ratio Algorithm of Malitskyi, 2017
    d) RAAR: Relaxed Averaged Alternating Reflection algorithm
    e) RAAR_expert: another instance of RAAR

There are other algorithms but these algortihms are either not working or have not been updated for some time. So it is not recommened to use the other algortihms at the moment. Some of them will be updated and some will be removed in future.
    
    Why are there two verison of AP, RAAR?

    There is a subclass SimpleAlgortihm of the Algortim class. AP and RAAR are instances of SimpleAlgorithm, while AP_expert and RAAR_expert are instances of the Algorithm class. The SimpleAlgortihm class handles the pre- and postprocessing and calculates and stores the diagnostics (change, gap, etc.). Therefore in AP and RAAR we only need to calculate the next iterate. This makes the algortim much easier to read and maintain, but some computations are done redundantly. Therefore the expert versions are more messy to read and maintain but should be faster.

2. Problems

    This folder contains the different problems to which the proxtoolbox can be applied. These problems are instances of the class Problem. At the moment two classes of problems are supported:
    a) CT
    b) Phase
    c) Sudoku
    You can find a list of the working demos at the top of this file.
    There are also ptychography.py and DRprotein.py. These files have not been updated for some time and are probably not working at the moment.
    
    What does the Problem folder contain:
    a) CT
        - ART.py contains the class for handeling ART problems.
        - *_in.py contain a dictionary which is used as input for our demos.
        - ART_graphics.py: graphic output
        - ART_data_processor.py: reads and prepares data
        - BLOCK_ART_data_processor.py: reads and prepares data, used with block wise ProxOperators
    b) Phase
        - phase.py  contains the class for handeling phase problems.
        - *_in.py contain a dictionary which is used as input for our demos.
        - *_graphics.py: graphic output
        - *_data_processor.py: reads and prepares data

3. ProxOperators
    Contains the ProxOperator class and various instances of ProxOperators.

4. Utilities
    Contains some Utilities, mostly for data processors.


II. InputData
When you first download the proxtoolbox this folder should be empty. When you call a ART/Phase demo the first time, you are asked if you want to automatically download the InputData for ART or Phase. If this does not work you need to download the data manually from:

http://vaopt.math.uni-goettingen.de/data/CT.tar.gz
http://vaopt.math.uni-goettingen.de/data/Phase.tar.gz

Extract the data and store it in
InputData/CT/
InputData/Phase/


III. TestSuite

Contains various demos, which you can run by

python3 demo_name.


This is the older README. Some of this information might be outdated.

So, you made it this far.  Congratulations. There are three paths to 
enlightenment:  the Way of Faith, the Way of Action and the Way 
of Contemplation.  This README follows the Way of Action to gain enlightenment
in python and proxtoolbox.  You've made the first step by just 
getting to the question.  

1. You need to get your confidence up by running a demo in python.  
So the first thing to do is open another X-terminal window and 
change directories to the same directory 
you are reading this from.  At the shell-prompt, type 

python3

This will open an instance of python3.##.  If your ipython interpreter
calls up version 3 of python, then you can type 

ipython

If you don't understand that last statement, skip it, it's not important.
Either way, if you see a different cursor prompt, then your in.  Congratulations, 
you have graduated to step 2!

2. Type 

import numpy as np

at the python prompt.   Don't ask why, just do it - it is not interesting.
But you should always do it without fail every time you open python for numerical 
mathematics.  It's like washing your hands before you eat, except the opposite, 
if that makes any sense, which it probably doesn't, so just do it and stop asking 
questions.  There will be a lot of these kind of command you will be told to type 
without explanation, so we'll just denote them JDI for ``just do it".  

3. Demo-specific JDI:  enter 

from proxtoolbox.Problems import Sudoku

at the python prompt.   

4. Demo-specific JDI: enter

ps=Sudoku()

at the python prompt.

5. Demo-specific JDI: type

ps.solve()

at the python prompt.

6. Demo-specific JDI: type
   
ps.show()

at the python prompt.

7. Did a graphics window appear with two Sudoku boards (one unsolved, one solved) and two graphs
underneath?  Congratulations, you've just taken your first step into a larger world!
 

Got that?  Great.  Let's try it again, with another demo.  First, to keep python running fast
lets clear out the current session by typing

exit()

at the python prompt.  Yup.  You're back in the X-window shell.  Repeat steps 1-2 above. 

3'.  Demo-specific JDI:  enter

from proxtoolbox.Problems import Ptychography

at the python prompt.

4'. Demo-specific JDI: enter

pp=Ptychography()

at the python prompt.

5'. Demo-specific JDI: enter

pp.solve()

at the python prompt.

6'. Demo-specific JDI: enter

pp.show()

at the python prompt.

7'. Did several graphics windows appear with lots of graphs and pictures? 
Congratulations, you've just taken your second step into a larger world! 
Ptyriffic isn't it? Ok, let's take a look at one final example... But first
don't forget to type

exit()

1''. For this example, we already done the hard work for you. This time we'll
run the ProxToolbox from using a "DRprotein_demo_demo.py" script. You can peek inside it
by openning it by your favourite editor, or even simpler type

cat DRprotein_demo_demo.py

See how its got all the algorithmic parameters already in the file. Pretty neat huh?

2''. Ok. Time to run the script! It going to take a bit longer... actually a lot longer
than the first two examples (hours maybe). We can work on speeding it up later. Type

python3 DRprotein_demo_demo.py

You should see some (boring) information about the problem instance, etc. While the algorithm
is running you a report of its progress even 10 iterations, and will stop when the relative
error drops below 10e-5. 

If you waited long enough "Jmol", a molecular viewer written in Java (which you should already
have installed) will open. Hopefully it looks like a molecule.

   
 
