# gas-properties
GUI to simply and easily calculate gas properties using cantera.

![image](https://user-images.githubusercontent.com/16088743/121784497-760c4700-cbb4-11eb-9342-9f23ec316f84.png)

This program was created when I tried to help someone with a problem by having them install
GasEq (Chris Morley, http://www.gaseq.co.uk/).
GasEq is a very useful GUI program for calculating gas-phase thermochemistry properties for those who are not comfortable with
more advanced tools such as cantera and chemkin. Unfortunately, I found that GasEq is no longer supported with modern
operating systems.

So, I set out to build a simple GUI in python that has similar functionality to GasEq

I recently learned of another GUI, ASALI (https://github.com/srebughini/ASALI). ASALI provides some functionality
that is missing from gas-properties, and vice versa. They also have different interface designs (of course).


Installation:
-------------
This program is written for python 3. Install it via www.python.org.
For a more user-friendly installation experience that includes many dependencies, install Anaconda: https://www.anaconda.com/distribution/
This program also depends on cantera. See http://www.cantera.org/docs/sphinx/html/install.html for installation instructions.
I think the easiest installation method is to open an anaconda prompt and type:

    conda install -c cantera cantera


To run:
-------
From a windows command line:

    python -i PATH\TO\gas_properties.py

Or, run gas-properties.py from your preffered python IDE (e.g. spyder, IDLE, etc.)

Usage:
In the top left, you can input molecules (H2, O2, N2, CH4, iC8H18, nC7H16, etc.) and their relative molar content.
Choose two state variables to define your initial state and click the calculate button. The variables you choose are in red,
calculated values are in black. This will also calculate equivalence ratio, which you can edit in order to change the mixture.

You can then calculate the final condition in three ways
  1. Chemically frozen: choose two state variables for final conditions, and the other properties will be calculated assuming
      no reaction has occured.
  2. Isentropic compression: input a compression ratio, and the final entropy and density will be fixed based
      on this compression ratio. Other properties are then calculated.
  3. Chemical Equilibrium: Choose two state variables for the final condition, and the other properties will be calculated
      after allowing the mixture to reach chemical equilibrium.

The "fixed" state variables (the red ones), are determined using the radio buttons (TP, UV, etc.).
If you want these state variables to be determined, but not held constant from the initial state, use "manual," which allows
the user to set two properties for the final state.
