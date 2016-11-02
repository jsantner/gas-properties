# gas-properties
GUI to simply and easily calculate gas properties using cantera. 

This program was created when I tried to help someone with a problem by having them install 
GasEq (Chris Morley, http://www.gaseq.co.uk/).
GasEq is a very useful GUI program for calculating gas-phase thermochemistry properties for those who are not comfortable with 
more advanced tools such as cantera and chemkin. Unfortunately, I found that GasEq is no longer supported with modern 
operating systems.

So, I set out to build a simple GUI in python that has similar functionality to GasEq


Installation:
This program is written for python 3. Install it via www.python.org.
For a more user-friendly installation experience that includes many dependencies, install Anaconda: https://www.continuum.io/downloads
This program also depends on cantera. See http://www.cantera.org/docs/sphinx/html/install.html for installation instructions


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

In order to change your "fixed" state variables (the red ones), delete all state variables except for the two that you 
want to fix.
