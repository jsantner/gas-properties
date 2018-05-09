""" GUI to calculate gas properties based on cantera.

Designed with similar functionality to GasEq

Created by Jeffrey Santner
"""

import os
import sys
import tkinter
from tkinter import filedialog
from tkinter import ttk
from tkinter import messagebox
import cantera


class GasProps(ttk.Frame):
    "Main window"
    fmt = '{:.6g}'  # Number format used throughout
    short_fmt = '{:.4g}'  # Number format used throughout

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.root = parent
        self.init_gui()

    def init_gui(self):
        "Build the GUI."
        self.root.title('Gas Property Calculator')
        self.grid(column=0, row=0, sticky='nsew')

        # Initialize Variables
        self.initial_type = None
        self.final_type = None

        # Mixture creation
        species = []
        species_vals = []
        num_species = 9  # Number of species to include
        for i in range(num_species):
            # Create species buttons
            specie = tkinter.StringVar()
            spec = ttk.Entry(self, width=10, textvariable=specie)
            spec.grid(column=0, row=i+1)
            val = tkinter.DoubleVar()
            vals = ttk.Entry(self, width=12, textvariable=val)
            vals.grid(column=1, row=i+1)
            species.append(specie)
            species_vals.append(val)
        self.species = species
        self.species_vals = species_vals
        int_row = i + 1  # Intermediate row
        self.phi = tkinter.StringVar()
        self.phi_entry = ttk.Entry(self, width=10, textvariable=self.phi)
        self.phi_entry.grid(column=1, row=int_row+1)
        self.phi.trace('w', self.change_phi)  # Run when value changes
        ttk.Label(self, text='Equivalence Ratio').grid(column=0, row=int_row+1)
        ttk.Label(self, text='Species').grid(column=0, row=0)
        ttk.Label(self, text='Relative Moles').grid(column=1, row=0)

        # Separator
        ttk.Separator(self, orient='horizontal').grid(column=0, row=int_row+2,
                                                      columnspan=7,
                                                      sticky='ew')
        int_row += 3  # Move reference row down a bit
        # Initial / Calculated Conditions
        self.T = tkinter.StringVar()
        self.P = tkinter.StringVar()
        self.H = tkinter.StringVar()
        self.S = tkinter.StringVar()
        self.D = tkinter.StringVar()
        self.U = tkinter.StringVar()

        self.T_entry = ttk.Entry(self, width=12, textvariable=self.T)
        self.P_entry = ttk.Entry(self, width=12, textvariable=self.P)
        self.H_entry = ttk.Entry(self, width=12, textvariable=self.H)
        self.S_entry = ttk.Entry(self, width=12, textvariable=self.S)
        self.D_entry = ttk.Entry(self, width=12, textvariable=self.D)
        self.U_entry = ttk.Entry(self, width=12, textvariable=self.U)
        self.initial_entries = [self.T_entry, self.P_entry, self.H_entry,
                                self.S_entry, self.D_entry, self.U_entry]

        self.T_entry.grid(column=0, row=int_row+1)
        self.P_entry.grid(column=0, row=int_row+2)
        self.H_entry.grid(column=0, row=int_row+3)
        self.S_entry.grid(column=0, row=int_row+4)
        self.D_entry.grid(column=0, row=int_row+5)
        self.U_entry.grid(column=0, row=int_row+6)

        ttk.Label(self, text='Initial Conditions').grid(column=0, row=int_row)

        self.initial_button = ttk.Button(
            self, text='Calc. Initial Conditions',
            command=self.calculate_initial)
        self.initial_button.grid(column=0, row=int_row+7)

        # Final Conditions
        ttk.Separator(self, orient='vertical').grid(column=3, row=0,
                                                    rowspan=25, sticky='ns')

        self.T_f = tkinter.StringVar()
        self.P_f = tkinter.StringVar()
        self.H_f = tkinter.StringVar()
        self.S_f = tkinter.StringVar()
        self.D_f = tkinter.StringVar()
        self.U_f = tkinter.StringVar()

        self.T_f_entry = ttk.Entry(self, width=12, textvariable=self.T_f)
        self.P_f_entry = ttk.Entry(self, width=12, textvariable=self.P_f)
        self.H_f_entry = ttk.Entry(self, width=12, textvariable=self.H_f)
        self.S_f_entry = ttk.Entry(self, width=12, textvariable=self.S_f)
        self.D_f_entry = ttk.Entry(self, width=12, textvariable=self.D_f)
        self.U_f_entry = ttk.Entry(self, width=12, textvariable=self.U_f)
        self.final_entries = [self.T_f_entry, self.P_f_entry, self.H_f_entry,
                              self.S_f_entry, self.D_f_entry, self.U_f_entry]

        self.T_f_entry.grid(column=2, row=int_row+1)
        self.P_f_entry.grid(column=2, row=int_row+2)
        self.H_f_entry.grid(column=2, row=int_row+3)
        self.S_f_entry.grid(column=2, row=int_row+4)
        self.D_f_entry.grid(column=2, row=int_row+5)
        self.U_f_entry.grid(column=2, row=int_row+6)

        ttk.Label(self, text='Final Conditions').grid(column=2, row=int_row)
        ttk.Label(self, text='Temperature (K)').grid(column=1, row=int_row+1)
        ttk.Label(self, text='Pressure (Pa)').grid(column=1, row=int_row+2)
        ttk.Label(self, text='Enthalpy (J/kg)').grid(column=1, row=int_row+3)
        ttk.Label(self, text='Entropy (J/kg/K)').grid(column=1, row=int_row+4)
        ttk.Label(self, text='Density (kg/m^3)').grid(column=1, row=int_row+5)
        ttk.Label(self, text='Energy (J/kg)').grid(column=1, row=int_row+6)

        self.final_button = ttk.Button(
            self, text='Calc. Final Conditions',
            command=self.calculate_final)
        self.final_button.grid(column=2, row=int_row+7)

        # Menu bar
        self.root.option_add('*tearOff', 'FALSE')
        self.menubar = tkinter.Menu(self.root)
        self.menu_file = tkinter.Menu(self.menubar)
        self.menu_file.add_command(label='Exit', command=self.on_quit)
        self.menu_file.add_command(label='Open Chem. File', command=self.openf)
        self.menu_edit = tkinter.Menu(self.menubar)
        self.menubar.add_cascade(menu=self.menu_file, label='File')
        self.root.config(menu=self.menubar)

        # Problem type radio buttons
        ttk.Label(self, text='Problem Type').grid(column=4, row=0)
        self.problem_type = tkinter.IntVar()
        ttk.Radiobutton(self, text='Chemically Frozen',
                        variable=self.problem_type, value=0,
                        command=self.change_ptype).grid(column=4, row=1)
        ttk.Radiobutton(self, text='Isentropic Compression. cr:',
                        variable=self.problem_type, value=1,
                        command=self.change_ptype).grid(column=4, row=2)
        ttk.Radiobutton(self, text='Chemical Equilibrium',
                        variable=self.problem_type, value=2,
                        command=self.change_ptype).grid(column=4, row=3)
        self.cr = tkinter.StringVar()
        self.cr_Entry = ttk.Entry(self, width=5, textvariable=self.cr)
        self.cr_Entry.grid(column=5, row=2)
        self.cr.trace('w', self.change_cr)

        self.final_type_int = tkinter.IntVar()
        self.finaltype_buttons = [
                ttk.Radiobutton(self, text='TP', variable=self.final_type_int,
                                command=self.change_ftype, value=0),
                ttk.Radiobutton(self, text='HP', variable=self.final_type_int,
                                command=self.change_ftype, value=1),
                ttk.Radiobutton(self, text='SP', variable=self.final_type_int,
                                command=self.change_ftype, value=2),
                ttk.Radiobutton(self, text='SD', variable=self.final_type_int,
                                command=self.change_ftype, value=3),
                ttk.Radiobutton(self, text='UV', variable=self.final_type_int,
                                command=self.change_ftype, value=4),
                ttk.Radiobutton(self, text='Manual',
                                variable=self.final_type_int,
                                command=self.change_ftype, value=5)
                ]

        count = 5
        for button in self.finaltype_buttons:
            button.grid(column=4, row=count)
#            button.state(statespec=['disabled'])
            count += 1

        # Calculated mixture properties
        self.calc_labels = ['Cp (J/kg/K)', 'Cv (J/kg/K', 'gamma',
                            'Gibbs (J/kg)', 'MW (kg/kmol)',
                            'IsoT compressibility (1/Pa)',
                            'Therm. expansion coeff (1/K)']

        self.calc_calls = ['cp', 'cv', 'GAMMA', 'g', 'mean_molecular_weight',
                           'isothermal_compressibility',
                           'thermal_expansion_coeff']
        i = int_row+8
        initials = []
        finals = []
        for label in self.calc_labels:
            init = tkinter.StringVar()
            final = tkinter.StringVar()
            ttk.Entry(self, width=12, textvariable=init,
                      state='readonly').grid(column=0, row=i)
            ttk.Entry(self, width=12, textvariable=final,
                      state='readonly').grid(column=2, row=i)
            ttk.Label(self, text=label).grid(column=1, row=i)
            i += 1
            initials.append(init)
            finals.append(final)
        self.initial_calc = initials
        self.final_calc = finals

        # Final Composition
        ttk.Label(self, text='Final Composition').grid(column=4, row=int_row,
                                                       columnspan=2)
        ttk.Label(self, text='Species').grid(column=4, row=int_row+1)
        ttk.Label(self, text='Mole Fraction').grid(column=5, row=int_row+1)
        final_spec = []
        final_mole_frac = []
        for ind in range(int_row+2, int_row+15):
            spec = tkinter.StringVar()
            val = tkinter.StringVar()
            ttk.Entry(self, width=10, textvariable=spec,
                      state='readonly').grid(column=4, row=ind)
            ttk.Entry(self, width=12, textvariable=val,
                      state='readonly').grid(column=5, row=ind)
            final_spec.append(spec)
            final_mole_frac.append(val)
        self.final_spec = final_spec
        self.final_mole_frac = final_mole_frac

        # Initialize chemistry file to therm.cti
        self.chemfile = os.path.join(sys.path[0], 'therm.cti')

    def openf(self):
        """ Opens chemistry file """
        # Allows user to choose thermochemistry file from Menu.
        self.chemfile = tkinter.filedialog.askopenfilename()

    def on_quit(self):
        """ Exits program. """
        quit()

    def change_ptype(self):
        """ Change the problem type. """
        problem_type = self.problem_type.get()

        if problem_type == 0:
            # Chemically frozen
            spec = 'normal'
        elif problem_type == 1:
            # Isentropic Compression
            spec = 'disabled'
            self.final_type_int.set(5)
        elif problem_type == 2:
            # Chemical Equilibrium
            spec = 'normal'

        for button in self.finaltype_buttons:
            button['state'] = spec

        self.change_ftype()  # Update

    def change_ftype(self):
        """ Change the final type (TP, HP, SP, SD, UV) """
        # Use radio buttons to write known final state variables
        final_type_int = self.final_type_int.get()

        # Reset all text to black
        [x.config(foreground='black') for x in self.final_entries]

        # Delete all final states
        for item in self.final_entries:
            item.delete(0, tkinter.END)

        if final_type_int == 0:
            self.T_f.set(self.T.get())
            self.P_f.set(self.P.get())
        if final_type_int == 1:
            self.H_f.set(self.H.get())
            self.P_f.set(self.P.get())
        if final_type_int == 2:
            self.S_f.set(self.S.get())
            self.P_f.set(self.P.get())
        if final_type_int == 3:
            self.S_f.set(self.S.get())
            self.D_f.set(self.D.get())
        if final_type_int == 4:
            self.U_f.set(self.U.get())
            self.D_f.set(self.D.get())

        # Reset calculated properties
        for var, attr in zip(self.final_calc, self.calc_calls):
            var.set('')
        for i in range(len(self.final_spec)):
            self.final_spec[i].set('')
            self.final_mole_frac[i].set('')
            

    def change_cr(self, *args):
        """ Update items when compression ratio is changed. """
        if self.cr.get()[-1] == 0:
            # Don't update if last character is zero.
            return
        try:
            self.S_f.set(self.S.get())
            self.D_f.set(float(self.D.get()) * float(self.cr.get()))
        except ValueError:
            self.S_f.set('')
            self.D_f.set('')

    def calculate_initial(self):
        """ Calculate initial conditions and display results """

        # Normalize and create mixture
        gas_list = []
        for spec, val in zip(self.species, self.species_vals):
            if spec.get() is not '' and val.get() is not '':
                gas_list.append([spec.get(), float(val.get())])
        total = sum([x[1] for x in gas_list])
        self.gas_list = [[x[0], x[1]/total] for x in gas_list]
        reactants = ', '.join(['{:}:{:.5f}'.format(*x) for x in gas_list])

        # Read thermodynamic state variables
        inputs = [self.T, self.P, self.H, self.S, self.D, self.U]
        outputs = []
        for inp in inputs:
            out = inp.get()
            try:
                out = float(out)
            except ValueError:
                out = None
            if out == 0:
                out = None
            outputs.append(out)
        T, P, H, S, D, U = outputs

        # Find the problem type
        if len([x for x in outputs if x is not None]) == 2:
            if T is not None and P is not None:
                self.initial_type = 'TP'
            elif H is not None and P is not None:
                self.initial_type = 'HP'
            elif S is not None and P is not None:
                self.initial_type = 'SP'
            elif S is not None and D is not None:
                self.initial_type = 'SD'
            elif T is not None and D is not None:
                self.initial_type = 'TD'
            elif U is not None and D is not None:
                self.initial_type = 'UV'
            else:
                tkinter.messagebox.showerror('Error', 'Must input exactly two '
                                             'state variables '
                                             '(TP, HP, SP, SD, UV)')
                return
        elif self.initial_type is None:
            tkinter.messagebox.showerror('Error', 'Must input exactly two '
                                         'state variables '
                                         '(TP, HP, SP, SD, UD)')
            return

        # Create the cantera Solution object
        gas = cantera.Solution(self.chemfile)

        # Reset all text to black
        [x.config(foreground='black') for x in self.initial_entries]

        try:
            if self.initial_type == 'TP':
                gas.TPX = T, P, reactants
                self.T_entry.config(foreground='red')
                self.P_entry.config(foreground='red')
            elif self.initial_type == 'HP':
                gas.HPX = H, P, reactants
                self.H_entry.config(foreground='red')
                self.P_entry.config(foreground='red')
            elif self.initial_type == 'SP':
                gas.SPX = S, P, reactants
                self.S_entry.config(foreground='red')
                self.P_entry.config(foreground='red')
            elif self.initial_type == 'SD':
                gas.SVX = S, 1/D, reactants
                self.S_entry.config(foreground='red')
                self.V_entry.config(foreground='red')
            elif self.initial_type == 'TD':
                gas.TDX = T, D, reactants
                self.T_entry.config(foreground='red')
                self.D_entry.config(foreground='red')
            elif self.initial_type == 'UV':
                gas.UVX = U, 1/D, reactants
                self.U_entry.config(foreground='red')
                self.D_entry.config(foreground='red')
        except RuntimeError as e:  # Cantera has a problem
            tkinter.messagebox.showerror('Error', e)
            [x.set('') for x in inputs]  # Reset all state variables
        else:
            self.gas = gas

            # Display Results
            self.phi.set(self.short_fmt.format(
                    equivalence_ratio(gas, gas_list)[0]))
            H, P = gas.HP
            T, D = gas.TD
            S = gas.s
            U = gas.u
            [x.set(self.fmt.format(y))
             for x, y in zip(inputs, [T, P, H, S, D, U])]

            # More calculated properties
            for var, attr in zip(self.initial_calc, self.calc_calls):
                if attr == 'GAMMA':
                    var.set(self.fmt.format(gas.cp/gas.cv))
                else:
                    var.set(self.fmt.format(getattr(gas, attr)))

    def calculate_final(self):
        """ Calculate the final conditions. """

        # Update Initial Conditions
        self.calculate_initial()

        problem_type = self.problem_type.get()

        # Read thermodynamic state variables
        inputs = [self.T_f, self.P_f, self.H_f, self.S_f, self.D_f, self.U_f]
        outputs = []
        for inp in inputs:
            out = inp.get()
            try:
                out = float(out)
            except ValueError:
                out = None
            if out == 0:
                out = None
            outputs.append(out)
        T, P, H, S, D, U = outputs

        # Find the problem type
        if problem_type == 1:  # Isentropic Compression
            T, P, H, S, D, U = [None]*6
            S = float(self.S.get())  # Final entropy = initial entropy
            try:
                # Change density due to compression / expansion
                D = float(self.D.get()) * float(self.cr.get())
            except ValueError:
                tkinter.messagebox.showerror('Error', 'compression ratio'
                                             'must be a number')
                return
        if len([x for x in [T, P, H, S, D, U] if x is not None]) == 2:
            # This block could be simplified and just look at final_type_int.
            # But, doing it this way serves as a check
            if T is not None and P is not None:
                self.final_type = 'TP'
            elif H is not None and P is not None:
                self.final_type = 'HP'
            elif S is not None and P is not None:
                self.final_type = 'SP'
            elif S is not None and D is not None:
                self.final_type = 'SD'
            elif T is not None and D is not None:
                self.final_type = 'TD'
            elif U is not None and D is not None:
                self.final_type = 'UV'
            else:
                tkinter.messagebox.showerror('Error', 'Must input exactly two '
                                             'state variables '
                                             '(TP, HP, SP, SD, UV)')
                return
        elif len([x for x in [T, P, H, S, D, U] if x is not None]) == 0:
            self.change_ptype()  # Update problem type.
        elif self.final_type is None:
            tkinter.messagebox.showerror('Error', 'Must input exactly two '
                                         'state variables '
                                         '(TP, HP, SP, SD, UV)')
            return

        # Reset all text to black
        [x.config(foreground='black') for x in self.final_entries]

        # Change state variables, re-calculate gas state
        try:
            if self.final_type == 'TP':
                self.gas.TP = T, P
                self.T_f_entry.config(foreground='red')
                self.P_f_entry.config(foreground='red')
            elif self.final_type == 'HP':
                self.gas.HP = H, P
                self.H_f_entry.config(foreground='red')
                self.P_f_entry.config(foreground='red')
            elif self.final_type == 'SP':
                self.gas.SP = S, P
                self.S_f_entry.config(foreground='red')
                self.P_f_entry.config(foreground='red')
            elif self.final_type == 'SD':
                self.gas.SV = S, 1/D
                self.S_f_entry.config(foreground='red')
                self.D_f_entry.config(foreground='red')
            elif self.final_type == 'TD':
                self.gas.TD = T, D
                self.T_f_entry.config(foreground='red')
                self.D_f_entry.config(foreground='red')
            elif self.final_type == 'UV':
                self.gas.UV = U, 1/D
                self.U_f_entry.config(foreground='red')
                self.D_f_entry.config(foreground='red')
            if problem_type == 2:  # Equilibrate
                equil_key = self.final_type.replace('D', 'V')
                self.gas.equilibrate(equil_key)

        except RuntimeError as e:  # Cantera has a problem
            tkinter.messagebox.showerror('Error', e)
            [x.set('') for x in inputs]  # Reset all state variables
        else:
            gas = self.gas
            # Display Results
            H, P = gas.HP
            T, D = gas.TD
            S = gas.s
            U = gas.u
            [x.set(self.fmt.format(y)) for x, y in
             zip(inputs, [T, P, H, S, D, U])]

            # More calculated properties
            for var, attr in zip(self.final_calc, self.calc_calls):
                if attr == 'GAMMA':
                    var.set(self.fmt.format(gas.cp/gas.cv))
                else:
                    var.set(self.fmt.format(getattr(gas, attr)))

            # Final mole fractions
            i = 0
            d = gas.mole_fraction_dict()
            for k in sorted(d, key=d.get, reverse=True):
                if i < len(self.final_spec):
                    self.final_spec[i].set(k)
                    self.final_mole_frac[i].set(self.fmt.format(d[k]))
                i += 1

    def change_phi(self, *args):
        """Change the equivalence ratio by changing fuel content,
        holding diluent/o2 ratio constant.
        """
        if self.phi.get() == '':
            return
        if self.phi.get()[-1] == 0:
            # Don't update if last character is 0.
            return
        try:
            new_phi = float(self.phi.get())
        except ValueError:
            return

        try:
            phi, fuels = equivalence_ratio(self.gas, self.gas_list)
        except AttributeError:
            tkinter.messagebox.showerror('Error', 'Must calculate initial'
                                         'conditions before altering'
                                         'equivalence ratio')
            return
        if phi == 0 or new_phi == 0:
            return
        if new_phi == phi:
            return

        for i in range(len(self.gas_list)):
            row = self.gas_list[i]
            if row[0] in fuels:
                row[1] *= new_phi / phi
            self.species[i].set(row[0])
            self.species_vals[i].set(self.fmt.format(row[1]))

        # Reset variables before re-calculating initial state
        if self.initial_type == 'TP':
            self.H.set(0)
            self.S.set(0)
            self.D.set(0)
        elif self.initial_type == 'HP':
            self.T.set(0)
            self.S.set(0)
            self.D.set(0)
        elif self.initial_type == 'SP':
            self.T.set(0)
            self.H.set(0)
            self.D.set(0)
        elif self.initial_type == 'SD':
            self.T.set(0)
            self.P.set(0)
            self.H.set(0)
        self.calculate_final()  # Re-calculate


def equivalence_ratio(gas, gas_list):
    """ Calculate the equivalence ratio of the mixture

    input: cantera Solution, mixture components in list format
    Returns the equivalence ratio and a list of hydrocarbons
    """
    alpha = 0
    mol_O2 = 1  # Initialize
    fuels = []
    for row in gas_list:
        species = row[0]
        C = gas.n_atoms(species, 'C')
        H = gas.n_atoms(species, 'H')
        O = gas.n_atoms(species, 'O')
        if C == 1 and O == 2 and H == 0:
            continue  # This is CO2, move on
        if C == 0 and O == 1 and H == 2:
            continue  # This is H2O, move on
        if C == 0 and O == 2 and H == 0:
            mol_O2 = row[1]
            continue  # This is O2, move on
        if C + H + O == 0:
            continue  # This is inert

        fuels.append(row[0])
        alpha += (C + H/4 - O/2) * row[1]
    return (alpha / mol_O2, fuels)

if __name__ == "__main__":
    root = tkinter.Tk()
    GasProps(root)
    root.mainloop()
