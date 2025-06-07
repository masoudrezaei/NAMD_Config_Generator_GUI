import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import ttkbootstrap as bstrap
from ttkbootstrap.tooltip import ToolTip
import json
import webbrowser
import math

class ScrollableFrame(ttk.Frame):
    """A reusable scrollable frame widget."""
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        
        # Create a canvas to act as the viewport
        canvas = tk.Canvas(self, borderwidth=0, background=self.master.cget('background'))
        
        # Create a scrollbar and link it to the canvas
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        
        # This is the frame that will hold the actual content and be scrolled
        self.scrollable_frame = ttk.Frame(canvas)
        
        # Update the scrollregion of the canvas whenever the size of the inner frame changes
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )
        
        # Place the inner frame within the canvas window
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        
        # Configure the canvas to use the scrollbar's command
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Pack the canvas and scrollbar into the main frame
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

class NamdConfigurator(bstrap.Window):
    """
    A professional GUI for generating NAMD configuration files with
    automated cell dimension calculation and a scrollable interface.
    """
    def __init__(self):
        super().__init__(themename="litera")
        self.title("Smart NAMD Configuration Generator")
        # Set a reasonable default size; the scrollbar will handle overflow.
        self.geometry("850x750")
        self.minsize(600, 500)

        # --- Data Storage using tkinter variables ---
        self.vars = {
            # Input
            "coordinates": tk.StringVar(value="complex.pdb"),
            "structure": tk.StringVar(value="complex.psf"),
            "parameters": tk.StringVar(value="param.prm"),
            "paratypecharmm": tk.BooleanVar(value=True),
            "bincoordinates": tk.StringVar(),
            "binvelocities": tk.StringVar(),
            # Output
            "output": tk.StringVar(value="decoy"),
            "dcdfreq": tk.IntVar(value=5000),
            "xstFreq": tk.IntVar(value=5000),
            "binaryoutput": tk.BooleanVar(value=True),
            "binaryrestart": tk.BooleanVar(value=True),
            "outputEnergies": tk.IntVar(value=1000),
            "restartFreq": tk.IntVar(value=5000),
            # Basic Dynamics
            "exclude": tk.StringVar(value="scaled1-4"),
            "scaling1-4": tk.DoubleVar(value=1.0),
            "com_motion": tk.BooleanVar(value=False),
            "dielectric": tk.DoubleVar(value=1.0),
            # Simulation Space
            "switching": tk.BooleanVar(value=True),
            "switchdist": tk.DoubleVar(value=9.0),
            "cutoff": tk.DoubleVar(value=10.0),
            "pairlistdist": tk.DoubleVar(value=11.5),
            # Timestepping
            "firsttimestep": tk.IntVar(value=0),
            "timestep": tk.DoubleVar(value=2.0),
            "stepspercycle": tk.IntVar(value=20),
            # Temperature
            "temperature": tk.DoubleVar(value=310.0),
            # Langevin
            "langevin": tk.BooleanVar(value=True),
            "langevinDamping": tk.DoubleVar(value=1.0),
            # PME
            "pme": tk.BooleanVar(value=True),
            "pme_grid_x": tk.IntVar(value=0),
            "pme_grid_y": tk.IntVar(value=0),
            "pme_grid_z": tk.IntVar(value=0),
            # Pressure
            "useGroupPressure": tk.BooleanVar(value=True),
            "langevinPiston": tk.BooleanVar(value=True),
            "langevinPistonTarget": tk.DoubleVar(value=1.01325),
            "langevinPistonPeriod": tk.DoubleVar(value=150.0),
            "langevinPistonDecay": tk.DoubleVar(value=90.0),
            # PBC
            "cell_vec1_x": tk.DoubleVar(value=0.0), "cell_vec1_y": tk.DoubleVar(value=0.0), "cell_vec1_z": tk.DoubleVar(value=0.0),
            "cell_vec2_x": tk.DoubleVar(value=0.0), "cell_vec2_y": tk.DoubleVar(value=0.0), "cell_vec2_z": tk.DoubleVar(value=0.0),
            "cell_vec3_x": tk.DoubleVar(value=0.0), "cell_vec3_y": tk.DoubleVar(value=0.0), "cell_vec3_z": tk.DoubleVar(value=0.0),
            "cell_origin_x": tk.DoubleVar(value=0.0), "cell_origin_y": tk.DoubleVar(value=0.0), "cell_origin_z": tk.DoubleVar(value=0.0),
            "wrapWater": tk.BooleanVar(value=True),
            "wrapNearest": tk.BooleanVar(value=False),
            "wrapAll": tk.BooleanVar(value=True),
            # Script
            "minimize_steps": tk.IntVar(value=1000),
            "run_steps": tk.IntVar(value=500000),
            # Cell Calculation Filters
            "calc_segment": tk.StringVar(value="all"),
            "calc_residue": tk.StringVar(value="all"),
            "calc_atomtype": tk.StringVar(value="all"),
            "calc_firstatom": tk.StringVar(value="all"),
            "calc_lastatom": tk.StringVar(value="all"),
        }

        self._create_menu()
        self._create_widgets()
        self._create_status_bar()

    # --- Core Logic: Cell Dimension Calculation ---
    def _calculate_cell_dimensions(self):
        pdbfile = self.vars["coordinates"].get()
        if not pdbfile:
            messagebox.showerror("Error", "Please specify a Coordinates (PDB) file first on the 'Input / Output' tab.")
            return None
        segment = self.vars["calc_segment"].get()
        residue = self.vars["calc_residue"].get()
        atomtype = self.vars["calc_atomtype"].get()
        firstatom = self.vars["calc_firstatom"].get()
        lastatom = self.vars["calc_lastatom"].get()
        try:
            with open(pdbfile, 'r') as f:
                file_contents = f.read().splitlines()
        except FileNotFoundError:
            messagebox.showerror("File Not Found", f"The file '{pdbfile}' was not found.")
            return None
        except Exception as e:
            messagebox.showerror("Error", f"Could not read the file '{pdbfile}':\n{e}")
            return None

        natoms, atoms_considered = 0, 0
        xmin, ymin, zmin = float('inf'), float('inf'), float('inf')
        xmax, ymax, zmax = float('-inf'), float('-inf'), float('-inf')

        for line in file_contents:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                natoms += 1
                consider = True
                if firstatom.isdigit() and natoms < int(firstatom): consider = False
                if lastatom.isdigit() and natoms > int(lastatom): consider = False
                if segment != "all" and line[72:76].strip() != segment: consider = False
                if residue != "all" and line[17:21].strip() != residue: consider = False
                if atomtype != "all" and line[12:16].strip() != atomtype: consider = False
                if consider:
                    try:
                        atoms_considered += 1
                        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                        xmin, ymin, zmin = min(xmin, x), min(ymin, y), min(zmin, z)
                        xmax, ymax, zmax = max(xmax, x), max(ymax, y), max(zmax, z)
                    except (ValueError, IndexError): continue
        if atoms_considered == 0:
            messagebox.showwarning("Warning", "No atoms matched the specified filters. Cell dimensions were not calculated.")
            return None
        return {
            "vec1": (xmax - xmin, 0, 0), "vec2": (0, ymax - ymin, 0), "vec3": (0, 0, zmax - zmin),
            "origin": ((xmax + xmin) / 2, (ymax + ymin) / 2, (zmax + zmin) / 2),
            "atoms_found": atoms_considered
        }

    # --- GUI Button Callback ---
    def _on_calculate_cell(self):
        self.status_var.set("Calculating cell dimensions...")
        self.update_idletasks()
        results = self._calculate_cell_dimensions()
        if results:
            v1x, _, _ = results["vec1"]; _, v2y, _ = results["vec2"]; _, _, v3z = results["vec3"]
            ox, oy, oz = results["origin"]
            self.vars["cell_vec1_x"].set(round(v1x, 5)); self.vars["cell_vec2_y"].set(round(v2y, 5)); self.vars["cell_vec3_z"].set(round(v3z, 5))
            self.vars["cell_origin_x"].set(round(ox, 5)); self.vars["cell_origin_y"].set(round(oy, 5)); self.vars["cell_origin_z"].set(round(oz, 5))
            self.vars["pme_grid_x"].set(math.ceil(v1x)); self.vars["pme_grid_y"].set(math.ceil(v2y)); self.vars["pme_grid_z"].set(math.ceil(v3z))
            msg = f"Calculation complete. Found {results['atoms_found']} matching atoms. GUI fields updated."
            self.status_var.set(msg)
            messagebox.showinfo("Success", msg)
        else:
            self.status_var.set("Calculation failed or was cancelled.")

    # --- GUI Construction ---
    def _create_menu(self):
        menu_bar = tk.Menu(self)
        self.config(menu=menu_bar)
        file_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Generate NAMD Config...", command=self.generate_config)
        file_menu.add_separator()
        file_menu.add_command(label="Save GUI Settings...", command=self.save_settings)
        file_menu.add_command(label="Load GUI Settings...", command=self.load_settings)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.quit)
        help_menu = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="NAMD Documentation", command=lambda: webbrowser.open_new(r"https://www.ks.uiuc.edu/Research/namd/2.14/ug/"))
        help_menu.add_command(label="About", command=lambda: messagebox.showinfo("About", "Smart NAMD Configuration Generator\n\nIntegrates automated cell dimension calculation."))

    def _create_widgets(self):
        # Create the main scrollable container
        scroll_container = ScrollableFrame(self)
        scroll_container.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        # Get the inner frame where all content will go
        main_frame = scroll_container.scrollable_frame
        
        notebook = ttk.Notebook(main_frame)
        notebook.pack(fill=tk.BOTH, expand=True, pady=5)
        
        io_tab = self._create_tab(notebook, "Input / Output")
        dyn_tab = self._create_tab(notebook, "Dynamics & Timestep")
        thermo_tab = self._create_tab(notebook, "Thermostat / Barostat")
        pbc_tab = self._create_tab(notebook, "PME / PBC")
        script_tab = self._create_tab(notebook, "Simulation Protocol")
        
        self._populate_io_tab(io_tab)
        self._populate_dyn_tab(dyn_tab)
        self._populate_thermo_tab(thermo_tab)
        self._populate_pbc_tab(pbc_tab)
        self._populate_script_tab(script_tab)
        
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(fill=tk.X, pady=(10, 0))
        bstrap.Button(button_frame, text="Generate NAMD Config File", command=self.generate_config, bootstyle="success").pack(side=tk.RIGHT, padx=5)
        bstrap.Button(button_frame, text="Quit", command=self.quit, bootstyle="danger-outline").pack(side=tk.RIGHT)

    def _create_tab(self, notebook, text):
        frame = ttk.Frame(notebook, padding="10")
        notebook.add(frame, text=text)
        return frame

    def _create_entry_row(self, parent, label_text, var, tooltip_text, is_file=False, file_types=None):
        row_frame = ttk.Frame(parent)
        row_frame.pack(fill=tk.X, pady=4)
        ttk.Label(row_frame, text=label_text, width=20).pack(side=tk.LEFT, anchor="w")
        entry = ttk.Entry(row_frame, textvariable=var)
        entry.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
        ToolTip(entry, text=tooltip_text)
        if is_file:
            browse_btn = ttk.Button(row_frame, text="Browse...", command=lambda v=var, ft=file_types: self._browse_file(v, ft))
            browse_btn.pack(side=tk.LEFT)
        return row_frame

    def _create_check_row(self, parent, label_text, var, tooltip_text):
        check = bstrap.Checkbutton(parent, text=label_text, variable=var, bootstyle="primary-round-toggle")
        check.pack(anchor="w", pady=4)
        ToolTip(check, text=tooltip_text)

    def _populate_pbc_tab(self, tab):
        calc_frame = ttk.LabelFrame(tab, text="Cell Dimension Calculator", padding="10")
        calc_frame.pack(fill=tk.X, pady=5)
        ttk.Label(calc_frame, text="Use 'all' to ignore a filter. Calculation is based on the PDB file from the Input tab.").pack(anchor='w', pady=(0, 5))
        filter_grid = ttk.Frame(calc_frame)
        filter_grid.pack(fill=tk.X, expand=True)
        self._create_entry_row(filter_grid, "Segment", self.vars["calc_segment"], "Filter by segment name (e.g., 'PROT').")
        self._create_entry_row(filter_grid, "Residue Type", self.vars["calc_residue"], "Filter by residue name (e.g., 'LYS').")
        self._create_entry_row(filter_grid, "Atom Type", self.vars["calc_atomtype"], "Filter by atom name (e.g., 'CA').")
        self._create_entry_row(filter_grid, "First Atom Index", self.vars["calc_firstatom"], "Consider atoms with index >= this value.")
        self._create_entry_row(filter_grid, "Last Atom Index", self.vars["calc_lastatom"], "Consider atoms with index <= this value.")
        calc_btn = bstrap.Button(calc_frame, text="Auto-calculate Cell & PME Grid", command=self._on_calculate_cell, bootstyle="info")
        calc_btn.pack(pady=10)
        ToolTip(calc_btn, "Calculates dimensions from the PDB file and fills the fields below.")
        pme_frame = ttk.LabelFrame(tab, text="Particle Mesh Ewald (PME)", padding="10")
        pme_frame.pack(fill=tk.X, pady=5)
        self._create_check_row(pme_frame, "Use PME", self.vars["pme"], "Enable PME for long-range electrostatic calculations.")
        grid_frame = ttk.Frame(pme_frame)
        grid_frame.pack(fill=tk.X, pady=4)
        ttk.Label(grid_frame, text="PME Grid Size", width=20).pack(side=tk.LEFT, anchor="w")
        for dim in ['x', 'y', 'z']:
            ttk.Label(grid_frame, text=f"{dim.upper()}:").pack(side=tk.LEFT, padx=(10, 2))
            entry = ttk.Entry(grid_frame, textvariable=self.vars[f"pme_grid_{dim}"], width=8)
            entry.pack(side=tk.LEFT)
            ToolTip(entry, f"Number of grid points in the {dim} dimension for PME. Can be auto-calculated.")
        pbc_frame = ttk.LabelFrame(tab, text="Periodic Boundary Conditions (Manual or Auto-calculated)", padding="10")
        pbc_frame.pack(fill=tk.X, pady=5)
        for i in range(1, 4):
            vec_frame = ttk.Frame(pbc_frame)
            vec_frame.pack(fill=tk.X, pady=2)
            ttk.Label(vec_frame, text=f"Cell Basis Vector {i}", width=20).pack(side=tk.LEFT, anchor="w")
            for dim in ['x', 'y', 'z']:
                entry = ttk.Entry(vec_frame, textvariable=self.vars[f"cell_vec{i}_{dim}"], width=12)
                entry.pack(side=tk.LEFT, padx=5)
                ToolTip(entry, f"Vector {i}, {dim}-component")
        origin_frame = ttk.Frame(pbc_frame)
        origin_frame.pack(fill=tk.X, pady=2)
        ttk.Label(origin_frame, text="Cell Origin", width=20).pack(side=tk.LEFT, anchor="w")
        for dim in ['x', 'y', 'z']:
            entry = ttk.Entry(origin_frame, textvariable=self.vars[f"cell_origin_{dim}"], width=12)
            entry.pack(side=tk.LEFT, padx=5)
            ToolTip(entry, f"Origin {dim}-component")
        wrap_frame = ttk.LabelFrame(tab, text="Atom Wrapping", padding="10")
        wrap_frame.pack(fill=tk.X, pady=5)
        self._create_check_row(wrap_frame, "Wrap Water", self.vars["wrapWater"], "Wrap water molecules back into the primary cell.")
        self._create_check_row(wrap_frame, "Wrap All", self.vars["wrapAll"], "Wrap all atoms back into the primary cell.")
        self._create_check_row(wrap_frame, "Wrap Nearest", self.vars["wrapNearest"], "Use for non-rectangular cells.")

    def _populate_io_tab(self, tab):
        input_frame = ttk.LabelFrame(tab, text="Input Files", padding="10")
        input_frame.pack(fill=tk.X, pady=5)
        self._create_entry_row(input_frame, "Coordinates", self.vars["coordinates"], "PDB file. Used for auto-calculating cell dimensions.", True, [("PDB files", "*.pdb"), ("All files", "*.*")])
        self._create_entry_row(input_frame, "Structure", self.vars["structure"], "PSF file describing molecular structure.", True, [("PSF files", "*.psf"), ("All files", "*.*")])
        self._create_entry_row(input_frame, "Parameters", self.vars["parameters"], "Parameter file for the force field.", True, [("Parameter files", "*.prm *.par"), ("All files", "*.*")])
        self._create_check_row(input_frame, "Parameters are CHARMM type", self.vars["paratypecharmm"], "Specifies that the parameter file uses the CHARMM format (X-PLOR otherwise).")
        restart_frame = ttk.LabelFrame(tab, text="Restart Files (Optional)", padding="10")
        restart_frame.pack(fill=tk.X, pady=5)
        self._create_entry_row(restart_frame, "Binary Coordinates", self.vars["bincoordinates"], "Restart file with binary coordinates (.coor).", True, [("Restart Coordinates", "*.coor"), ("All files", "*.*")])
        self._create_entry_row(restart_frame, "Binary Velocities", self.vars["binvelocities"], "Restart file with binary velocities (.vel).", True, [("Restart Velocities", "*.vel"), ("All files", "*.*")])
        output_frame = ttk.LabelFrame(tab, text="Output Files", padding="10")
        output_frame.pack(fill=tk.X, pady=5)
        self._create_entry_row(output_frame, "Output Name", self.vars["output"], "Base name for all output files (e.g., 'decoy').")
        self._create_entry_row(output_frame, "DCD Frequency", self.vars["dcdfreq"], "Frequency in steps to write to the DCD trajectory file.")
        self._create_entry_row(output_frame, "XST Frequency", self.vars["xstFreq"], "Frequency in steps to write cell information to the XST file.")
        self._create_entry_row(output_frame, "Energy Frequency", self.vars["outputEnergies"], "Frequency in steps to print energies to the log file.")
        self._create_entry_row(output_frame, "Restart Frequency", self.vars["restartFreq"], "Frequency in steps to write restart files.")
        self._create_check_row(output_frame, "Use Binary Output", self.vars["binaryoutput"], "Generate smaller, binary output files.")
        self._create_check_row(output_frame, "Use Binary Restart", self.vars["binaryrestart"], "Generate smaller, binary restart files.")

    def _populate_dyn_tab(self, tab):
        dynamics_frame = ttk.LabelFrame(tab, text="Basic Dynamics", padding="10")
        dynamics_frame.pack(fill=tk.X, pady=5)
        row = ttk.Frame(dynamics_frame)
        row.pack(fill=tk.X, pady=4)
        ttk.Label(row, text="Exclusions", width=20).pack(side=tk.LEFT, anchor="w")
        combo = ttk.Combobox(row, textvariable=self.vars["exclude"], values=["none", "scaled1-4"])
        combo.pack(side=tk.LEFT, expand=True, fill=tk.X, padx=5)
        ToolTip(combo, "Policy for excluding non-bonded interactions.")
        self._create_entry_row(dynamics_frame, "1-4 Scaling", self.vars["scaling1-4"], "Scaling factor for 1-4 electrostatic and VdW interactions.")
        self._create_entry_row(dynamics_frame, "Dielectric Constant", self.vars["dielectric"], "Dielectric constant for electrostatic calculations.")
        self._create_check_row(dynamics_frame, "Remove COM Motion", self.vars["com_motion"], "Periodically remove center-of-mass motion (set to 'no' for constant pressure).")
        space_frame = ttk.LabelFrame(tab, text="Simulation Space Partitioning", padding="10")
        space_frame.pack(fill=tk.X, pady=5)
        self._create_check_row(space_frame, "Use Switching", self.vars["switching"], "Use a switching function to smoothly turn off interactions at the cutoff.")
        self._create_entry_row(space_frame, "Switch Distance", self.vars["switchdist"], "Distance at which the switching function begins.")
        self._create_entry_row(space_frame, "Cutoff Distance", self.vars["cutoff"], "The non-bonded interaction cutoff distance.")
        self._create_entry_row(space_frame, "Pairlist Distance", self.vars["pairlistdist"], "Distance for generating the pairlist, should be > cutoff.")
        time_frame = ttk.LabelFrame(tab, text="Multiple Timestepping", padding="10")
        time_frame.pack(fill=tk.X, pady=5)
        self._create_entry_row(time_frame, "First Timestep", self.vars["firsttimestep"], "Initial timestep number for the simulation.")
        self._create_entry_row(time_frame, "Timestep (fs)", self.vars["timestep"], "Simulation timestep in femtoseconds.")
        self._create_entry_row(time_frame, "Steps per Cycle", self.vars["stepspercycle"], "Number of timesteps in a cycle. Long-range forces are evaluated once per cycle.")

    def _populate_thermo_tab(self, tab):
        temp_frame = ttk.LabelFrame(tab, text="Temperature Control (Thermostat)", padding="10")
        temp_frame.pack(fill=tk.X, pady=5)
        self._create_entry_row(temp_frame, "Target Temperature (K)", self.vars["temperature"], "The target temperature for the system in Kelvin.")
        self._create_check_row(temp_frame, "Use Langevin Dynamics", self.vars["langevin"], "Enable the Langevin thermostat for temperature control.")
        self._create_entry_row(temp_frame, "Langevin Damping (1/ps)", self.vars["langevinDamping"], "Langevin dynamics damping coefficient.")
        pressure_frame = ttk.LabelFrame(tab, text="Pressure Control (Barostat)", padding="10")
        pressure_frame.pack(fill=tk.X, pady=5)
        self._create_check_row(pressure_frame, "Use Langevin Piston", self.vars["langevinPiston"], "Enable the Langevin piston for pressure control.")
        self._create_entry_row(pressure_frame, "Piston Target (bar)", self.vars["langevinPistonTarget"], "Target pressure in bar (will be converted to NAMD units).")
        self._create_entry_row(pressure_frame, "Piston Period (fs)", self.vars["langevinPistonPeriod"], "Oscillation period for the Langevin piston.")
        self._create_entry_row(pressure_frame, "Piston Decay (fs)", self.vars["langevinPistonDecay"], "Damping time scale for the Langevin piston.")
        self._create_check_row(pressure_frame, "Use Group Pressure", self.vars["useGroupPressure"], "Calculate pressure based on atomic groups for better stability.")

    def _populate_script_tab(self, tab):
        script_frame = ttk.LabelFrame(tab, text="Simulation Protocol", padding="10")
        script_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        ttk.Label(script_frame, text="This tab defines the sequence of simulation steps.").pack(anchor='w', pady=(0,10))
        self._create_entry_row(script_frame, "Minimization Steps", self.vars["minimize_steps"], "Number of conjugate gradient energy minimization steps to perform first.")
        self._create_entry_row(script_frame, "MD Run Steps", self.vars["run_steps"], "Number of molecular dynamics steps to run after minimization/reinitialization.")
        total_time = (self.vars["run_steps"].get() * self.vars["timestep"].get()) / 1000
        time_label_text = f"This corresponds to a simulation time of {total_time:.2f} ps ({total_time/1000:.3f} ns)."
        self.time_label = ttk.Label(script_frame, text=time_label_text)
        self.time_label.pack(anchor='w', pady=5, padx=5)
        self.vars["run_steps"].trace_add("write", self._update_sim_time)
        self.vars["timestep"].trace_add("write", self._update_sim_time)

    def _update_sim_time(self, *args):
        try:
            total_time = (self.vars["run_steps"].get() * self.vars["timestep"].get()) / 1000
            self.time_label.config(text=f"This corresponds to a simulation time of {total_time:.2f} ps ({total_time/1000:.3f} ns).")
        except (tk.TclError, ValueError):
            self.time_label.config(text="Invalid number format.")

    # --- Utility and File I/O Methods ---
    def _browse_file(self, var, file_types):
        filename = filedialog.askopenfilename(filetypes=file_types)
        if filename: var.set(filename)

    def _create_status_bar(self):
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(self, textvariable=self.status_var, anchor=tk.W, relief=tk.SUNKEN, padding=5)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)

    def save_settings(self):
        filepath = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json"), ("All files", "*.*")], title="Save GUI Settings")
        if not filepath: return
        settings_to_save = {key: var.get() for key, var in self.vars.items()}
        try:
            with open(filepath, 'w') as f: json.dump(settings_to_save, f, indent=4)
            self.status_var.set(f"Settings saved to {filepath}")
        except Exception as e: messagebox.showerror("Error", f"Failed to save settings: {e}")

    def load_settings(self):
        filepath = filedialog.askopenfilename(filetypes=[("JSON files", "*.json"), ("All files", "*.*")], title="Load GUI Settings")
        if not filepath: return
        try:
            with open(filepath, 'r') as f: loaded_settings = json.load(f)
            for key, value in loaded_settings.items():
                if key in self.vars: self.vars[key].set(value)
            self.status_var.set(f"Settings loaded from {filepath}")
        except Exception as e: messagebox.showerror("Error", f"Failed to load settings: {e}")

    def generate_config(self):
        filepath = filedialog.asksaveasfilename(defaultextension=".namd", filetypes=[("NAMD Config files", "*.namd"), ("All files", "*.*")], title="Save NAMD Configuration File")
        if not filepath: return
        on_off = lambda b: "on" if b else "off"; yes_no = lambda b: "yes" if b else "no"
        pressure_atm = self.vars['langevinPistonTarget'].get() * 0.986923
        config_str = f"""# NAMD Configuration File generated by Smart Python GUI
# ----- Input Files -----
coordinates\t\t{self.vars['coordinates'].get()}
structure\t\t{self.vars['structure'].get()}
parameters\t\t{self.vars['parameters'].get()}
paratypecharmm\t\t{yes_no(self.vars['paratypecharmm'].get())}
# bincoordinates\t\t{self.vars['bincoordinates'].get() or 'previous.restart.coor'}
# binvelocities\t\t{self.vars['binvelocities'].get() or 'previous.restart.vel'}
# ----- Output -----
set output\t\t{self.vars['output'].get()}
outputname\t\t${{output}}
dcdfile\t\t\t${{output}}.dcd
xstFile\t\t\t${{output}}.xst
dcdfreq\t\t\t{self.vars['dcdfreq'].get()}
xstFreq\t\t\t{self.vars['xstFreq'].get()}
binaryoutput\t\t{yes_no(self.vars['binaryoutput'].get())}
binaryrestart\t\t{yes_no(self.vars['binaryrestart'].get())}
outputEnergies\t\t{self.vars['outputEnergies'].get()}
restartFreq\t\t{self.vars['restartFreq'].get()}
# ----- Basic Dynamics -----
exclude\t\t\t{self.vars['exclude'].get()}
1-4scaling\t\t{self.vars['scaling1-4'].get()}
COMmotion\t\t{yes_no(self.vars['com_motion'].get())}
dielectric\t\t{self.vars['dielectric'].get()}
# ----- Simulation Space Partitioning -----
switching\t\t{on_off(self.vars['switching'].get())}
switchdist\t\t{self.vars['switchdist'].get()}
cutoff\t\t\t{self.vars['cutoff'].get()}
pairlistdist\t\t{self.vars['pairlistdist'].get()}
# ----- Multiple Timestepping -----
firsttimestep\t\t{self.vars['firsttimestep'].get()}
timestep\t\t{self.vars['timestep'].get()}
stepspercycle\t\t{self.vars['stepspercycle'].get()}
# ----- Temperature Control -----
set temperature\t\t{self.vars['temperature'].get()}
temperature\t\t$temperature
langevin\t\t{on_off(self.vars['langevin'].get())}
langevinDamping\t\t{self.vars['langevinDamping'].get()}
langevinTemp\t\t$temperature
# ----- Pressure and Periodic Boundary Conditions -----
PME\t\t\t{on_off(self.vars['pme'].get())}
PMEGridSizeX\t\t{self.vars['pme_grid_x'].get()}
PMEGridSizeY\t\t{self.vars['pme_grid_y'].get()}
PMEGridSizeZ\t\t{self.vars['pme_grid_z'].get()}
useGroupPressure\t{yes_no(self.vars['useGroupPressure'].get())}
LangevinPiston\t\t{on_off(self.vars['langevinPiston'].get())}
LangevinPistonTarget\t{pressure_atm:.5f} ;# {self.vars['langevinPistonTarget'].get()} bar
LangevinPistonPeriod\t{self.vars['langevinPistonPeriod'].get()}
LangevinPistonDecay\t{self.vars['langevinPistonDecay'].get()}
LangevinPistonTemp\t$temperature
cellBasisVector1\t{self.vars['cell_vec1_x'].get()} {self.vars['cell_vec1_y'].get()} {self.vars['cell_vec1_z'].get()}
cellBasisVector2\t{self.vars['cell_vec2_x'].get()} {self.vars['cell_vec2_y'].get()} {self.vars['cell_vec2_z'].get()}
cellBasisVector3\t{self.vars['cell_vec3_x'].get()} {self.vars['cell_vec3_y'].get()} {self.vars['cell_vec3_z'].get()}
cellOrigin\t\t{self.vars['cell_origin_x'].get()} {self.vars['cell_origin_y'].get()} {self.vars['cell_origin_z'].get()}
wrapWater\t\t{on_off(self.vars['wrapWater'].get())}
wrapNearest\t\t{on_off(self.vars['wrapNearest'].get())}
wrapAll\t\t\t{on_off(self.vars['wrapAll'].get())}
# ----- Simulation Protocol -----
minimize\t\t{self.vars['minimize_steps'].get()}
reinitvels\t\t$temperature
run\t\t\t{self.vars['run_steps'].get()};
END
"""
        try:
            with open(filepath, 'w') as f: f.write(config_str)
            messagebox.showinfo("Success", f"NAMD configuration file saved successfully to:\n{filepath}")
            self.status_var.set(f"Config file saved: {filepath}")
        except Exception as e: messagebox.showerror("Error", f"Could not save file: {e}")

if __name__ == "__main__":
    app = NamdConfigurator()
    app.mainloop()