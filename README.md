# NAMD_Config_Generator_GUI
A Python GUI for generating NAMD configuration files with ease and accuracy
This tool is designed to accelerate the setup process for molecular dynamics (MD) simulations by automating tedious tasks and providing a clear, user-friendly interface.

![GUI Screenshot](screenshot.png)

## Overview

Setting up a NAMD simulation involves creating a detailed configuration file with dozens of parameters. This process is often done manually, which is time-consuming and prone to human error. The Smart NAMD Configuration Generator solves this problem by:

-   Providing an intuitive graphical user interface (GUI) to manage all parameters.
-   **Automating the calculation** of periodic boundary conditions and PME grid sizes directly from your PDB file.
-   Allowing users to save and load entire configuration profiles, ensuring reproducibility and speed for future simulations.

The application is built with Python, using `ttkbootstrap` for a modern look and feel.


![image](https://github.com/user-attachments/assets/0262be4d-d165-438b-8bdf-5d46bceaa537)

## Key Features

-   **Automated Cell & PME Grid Calculation:** Reads your input PDB file and automatically calculates and populates the `cellBasisVector`, `cellOrigin`, and `PMEGridSize` fields.
-   **Advanced Filtering:** Fine-tune the cell calculation by filtering atoms based on segment, residue type, atom type, or index.
-   **Intuitive Tabbed Interface:** Parameters are logically grouped into sections like Input/Output, Dynamics, Thermostat/Barostat, PME/PBC, and Simulation Protocol.
-   **Save & Load Settings:** Save your entire GUI configuration to a JSON file and load it for later use.
-   **Interactive File Browsing:** "Browse..." buttons for all file path inputs.
-   **Helpful Tooltips:** Hover over any input field to get a brief explanation of the NAMD parameter.
-   **Scrollable Interface:** The GUI is fully scrollable, ensuring all features are accessible even on smaller screens.

## Installation

To get started, you'll need Python 3. Follow these steps to set up the application.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/masoudrezaei/NAMD_Config_Generator_GUI.git
    cd NAMD_Config_Generator_GUI
    ```

2.  **Install the required packages:**
    Install the only external dependency using pip:
    
    ```bash
    pip install ttkbootstrap
    ```
## Acknowledgements
The initial concept for a GUI-based NAMD configuration tool was inspired by the work of shmoe6 on the [NAMD-Config-Generator repository](https://github.com/shmoe6/NAMD-Config-Generator). This project builds on that idea with added a modernized toolkit.

## How to Run

Once the requirement is installed, you can run the application with a single command:

```bash
python namd_conf_gui.py

