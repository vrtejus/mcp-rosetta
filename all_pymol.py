PYMOL_COMMANDS = {
    # MOLECULAR VISUALIZATION
    "show": {
        "description": "Shows a representation for the specified selection",
        "pattern": r"^show\s+([\w.]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "representation", "required": True, "options": [
                "lines", "sticks", "spheres", "surface", "mesh", "dots", 
                "ribbon", "cartoon", "labels", "nonbonded", "nb_spheres", 
                "ellipsoids", "volume", "slice", "extent", "dots_as_spheres", 
                "cell", "cgo", "everything", "dashes", "angles", "dihedrals", 
                "licorice", "spheres", "putty"
            ]},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "hide": {
        "description": "Hides a representation for the specified selection",
        "pattern": r"^hide\s+([\w.]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "representation", "required": True, "options": [
                "lines", "sticks", "spheres", "surface", "mesh", "dots", 
                "ribbon", "cartoon", "labels", "nonbonded", "nb_spheres", 
                "ellipsoids", "volume", "slice", "extent", "dots_as_spheres", 
                "cell", "cgo", "everything", "dashes", "angles", "dihedrals", 
                "licorice", "spheres", "putty"
            ]},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "color": {
        "description": "Sets the color for the specified selection",
        "pattern": r"^color\s+([\w.]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "color", "required": True},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "as": {
        "description": "Shows one representation while hiding all others for the specified selection",
        "pattern": r"^as\s+([\w.]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "representation", "required": True},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True 
    },
    "set": {
        "description": "Sets a PyMOL setting to a specified value",
        "pattern": r"^set\s+([\w.]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "setting", "required": True},
            {"name": "value", "required": True},
            {"name": "selection", "required": False}
        ],
        "check_selection": False
    },
    "cartoon": {
        "description": "Sets the cartoon type for the specified selection",
        "pattern": r"^cartoon\s+([\w.]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "type", "required": True, "options": [
                "automatic", "loop", "rectangle", "oval", "tube", "arrow", "dumbbell", "putty"
            ]},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "spectrum": {
        "description": "Colors selection in a spectrum",
        "pattern": r"^spectrum\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "expression", "required": True},
            {"name": "palette", "required": False, "default": "rainbow"},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "label": {
        "description": "Adds labels to atoms in the selection",
        "pattern": r"^label\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "selection", "required": True},
            {"name": "expression", "required": False, "default": "name"}
        ],
        "check_selection": True
    },
    "distance": {
        "description": "Measures the distance between two selections",
        "pattern": r"^distance(?:\s+([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?$",
        "parameters": [
            {"name": "name", "required": False},
            {"name": "selection1", "required": False, "default": "(pk1)"},
            {"name": "selection2", "required": False, "default": "(pk2)"}
        ],
        "check_selection": True
    },
    "angle": {
        "description": "Measures the angle between three selections",
        "pattern": r"^angle(?:\s+([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?$",
        "parameters": [
            {"name": "name", "required": False},
            {"name": "selection1", "required": False, "default": "(pk1)"},
            {"name": "selection2", "required": False, "default": "(pk2)"},
            {"name": "selection3", "required": False, "default": "(pk3)"}
        ],
        "check_selection": True
    },
    "dihedral": {
        "description": "Measures the dihedral angle between four selections",
        "pattern": r"^dihedral(?:\s+([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?$",
        "parameters": [
            {"name": "name", "required": False},
            {"name": "selection1", "required": False, "default": "(pk1)"},
            {"name": "selection2", "required": False, "default": "(pk2)"},
            {"name": "selection3", "required": False, "default": "(pk3)"},
            {"name": "selection4", "required": False, "default": "(pk4)"}
        ],
        "check_selection": True
    },
    
    # VIEWING OPERATIONS
    "center": {
        "description": "Centers the view on a selection",
        "pattern": r"^center(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "orient": {
        "description": "Orients the view to align with principal axes of the selection",
        "pattern": r"^orient(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "zoom": {
        "description": "Zooms the view on a selection",
        "pattern": r"^zoom(?:\s+([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"},
            {"name": "buffer", "required": False, "default": "5"}
        ],
        "check_selection": True
    },
    "reset": {
        "description": "Resets the view, optionally resetting an object's matrix",
        "pattern": r"^reset(?:\s+(.+))?$",
        "parameters": [
            {"name": "object", "required": False}
        ],
        "check_selection": False
    },
    "turn": {
        "description": "Rotates the camera around an axis",
        "pattern": r"^turn\s+([xyz])(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "axis", "required": True, "options": ["x", "y", "z"]},
            {"name": "angle", "required": False, "default": "90"}
        ],
        "check_selection": False
    },
    "move": {
        "description": "Moves the camera along an axis",
        "pattern": r"^move\s+([xyz])(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "axis", "required": True, "options": ["x", "y", "z"]},
            {"name": "distance", "required": False, "default": "1"}
        ],
        "check_selection": False
    },
    "clip": {
        "description": "Adjusts the clipping planes",
        "pattern": r"^clip\s+([\w.]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "mode", "required": True, "options": ["near", "far", "slab", "atoms", "near_slab", "far_slab"]},
            {"name": "distance", "required": False, "default": "1"}
        ],
        "check_selection": False
    },
    
    # FILE OPERATIONS
    "load": {
        "description": "Loads a file into PyMOL",
        "pattern": r"^load\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "filename", "required": True},
            {"name": "object", "required": False},
            {"name": "options", "required": False}
        ],
        "check_selection": False
    },
    "fetch": {
        "description": "Fetches a structure from a database (e.g., PDB)",
        "pattern": r"^fetch\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "code", "required": True},
            {"name": "name", "required": False},
            {"name": "options", "required": False}
        ],
        "check_selection": False
    },
    "save": {
        "description": "Saves data to a file",
        "pattern": r"^save\s+([^,]+)(?:\s*,\s*(.+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "filename", "required": True},
            {"name": "selection", "required": False, "default": "all"},
            {"name": "state", "required": False, "default": "-1"}
        ],
        "check_selection": True
    },
    "png": {
        "description": "Saves a PNG image",
        "pattern": r"^png\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "filename", "required": True},
            {"name": "options", "required": False}
        ],
        "check_selection": False
    },
    
    # SELECTION OPERATIONS
    "select": {
        "description": "Creates a named selection",
        "pattern": r"^select\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "name", "required": True},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": False
    },
    "deselect": {
        "description": "Clears the current selection",
        "pattern": r"^deselect$",
        "parameters": [],
        "check_selection": False
    },
    
    # OBJECT MANIPULATION
    "create": {
        "description": "Creates a new object from a selection",
        "pattern": r"^create\s+([^,]+)(?:\s*,\s*(.+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "name", "required": True},
            {"name": "selection", "required": False, "default": "all"},
            {"name": "source_state", "required": False, "default": "1"}
        ],
        "check_selection": True
    },
    "extract": {
        "description": "Extracts a selection to a new object",
        "pattern": r"^extract\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "name", "required": True},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "delete": {
        "description": "Deletes objects or selections",
        "pattern": r"^delete\s+(.+)$",
        "parameters": [
            {"name": "name", "required": True}
        ],
        "check_selection": False
    },
    "remove": {
        "description": "Removes atoms in a selection",
        "pattern": r"^remove\s+(.+)$",
        "parameters": [
            {"name": "selection", "required": True}
        ],
        "check_selection": True
    },
    "align": {
        "description": "Aligns one selection to another",
        "pattern": r"^align\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "mobile", "required": True},
            {"name": "target", "required": False, "default": "all"},
            {"name": "options", "required": False}
        ],
        "check_selection": True
    },
    "super": {
        "description": "Superimposes one selection onto another",
        "pattern": r"^super\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "mobile", "required": True},
            {"name": "target", "required": False, "default": "all"},
            {"name": "options", "required": False}
        ],
        "check_selection": True
    },
    "intra_fit": {
        "description": "Fits all states within an object",
        "pattern": r"^intra_fit\s+(.+)$",
        "parameters": [
            {"name": "selection", "required": True}
        ],
        "check_selection": True
    },
    "intra_rms": {
        "description": "Calculates RMSD between states within an object",
        "pattern": r"^intra_rms\s+(.+)$",
        "parameters": [
            {"name": "selection", "required": True}
        ],
        "check_selection": True
    },
    
    # UTILITY AND MODIFICATION
    "alter": {
        "description": "Alters atomic properties in a selection",
        "pattern": r"^alter\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "selection", "required": True},
            {"name": "expression", "required": True}
        ],
        "check_selection": True
    },
    "alter_state": {
        "description": "Alters atomic coordinates in a state",
        "pattern": r"^alter_state\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "state", "required": True},
            {"name": "selection", "required": True},
            {"name": "expression", "required": True}
        ],
        "check_selection": True
    },
    "h_add": {
        "description": "Adds hydrogens to a selection",
        "pattern": r"^h_add(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "h_fill": {
        "description": "Adds hydrogens and adjusts valences",
        "pattern": r"^h_fill(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "bond": {
        "description": "Creates a bond between two atoms",
        "pattern": r"^bond\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "atom1", "required": True},
            {"name": "atom2", "required": True},
            {"name": "order", "required": False, "default": "1"}
        ],
        "check_selection": True
    },
    "unbond": {
        "description": "Removes a bond between two atoms",
        "pattern": r"^unbond\s+([^,]+)(?:\s*,\s*([^,]+))?$",
        "parameters": [
            {"name": "atom1", "required": True},
            {"name": "atom2", "required": True}
        ],
        "check_selection": True
    },
    "rebuild": {
        "description": "Regenerates all displayed geometry",
        "pattern": r"^rebuild(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": False
    },
    "refresh": {
        "description": "Refreshes the display",
        "pattern": r"^refresh$",
        "parameters": [],
        "check_selection": False
    },
    
    # UTILITY FUNCTIONS
    "util.cbc": {
        "description": "Colors by chain (Color By Chain)",
        "pattern": r"^util\.cbc(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbaw": {
        "description": "Colors by atom, white carbons (Color By Atom, White)",
        "pattern": r"^util\.cbaw(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbag": {
        "description": "Colors by atom, green carbons (Color By Atom, Green)",
        "pattern": r"^util\.cbag(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbac": {
        "description": "Colors by atom, cyan carbons (Color By Atom, Cyan)",
        "pattern": r"^util\.cbac(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbam": {
        "description": "Colors by atom, magenta carbons (Color By Atom, Magenta)",
        "pattern": r"^util\.cbam(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbay": {
        "description": "Colors by atom, yellow carbons (Color By Atom, Yellow)",
        "pattern": r"^util\.cbay(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbas": {
        "description": "Colors by atom, salmon carbons (Color By Atom, Salmon)",
        "pattern": r"^util\.cbas(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbab": {
        "description": "Colors by atom, slate carbons (Color By Atom, slateBLue)",
        "pattern": r"^util\.cbab(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbao": {
        "description": "Colors by atom, orange carbons (Color By Atom, Orange)",
        "pattern": r"^util\.cbao(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbap": {
        "description": "Colors by atom, purple carbons (Color By Atom, Purple)",
        "pattern": r"^util\.cbap(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.cbak": {
        "description": "Colors by atom, pink carbons (Color By Atom, pinK)",
        "pattern": r"^util\.cbak(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.chainbow": {
        "description": "Colors chains in rainbow gradient (CHAINs in rainBOW)",
        "pattern": r"^util\.chainbow(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.rainbow": {
        "description": "Colors residues in rainbow from N to C terminus",
        "pattern": r"^util\.rainbow(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.ss": {
        "description": "Colors by secondary structure",
        "pattern": r"^util\.ss(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.color_by_element": {
        "description": "Colors atoms by their element",
        "pattern": r"^util\.color_by_element(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "util.color_secondary": {
        "description": "Colors secondary structure elements",
        "pattern": r"^util\.color_secondary(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    
    # MOLECULAR DYNAMICS AND ANALYSIS
    "spheroid": {
        "description": "Displays atoms as smooth spheres",
        "pattern": r"^spheroid(?:\s+(.+))?$",
        "parameters": [
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "isomesh": {
        "description": "Creates a mesh isosurface",
        "pattern": r"^isomesh\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "name", "required": True},
            {"name": "map_object", "required": True},
            {"name": "level", "required": True},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "isosurface": {
        "description": "Creates a solid isosurface",
        "pattern": r"^isosurface\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "name", "required": True},
            {"name": "map_object", "required": True},
            {"name": "level", "required": True},
            {"name": "selection", "required": False, "default": "all"}
        ],
        "check_selection": True
    },
    "sculpt_activate": {
        "description": "Activates sculpting mode for an object",
        "pattern": r"^sculpt_activate\s+(.+)$",
        "parameters": [
            {"name": "object", "required": True}
        ],
        "check_selection": False
    },
    "sculpt_deactivate": {
        "description": "Deactivates sculpting mode for an object",
        "pattern": r"^sculpt_deactivate\s+(.+)$",
        "parameters": [
            {"name": "object", "required": True}
        ],
        "check_selection": False
    },
    "sculpt_iterate": {
        "description": "Performs sculpting iterations",
        "pattern": r"^sculpt_iterate\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "iterations", "required": True},
            {"name": "object", "required": False, "default": "all"}
        ],
        "check_selection": False
    },
    
    # SCENES AND MOVIES
    "scene": {
        "description": "Manages scenes for later recall",
        "pattern": r"^scene\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "key", "required": True},
            {"name": "action", "required": False, "default": "recall"}
        ],
        "check_selection": False
    },
    "scene_order": {
        "description": "Sets the order of scenes",
        "pattern": r"^scene_order\s+(.+)$",
        "parameters": [
            {"name": "scene_list", "required": True}
        ],
        "check_selection": False
    },
    "mset": {
        "description": "Defines a sequence of states for movie playback",
        "pattern": r"^mset\s+(.+)$",
        "parameters": [
            {"name": "specification", "required": True}
        ],
        "check_selection": False
    },
    "mplay": {
        "description": "Starts playing the movie",
        "pattern": r"^mplay$",
        "parameters": [],
        "check_selection": False
    },
    "mstop": {
        "description": "Stops the movie",
        "pattern": r"^mstop$",
        "parameters": [],
        "check_selection": False
    },
    "frame": {
        "description": "Sets or queries the current frame",
        "pattern": r"^frame(?:\s+(.+))?$",
        "parameters": [
            {"name": "frame_number", "required": False}
        ],
        "check_selection": False
    },
    "forward": {
        "description": "Advances one frame",
        "pattern": r"^forward$",
        "parameters": [],
        "check_selection": False
    },
    "backward": {
        "description": "Goes back one frame",
        "pattern": r"^backward$",
        "parameters": [],
        "check_selection": False
    },
    "rock": {
        "description": "Toggles a rocking animation",
        "pattern": r"^rock$",
        "parameters": [],
        "check_selection": False
    },
    
    # RENDERING
    "ray": {
        "description": "Performs ray-tracing",
        "pattern": r"^ray(?:\s+([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "width", "required": False},
            {"name": "height", "required": False}
        ],
        "check_selection": False
    },
    "draw": {
        "description": "Uses OpenGL renderer (faster but lower quality)",
        "pattern": r"^draw(?:\s+([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "width", "required": False},
            {"name": "height", "required": False}
        ],
        "check_selection": False
    },
    "mpng": {
        "description": "Saves a series of PNG images for movie frames",
        "pattern": r"^mpng\s+(.+)$",
        "parameters": [
            {"name": "prefix", "required": True}
        ],
        "check_selection": False
    },
    
    # CRYSTALLOGRAPHY
    "symexp": {
        "description": "Generates symmetry-related copies",
        "pattern": r"^symexp\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "prefix", "required": True},
            {"name": "selection", "required": True},
            {"name": "cutoff", "required": False, "default": "20"},
            {"name": "segi", "required": False}
        ],
        "check_selection": True
    },
    "symexp": {
        "description": "Generates symmetry-related copies",
        "pattern": r"^symexp\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "prefix", "required": True},
            {"name": "selection", "required": True},
            {"name": "cutoff", "required": False, "default": "20"},
            {"name": "segi", "required": False}
        ],
        "check_selection": True
    },
    "set_symmetry": {
        "description": "Sets symmetry parameters for an object",
        "pattern": r"^set_symmetry\s+([^,]+)(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?(?:\s*,\s*([^,]+))?$",
        "parameters": [
            {"name": "selection", "required": True},
            {"name": "a", "required": True},
            {"name": "b", "required": True},
            {"name": "c", "required": True},
            {"name": "alpha", "required": True},
            {"name": "beta", "required": True},
            {"name": "gamma", "required": True}
        ],
        "check_selection": True
    },
    
    # OTHER
    "fab": {
        "description": "Creates a peptide chain from a sequence",
        "pattern": r"^fab\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "sequence", "required": True},
            {"name": "options", "required": False}
        ],
        "check_selection": False
    },
    "fragment": {
        "description": "Loads a molecular fragment",
        "pattern": r"^fragment\s+(.+)$",
        "parameters": [
            {"name": "name", "required": True}
        ],
        "check_selection": False
    },
    "full_screen": {
        "description": "Toggles fullscreen mode",
        "pattern": r"^full_screen$",
        "parameters": [],
        "check_selection": False
    },
    "viewport": {
        "description": "Sets the viewport size",
        "pattern": r"^viewport\s+([^,]+)(?:\s*,\s*(.+))?$",
        "parameters": [
            {"name": "width", "required": True},
            {"name": "height", "required": True}
        ],
        "check_selection": False
    },
    "cd": {
        "description": "Changes the current directory",
        "pattern": r"^cd\s+(.+)$",
        "parameters": [
            {"name": "path", "required": True}
        ],
        "check_selection": False
    },
    "pwd": {
        "description": "Prints the current directory",
        "pattern": r"^pwd$",
        "parameters": [],
        "check_selection": False
    },
    "ls": {
        "description": "Lists files in the current directory",
        "pattern": r"^ls(?:\s+(.+))?$",
        "parameters": [
            {"name": "path", "required": False}
        ],
        "check_selection": False
    },
    "system": {
        "description": "Executes a system command",
        "pattern": r"^system\s+(.+)$",
        "parameters": [
            {"name": "command", "required": True}
        ],
        "check_selection": False
    },
    "help": {
        "description": "Shows help for a command",
        "pattern": r"^help(?:\s+(.+))?$",
        "parameters": [
            {"name": "command", "required": False}
        ],
        "check_selection": False
    }
}

# Common error patterns in PyMOL
ERROR_PATTERNS = {
    "SYNTAX_ERROR": [
        r"Syntax error",
        r"invalid syntax",
        r"Unknown command"
    ],
    "SELECTION_ERROR": [
        r"Invalid selection",
        r"No atoms selected",
        r"Selection not found",
        r"Selection \S+ doesn't exist"
    ],
    "OBJECT_NOT_FOUND": [
        r"object \S+ not found",
        r"Object \S+ does not exist",
        r"Unable to find object named \S+"
    ],
    "ATOM_NOT_FOUND": [
        r"No atoms matched",
        r"No atoms in selection",
        r"Atom not found"
    ],
    "FILE_ERROR": [
        r"Unable to open file",
        r"No such file",
        r"Permission denied",
        r"Error reading file",
        r"Error writing file"
    ],
    "CONNECTION_ERROR": [
        r"Connection refused",
        r"Network error",
        r"Timeout",
        r"Failed to fetch"
    ],
    "PARAMETER_ERROR": [
        r"Incorrect number of parameters",
        r"Invalid parameter",
        r"Parameter out of range"
    ]
}