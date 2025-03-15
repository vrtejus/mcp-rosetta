from mcp.server.fastmcp import FastMCP, Context
import os
import tempfile
import asyncio
import shutil
import glob
import json
import uuid
from typing import Optional, List, Dict, Any, Tuple
from enum import Enum
import subprocess
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("rosetta-design-server")

# Create an MCP server
mcp = FastMCP("RosettaDesignServer")

# Default directories for input, output, and temporary files
BASE_DIR = os.environ.get(
    "ROSETTA_BASE_DIR",
    os.path.expanduser("~/rosetta_design_mcp")
)
INPUT_DIR = os.path.join(BASE_DIR, "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
TEMP_DIR = os.path.join(BASE_DIR, "temp")
DB_DIR = os.path.join(BASE_DIR, "database")

# Create directories if they don't exist
os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)
os.makedirs(DB_DIR, exist_ok=True)

# Path to Rosetta scripts and applications
# These should be set in environment variables or overridden by the user
ROSETTA_BIN_PATH = os.environ.get("ROSETTA_BIN_PATH", "")
ROSETTA_DB_PATH = os.environ.get("ROSETTA_DB_PATH", "")

# Detect operating system and set appropriate executable extensions
import platform
if platform.system() == "Darwin":  # macOS
    ROSETTA_SCRIPTS_PATH = os.path.join(ROSETTA_BIN_PATH, "rosetta_scripts.static.macosclangrelease")
    ANTIBODY_APP_PATH = os.path.join(ROSETTA_BIN_PATH, "antibody.static.macosclangrelease")
    logger.info("Using macOS executables")
else:  # Linux and others
    ROSETTA_SCRIPTS_PATH = os.path.join(ROSETTA_BIN_PATH, "rosetta_scripts.default.linuxgccrelease")
    ANTIBODY_APP_PATH = os.path.join(ROSETTA_BIN_PATH, "antibody.default.linuxgccrelease")
    logger.info("Using Linux executables")

# Project state tracking
class ProjectType(Enum):
    ANTIBODY = "antibody"
    SMALL_MOLECULE = "small_molecule"

# Dictionary to store active projects
active_projects = {}

class Project:
    def __init__(self, project_id: str, project_type: ProjectType, name: str):
        self.id = project_id
        self.type = project_type
        self.name = name
        self.steps_completed = []
        self.current_step = None
        self.input_files = {}
        self.output_files = {}
        self.parameters = {}
        self.working_dir = os.path.join(OUTPUT_DIR, project_id)
        
        # Create project directory
        os.makedirs(self.working_dir, exist_ok=True)
    
    def add_input_file(self, file_type: str, file_path: str):
        self.input_files[file_type] = file_path
    
    def add_output_file(self, file_type: str, file_path: str):
        self.output_files[file_type] = file_path
    
    def add_parameter(self, key: str, value: Any):
        self.parameters[key] = value
    
    def complete_step(self, step_name: str):
        self.steps_completed.append(step_name)
        self.current_step = None
    
    def start_step(self, step_name: str):
        self.current_step = step_name
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert project to dictionary for json serialization"""
        return {
            "id": self.id,
            "type": self.type.value,
            "name": self.name,
            "steps_completed": self.steps_completed,
            "current_step": self.current_step,
            "input_files": self.input_files,
            "output_files": self.output_files,
            "parameters": self.parameters,
            "working_dir": self.working_dir
        }

def save_project_state():
    """Save the current state of all projects to a JSON file"""
    state_file = os.path.join(BASE_DIR, "project_state.json")
    state_data = {
        project_id: project.to_dict() 
        for project_id, project in active_projects.items()
    }
    with open(state_file, "w") as f:
        json.dump(state_data, f, indent=2)

def load_project_state():
    """Load project state from JSON file if it exists"""
    state_file = os.path.join(BASE_DIR, "project_state.json")
    if os.path.exists(state_file):
        with open(state_file, "r") as f:
            state_data = json.load(f)
        
        for project_id, project_data in state_data.items():
            project = Project(
                project_id=project_data["id"],
                project_type=ProjectType(project_data["type"]),
                name=project_data["name"]
            )
            project.steps_completed = project_data["steps_completed"]
            project.current_step = project_data["current_step"]
            project.input_files = project_data["input_files"]
            project.output_files = project_data["output_files"]
            project.parameters = project_data["parameters"]
            project.working_dir = project_data["working_dir"]
            
            active_projects[project_id] = project

# Helper functions
async def run_command(cmd: List[str], cwd: Optional[str] = None, env: Optional[Dict[str, str]] = None) -> Tuple[str, str, int]:
    """Run a command asynchronously and return stdout, stderr, and return code"""
    logger.info(f"Running command: {' '.join(cmd)}")
    
    # Merge environment variables
    process_env = os.environ.copy()
    if env:
        process_env.update(env)
    
    # Add Rosetta database path to environment
    if ROSETTA_DB_PATH:
        process_env["ROSETTA_DATABASE_PATH"] = ROSETTA_DB_PATH
    
    process = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        cwd=cwd,
        env=process_env
    )
    
    stdout, stderr = await process.communicate()
    return stdout.decode(), stderr.decode(), process.returncode

def create_rosetta_options_file(options: Dict[str, Any], output_path: str) -> str:
    # Set of keys that should be separated by a space rather than a colon.
    space_separated_keys = {"database", "fasta", "fasta_heavy", "fasta_light", "in:file:s", "in:file:extra_res_fa", "in:file:native"}
    
    with open(output_path, "w") as f:
        for key, value in options.items():
            if value is None:
                f.write(f"-{key}\n")
            else:
                if key in space_separated_keys:
                    f.write(f"-{key} {value}\n")
                else:
                    f.write(f"-{key}:{value}\n")
    return output_path


# =============================================================================
# Project Management Tools
# =============================================================================

@mcp.tool()
async def list_projects(ctx: Context = None) -> str:
    """List all active protein design projects
    
    Returns:
        A formatted list of all projects with their IDs, names, and status
    """
    if not active_projects:
        return "No active projects found."
    
    result = "Active Projects:\n\n"
    for project_id, project in active_projects.items():
        result += f"ID: {project_id}\n"
        result += f"Name: {project.name}\n"
        result += f"Type: {project.type.value}\n"
        result += f"Steps Completed: {', '.join(project.steps_completed) if project.steps_completed else 'None'}\n"
        result += f"Current Step: {project.current_step or 'None'}\n\n"
    
    return result

@mcp.tool()
async def create_project(
    name: str,
    project_type: str,
    rosetta_bin_path: str = None,
    rosetta_db_path: str = None,
    ctx: Context = None
) -> str:
    """Create a new protein design project
    
    Args:
        name: Descriptive name for the project
        project_type: Type of project (either 'antibody' or 'small_molecule')
        
    Returns:
        Project ID and confirmation message
    """
    # Validate project type
    try:
        project_type_enum = ProjectType(project_type.lower())
    except ValueError:
        return f"Error: Invalid project type '{project_type}'. Use 'antibody' or 'small_molecule'."
    
    # Create project ID
    project_id = str(uuid.uuid4())[:8]
    
    # Set custom Rosetta paths if provided
    global ROSETTA_BIN_PATH, ROSETTA_DB_PATH, ROSETTA_SCRIPTS_PATH, ANTIBODY_APP_PATH
    
    if rosetta_bin_path:
        ROSETTA_BIN_PATH = rosetta_bin_path
        # Update dependent paths
        if platform.system() == "Darwin":  # macOS
            ROSETTA_SCRIPTS_PATH = os.path.join(ROSETTA_BIN_PATH, "rosetta_scripts.static.macosclangrelease")
            ANTIBODY_APP_PATH = os.path.join(ROSETTA_BIN_PATH, "antibody.static.macosclangrelease")
        else:  # Linux and others
            ROSETTA_SCRIPTS_PATH = os.path.join(ROSETTA_BIN_PATH, "rosetta_scripts.default.linuxgccrelease")
            ANTIBODY_APP_PATH = os.path.join(ROSETTA_BIN_PATH, "antibody.default.linuxgccrelease")
    
    if rosetta_db_path:
        ROSETTA_DB_PATH = rosetta_db_path
    
    # Create project
    project = Project(
        project_id=project_id,
        project_type=project_type_enum,
        name=name
    )
    
    # Store Rosetta paths in project parameters
    project.add_parameter("rosetta_bin_path", ROSETTA_BIN_PATH)
    project.add_parameter("rosetta_db_path", ROSETTA_DB_PATH)
    
    # Add to active projects
    active_projects[project_id] = project
    
    # Save project state
    save_project_state()
    
    return f"Project '{name}' created with ID: {project_id}"

@mcp.tool()
async def get_project_details(
    project_id: str,
    ctx: Context = None
) -> str:
    """Get detailed information about a specific project
    
    Args:
        project_id: ID of the project to examine
        
    Returns:
        Detailed project information including files and parameters
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    result = f"Project: {project.name} (ID: {project.id})\n"
    result += f"Type: {project.type.value}\n"
    result += f"Working Directory: {project.working_dir}\n\n"
    
    result += "Steps Completed:\n"
    if project.steps_completed:
        for i, step in enumerate(project.steps_completed, 1):
            result += f"  {i}. {step}\n"
    else:
        result += "  None\n"
    
    if project.current_step:
        result += f"\nCurrent Step: {project.current_step}\n"
    
    result += "\nInput Files:\n"
    if project.input_files:
        for file_type, file_path in project.input_files.items():
            result += f"  {file_type}: {file_path}\n"
    else:
        result += "  None\n"
    
    result += "\nOutput Files:\n"
    if project.output_files:
        for file_type, file_path in project.output_files.items():
            result += f"  {file_type}: {file_path}\n"
    else:
        result += "  None\n"
    
    result += "\nParameters:\n"
    if project.parameters:
        for key, value in project.parameters.items():
            result += f"  {key}: {value}\n"
    else:
        result += "  None\n"
    
    return result

@mcp.tool()
async def upload_target_pdb(
    project_id: str,
    pdb_content: str = None,
    file_name: str = "target.pdb",
    pdb_path: str = None,
    ctx: Context = None
) -> str:
    """Upload a target PDB file for a project
    
    Args:
        project_id: ID of the project
        pdb_content: Content of the PDB file as a string (optional if pdb_path is provided)
        file_name: Name to save the file as (default: target.pdb)
        pdb_path: Path to an existing PDB file (optional if pdb_content is provided)
        
    Returns:
        Confirmation message with the file path
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    # Create the target file path
    target_file_path = os.path.join(project.working_dir, file_name)
    
    # Check if we're using an existing file path or direct content
    if pdb_path:
        # Validate that the path exists
        if not os.path.exists(pdb_path):
            return f"Error: PDB file at path '{pdb_path}' not found."
        
        # Copy the file to the project directory
        try:
            shutil.copy(pdb_path, target_file_path)
            logger.info(f"Copied PDB from {pdb_path} to {target_file_path}")
        except Exception as e:
            return f"Error copying PDB file: {str(e)}"
    elif pdb_content:
        # Write the PDB content to the file
        with open(target_file_path, "w") as f:
            f.write(pdb_content)
        logger.info(f"Wrote PDB content to {target_file_path}")
    else:
        return "Error: Either pdb_content or pdb_path must be provided."
    
    # Update project with input file
    project.add_input_file("target_pdb", target_file_path)
    
    # Save project state
    save_project_state()
    
    return f"Target PDB saved to {target_file_path}"

# =============================================================================
# Antibody Modeling Tools
# =============================================================================

@mcp.tool()
async def antibody_setup(
    project_id: str,
    heavy_chain: str = None,
    light_chain: str = None,
    heavy_chain_file: str = None,
    light_chain_file: str = None,
    ctx: Context = None
) -> str:
    """Set up an antibody modeling project with heavy and light chain sequences
    
    Args:
        project_id: ID of the project
        heavy_chain: Amino acid sequence of the heavy chain (optional if heavy_chain_file is provided)
        light_chain: Amino acid sequence of the light chain (optional if light_chain_file is provided)
        heavy_chain_file: Path to a FASTA file containing the heavy chain sequence (optional)
        light_chain_file: Path to a FASTA file containing the light chain sequence (optional)
        
    Returns:
        Confirmation message
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.ANTIBODY:
        return f"Error: Project '{project_id}' is not an antibody project."
    
    # Handle heavy chain
    heavy_file_path = os.path.join(project.working_dir, "heavy.fasta")
    
    if heavy_chain_file:
        # Validate that the path exists
        if not os.path.exists(heavy_chain_file):
            return f"Error: Heavy chain file at path '{heavy_chain_file}' not found."
        
        try:
            shutil.copy(heavy_chain_file, heavy_file_path)
            logger.info(f"Copied heavy chain file from {heavy_chain_file} to {heavy_file_path}")
            
            # Extract sequence length for reporting
            with open(heavy_file_path, 'r') as f:
                content = f.read()
                seq_lines = [line.strip() for line in content.split('\n') if line.strip() and not line.startswith('>')]
                heavy_chain_length = sum(len(line) for line in seq_lines)
        except Exception as e:
            return f"Error copying heavy chain file: {str(e)}"
    elif heavy_chain:
        # Write the sequence to a FASTA file
        with open(heavy_file_path, "w") as f:
            f.write(f">heavy_chain\n{heavy_chain}\n")
        heavy_chain_length = len(heavy_chain)
    else:
        return "Error: Either heavy_chain or heavy_chain_file must be provided."
    
    # Handle light chain
    light_file_path = os.path.join(project.working_dir, "light.fasta")
    
    if light_chain_file:
        # Validate that the path exists
        if not os.path.exists(light_chain_file):
            return f"Error: Light chain file at path '{light_chain_file}' not found."
        
        try:
            shutil.copy(light_chain_file, light_file_path)
            logger.info(f"Copied light chain file from {light_chain_file} to {light_file_path}")
            
            # Extract sequence length for reporting
            with open(light_file_path, 'r') as f:
                content = f.read()
                seq_lines = [line.strip() for line in content.split('\n') if line.strip() and not line.startswith('>')]
                light_chain_length = sum(len(line) for line in seq_lines)
        except Exception as e:
            return f"Error copying light chain file: {str(e)}"
    elif light_chain:
        # Write the sequence to a FASTA file
        with open(light_file_path, "w") as f:
            f.write(f">light_chain\n{light_chain}\n")
        light_chain_length = len(light_chain)
    else:
        return "Error: Either light_chain or light_chain_file must be provided."
    
    # Update project with input files
    project.add_input_file("heavy_chain", heavy_file_path)
    project.add_input_file("light_chain", light_file_path)
    
    # Save project state
    save_project_state()
    
    return f"Antibody chains set up successfully:\nHeavy chain: {heavy_chain_length} residues\nLight chain: {light_chain_length} residues"

@mcp.tool()
async def run_antibody_modeling(
    project_id: str,
    num_models: int = 10,
    use_homology_models: bool = True,
    refine_models: bool = True,
    custom_antibody_path: str = None,
    ctx: Context = None
) -> str:
    """Run the RosettaAntibody modeling pipeline
    
    Args:
        project_id: ID of the project
        num_models: Number of models to generate (default: 10)
        use_homology_models: Whether to use homology models for CDR loops (default: True)
        refine_models: Whether to refine the final models (default: True)
        
    Returns:
        Summary of the modeling run
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.ANTIBODY:
        return f"Error: Project '{project_id}' is not an antibody project."
    
    # Check if required input files exist
    if "heavy_chain" not in project.input_files or "light_chain" not in project.input_files:
        return "Error: Heavy and light chain sequences are required. Run antibody_setup first."
    
    # Update project parameters
    project.add_parameter("num_models", num_models)
    project.add_parameter("use_homology_models", use_homology_models)
    project.add_parameter("refine_models", refine_models)
    
    # Start modeling step
    project.start_step("antibody_modeling")
    
    # Create options file
    options = {
        "database": ROSETTA_DB_PATH,
        "fasta": None,  # Using flags without values
        "fasta_heavy": project.input_files["heavy_chain"],
        "fasta_light": project.input_files["light_chain"],
        "out:path:all": project.working_dir,
        "out:prefix": "model_",
        "nstruct": num_models,
        "antibody:auto_generate_kink_constraint": None,
        "antibody:auto_generate_length_constraint": None
    }
    
    if not use_homology_models:
        options["antibody:do_homology_modeling"] = "false"
    
    if refine_models:
        options["antibody:refine"] = None
    
    options_file = os.path.join(project.working_dir, "antibody_options.txt")
    create_rosetta_options_file(options, options_file)
    
    # Set antibody executable path
    antibody_path = "/Users/ammachi/Downloads/rosetta-binary/main/source/bin/antibody.static.macosclangrelease"

    logger.info(f"Using antibody executable at: {antibody_path}")
    if ctx:
        ctx.info(f"Using antibody executable at: {antibody_path}")
    
    # Verify the antibody executable exists
    if not os.path.exists(antibody_path):
        available_executables = glob.glob(os.path.join(ROSETTA_BIN_PATH, "antibody*.macosclangrelease")) + \
                              glob.glob(os.path.join(ROSETTA_BIN_PATH, "antibody*linuxgccrelease"))
        
        error_msg = f"Error: Antibody executable not found at {antibody_path}.\n"
        
        if available_executables:
            error_msg += "Available antibody executables:\n"
            for exe in available_executables:
                error_msg += f"- {os.path.basename(exe)}\n"
            error_msg += "\nPlease specify one of these using the custom_antibody_path parameter."
        else:
            error_msg += f"No antibody executables found in {ROSETTA_BIN_PATH}."
        
        return error_msg
    
    # Execute RosettaAntibody command
    cmd = [
        antibody_path,
        "@" + options_file
    ]
    
    if ctx:
        ctx.info(f"Running antibody modeling with command: {' '.join(cmd)}")
    
    try:
        stdout, stderr, return_code = await run_command(cmd, cwd=project.working_dir)
        
        # Log output
        log_file = os.path.join(project.working_dir, "antibody_modeling.log")
        with open(log_file, "w") as f:
            f.write("STDOUT:\n")
            f.write(stdout)
            f.write("\nSTDERR:\n")
            f.write(stderr)
        
        if return_code != 0:
            project.current_step = None
            return f"Error: Antibody modeling failed with return code {return_code}. See log file: {log_file}"
        
        # Find output PDB files
        pdb_files = glob.glob(os.path.join(project.working_dir, "model_*.pdb"))
        
        # Store output file paths
        for i, pdb_file in enumerate(pdb_files):
            project.add_output_file(f"model_{i+1}", pdb_file)
        
        # Complete step
        project.complete_step("antibody_modeling")
        
        # Save project state
        save_project_state()
        
        return f"Antibody modeling completed successfully. Generated {len(pdb_files)} models."
    
    except Exception as e:
        project.current_step = None
        error_msg = f"Error running antibody modeling: {str(e)}"
        logger.error(error_msg)
        return error_msg

@mcp.tool()
async def analyze_antibody_models(
    project_id: str,
    target_pdb_path: Optional[str] = None,
    ctx: Context = None
) -> str:
    """Analyze antibody models to rank them and identify the best ones
    
    Args:
        project_id: ID of the project
        target_pdb_path: Optional path to a target protein PDB for interface analysis
        
    Returns:
        Analysis results with rankings and scores
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.ANTIBODY:
        return f"Error: Project '{project_id}' is not an antibody project."
    
    # Check if modeling was completed
    if "antibody_modeling" not in project.steps_completed:
        return "Error: Antibody modeling step has not been completed. Run run_antibody_modeling first."
    
    # Find model files
    model_files = []
    for key, path in project.output_files.items():
        if key.startswith("model_") and path.endswith(".pdb"):
            model_files.append(path)
    
    if not model_files:
        return "Error: No antibody models found."
    
    # Start analysis step
    project.start_step("antibody_analysis")
    
    # Create a simple script to score models
    score_file = os.path.join(project.working_dir, "score.xml")
    
    with open(score_file, "w") as f:
        f.write("""<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref2015" weights="ref2015"/>
    </SCOREFXNS>
    <PROTOCOLS>
        <Add mover="scoremover"/>
    </PROTOCOLS>
    <MOVERS>
        <ScoreInterface name="scoremover" scorefxn="ref2015"/>
    </MOVERS>
</ROSETTASCRIPTS>""")
    
    # Analyze each model
    results = []
    
    for model_file in model_files:
        model_name = os.path.basename(model_file)
        
        # Create options file
        options = {
            "database": ROSETTA_DB_PATH,
            "in:file:s": model_file,
            "parser:protocol": score_file,
            "out:file:score_only": f"{model_name}.score",
            "out:file:silent": os.path.join(project.working_dir, f"{model_name}.silent"),
            "out:file:scorefile": os.path.join(project.working_dir, "scores.sc")
        }
        
        # If target PDB is provided, add interface analysis
        if target_pdb_path:
            options["in:file:native"] = target_pdb_path
            options["out:file:silent"] = os.path.join(project.working_dir, f"{model_name}_interface.silent")
        
        options_file = os.path.join(project.working_dir, f"{model_name}_options.txt")
        create_rosetta_options_file(options, options_file)
        
        # Execute Rosetta Scripts command for scoring
        cmd = [
            ROSETTA_SCRIPTS_PATH,
            "@" + options_file
        ]
        
        if ctx:
            ctx.info(f"Running analysis on {model_name} with command: {' '.join(cmd)}")
        
        try:
            stdout, stderr, return_code = await run_command(cmd, cwd=project.working_dir)
            
            # Parse scores from stdout
            score_lines = [line for line in stdout.split("\n") if "total_score" in line or model_name in line]
            score_data = {}
            
            for line in score_lines:
                if ":" in line:
                    key, value = line.split(":", 1)
                    score_data[key.strip()] = value.strip()
            
            results.append({
                "model": model_name,
                "score": score_data.get("total_score", "N/A"),
                "rmsd": score_data.get("rmsd", "N/A") if target_pdb_path else "N/A"
            })
            
        except Exception as e:
            return f"Error analyzing model {model_name}: {str(e)}"
    
    # Sort results by score (lower is better)
    sorted_results = sorted(results, key=lambda x: float(x["score"]) if x["score"] != "N/A" else float('inf'))
    
    # Generate report
    analysis_file = os.path.join(project.working_dir, "antibody_analysis.txt")
    
    with open(analysis_file, "w") as f:
        f.write("Antibody Model Analysis\n")
        f.write("======================\n\n")
        f.write("Ranked by Rosetta Score (lower is better):\n\n")
        
        for i, result in enumerate(sorted_results, 1):
            f.write(f"Rank {i}: {result['model']}\n")
            f.write(f"  Score: {result['score']}\n")
            if target_pdb_path:
                f.write(f"  RMSD to target: {result['rmsd']}\n")
            f.write("\n")
    
    # Save best model
    best_model = sorted_results[0]["model"]
    best_model_path = [path for path in model_files if os.path.basename(path) == best_model][0]
    best_model_output = os.path.join(project.working_dir, "best_model.pdb")
    shutil.copy(best_model_path, best_model_output)
    
    project.add_output_file("best_model", best_model_output)
    project.add_output_file("analysis", analysis_file)
    
    # Complete step
    project.complete_step("antibody_analysis")
    
    # Save project state
    save_project_state()
    
    # Format results for return
    result_str = "Antibody Model Analysis Results:\n\n"
    result_str += f"Total models analyzed: {len(results)}\n"
    result_str += f"Best model: {best_model} (Score: {sorted_results[0]['score']})\n\n"
    
    result_str += "Top 5 Models:\n"
    for i, result in enumerate(sorted_results[:5], 1):
        result_str += f"{i}. {result['model']} - Score: {result['score']}"
        if target_pdb_path and result['rmsd'] != "N/A":
            result_str += f", RMSD: {result['rmsd']}"
        result_str += "\n"
    
    result_str += f"\nDetailed analysis saved to: {analysis_file}"
    result_str += f"\nBest model saved to: {best_model_output}"
    
    return result_str

# =============================================================================
# Small Molecule Design Tools
# =============================================================================

@mcp.tool()
async def small_molecule_setup(
    project_id: str,
    ligand_sdf: str = None,
    binding_site_residues: str = None,
    ligand_sdf_path: str = None,
    ctx: Context = None
) -> str:
    """Set up a small molecule interface design project
    
    Args:
        project_id: ID of the project
        ligand_sdf: SDF content of the ligand (optional if ligand_sdf_path is provided)
        binding_site_residues: Comma-separated list of binding site residue IDs (optional)
        ligand_sdf_path: Path to an existing ligand SDF file (optional if ligand_sdf is provided)
        
    Returns:
        Confirmation message
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.SMALL_MOLECULE:
        return f"Error: Project '{project_id}' is not a small molecule design project."
    
    # Check if target PDB has been uploaded
    if "target_pdb" not in project.input_files:
        return "Error: Target PDB file is required. Run upload_target_pdb first."
    
    # Handle ligand SDF
    if ligand_sdf_path:
        # Validate that the path exists
        if not os.path.exists(ligand_sdf_path):
            return f"Error: Ligand SDF file at path '{ligand_sdf_path}' not found."
        
        # Copy the file to the project directory
        ligand_file_path = os.path.join(project.working_dir, "ligand.sdf")
        try:
            shutil.copy(ligand_sdf_path, ligand_file_path)
            logger.info(f"Copied ligand SDF from {ligand_sdf_path} to {ligand_file_path}")
            project.add_input_file("ligand_sdf", ligand_file_path)
        except Exception as e:
            return f"Error copying ligand SDF file: {str(e)}"
    elif ligand_sdf:
        # Write the SDF content to the file
        ligand_file_path = os.path.join(project.working_dir, "ligand.sdf")
        with open(ligand_file_path, "w") as f:
            f.write(ligand_sdf)
        project.add_input_file("ligand_sdf", ligand_file_path)
    
    # If binding site residues are provided, save them
    if binding_site_residues:
        binding_site_file = os.path.join(project.working_dir, "binding_site.txt")
        with open(binding_site_file, "w") as f:
            f.write(binding_site_residues)
        project.add_input_file("binding_site", binding_site_file)
        project.add_parameter("binding_site_residues", binding_site_residues)
    
    # Save project state
    save_project_state()
    
    result = "Small molecule design project set up successfully.\n"
    if ligand_sdf or ligand_sdf_path:
        result += "Ligand SDF file saved.\n"
    if binding_site_residues:
        result += f"Binding site residues saved: {binding_site_residues}\n"
    
    return result

@mcp.tool()
async def prepare_ligand_params(
    project_id: str,
    ligand_name: str = "LIG",
    ctx: Context = None
) -> str:
    """Generate Rosetta parameter files for a ligand
    
    Args:
        project_id: ID of the project
        ligand_name: Three-letter code for the ligand (default: LIG)
        
    Returns:
        Status of parameter file generation
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.SMALL_MOLECULE:
        return f"Error: Project '{project_id}' is not a small molecule design project."
    
    # Check if ligand SDF has been uploaded
    if "ligand_sdf" not in project.input_files:
        return "Error: Ligand SDF file is required. Run small_molecule_setup with a ligand_sdf first."
    
    # Start step
    project.start_step("ligand_params")
    
    # Paths
    ligand_sdf_path = project.input_files["ligand_sdf"]
    output_params_path = os.path.join(project.working_dir, f"{ligand_name}.params")
    
    # Create script for generating params
    script_path = os.path.join(project.working_dir, "generate_params.py")
    
    with open(script_path, "w") as f:
        f.write("""
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem

# Input SDF file
sdf_file = sys.argv[1]
ligand_name = sys.argv[2]
output_params = sys.argv[3]

# Read the molecule
mol = Chem.SDMolSupplier(sdf_file)[0]
if mol is None:
    print("Error: Could not read molecule from SDF file")
    sys.exit(1)

# Add hydrogen atoms
mol = Chem.AddHs(mol)

# Generate 3D coordinates if not present
if not mol.GetNumConformers():
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

# Write the molecule to PDB format
pdb_file = os.path.splitext(sdf_file)[0] + ".pdb"
Chem.MolToPDBFile(mol, pdb_file)

# Run molfile_to_params.py script
os.system(f"python {os.environ.get('ROSETTA_SCRIPTS_PATH', '')}/python/apps/public/molfile_to_params.py -n {ligand_name} -p {output_params} {pdb_file}")

print(f"Params file generated: {output_params}")
""")
    
    # Execute script
    cmd = [
        "python",
        script_path,
        ligand_sdf_path,
        ligand_name,
        output_params_path
    ]
    
    if ctx:
        ctx.info(f"Generating ligand parameters with command: {' '.join(cmd)}")
    
    try:
        stdout, stderr, return_code = await run_command(cmd, cwd=project.working_dir)
        
        if return_code != 0:
            project.current_step = None
            return f"Error: Ligand parameter generation failed with return code {return_code}.\nError: {stderr}"
        
        # If successful, add params file to project
        if os.path.exists(output_params_path):
            project.add_output_file("ligand_params", output_params_path)
            project.add_parameter("ligand_name", ligand_name)
            
            # Complete step
            project.complete_step("ligand_params")
            
            # Save project state
            save_project_state()
            
            return f"Ligand parameters generated successfully: {output_params_path}"
        else:
            project.current_step = None
            return "Error: Parameter file was not created."
    
    except Exception as e:
        project.current_step = None
        error_msg = f"Error generating ligand parameters: {str(e)}"
        logger.error(error_msg)
        return error_msg

@mcp.tool()
async def run_interface_design(
    project_id: str,
    num_designs: int = 10,
    design_shell: float = 8.0,
    design_cycles: int = 3,
    ctx: Context = None
) -> str:
    """Run interface design to optimize protein-ligand binding
    
    Args:
        project_id: ID of the project
        num_designs: Number of design models to generate (default: 10)
        design_shell: Distance from ligand to define design shell in Angstroms (default: 8.0)
        design_cycles: Number of design cycles to run (default: 3)
        
    Returns:
        Status of interface design run
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.SMALL_MOLECULE:
        return f"Error: Project '{project_id}' is not a small molecule design project."
    
    # Check required files
    if "target_pdb" not in project.input_files:
        return "Error: Target PDB file is required. Run upload_target_pdb first."
    
    if "ligand_params" not in project.output_files:
        return "Error: Ligand parameter file is required. Run prepare_ligand_params first."
    
    # Start step
    project.start_step("interface_design")
    
    # Save parameters
    project.add_parameter("num_designs", num_designs)
    project.add_parameter("design_shell", design_shell)
    project.add_parameter("design_cycles", design_cycles)
    
    # Create RosettaScripts XML for interface design
    ligand_name = project.parameters.get("ligand_name", "LIG")
    binding_site_residues = project.parameters.get("binding_site_residues", "")
    
    # Create binding site selection if provided
    binding_site_selector = ""
    if binding_site_residues:
        binding_site_selector = f"""<ResidueSelector name="binding_site" residue_numbers="{binding_site_residues}"/>"""
    else:
        binding_site_selector = f"""<NeighborhoodResidueSelector name="binding_site" focus_selector="ligand" distance="{design_shell}"/>"""
    
    xml_script = f"""<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref2015" weights="ref2015"/>
        <ScoreFunction name="ref2015_soft" weights="ref2015">
            <Reweight scoretype="fa_rep" weight="0.5"/>
        </ScoreFunction>
    </SCOREFXNS>
    
    <RESIDUE_SELECTORS>
        <ResidueSelector name="ligand" chains="X"/>
        {binding_site_selector}
        <NotResidueSelector name="not_ligand" selector="ligand"/>
        <AndResidueSelector name="design_shell" selectors="binding_site,not_ligand"/>
        <NeighborhoodResidueSelector name="near_shell" focus_selector="design_shell" distance="4.0"/>
        <OrResidueSelector name="all_shell" selectors="design_shell,near_shell"/>
    </RESIDUE_SELECTORS>
    
    <TASKOPERATIONS>
        <InitializeFromCommandline name="init"/>
        <DesignRestrictions name="design_shell">
            <Action selector="design_shell" operation="ALLAA" />
            <Action selector="near_shell" operation="NATRO" />
            <Action selector="ligand" operation="NATRO" />
        </DesignRestrictions>
        <RestrictToResidueProperties name="polar_only">
            <Property property="POLAR" />
            <Action selector="design_shell" operation="RESTRICT_TO" />
        </RestrictToResidueProperties>
        <RestrictToRepacking name="restrict_repack" />
    </TASKOPERATIONS>
    
    <MOVERS>
        <AddMissingAtoms name="add_atoms"/>
        <MinMover name="min" scorefxn="ref2015" chi="1" bb="0" type="lbfgs_armijo_nonmonotone" tolerance="0.001"/>
        
        <PackRotamersMover name="pack" scorefxn="ref2015_soft" task_operations="init,design_shell"/>
        
        <ParsedProtocol name="design_cycle">
            <Add mover="pack"/>
            <Add mover="min"/>
        </ParsedProtocol>
        
        <ParsedProtocol name="design">
            <Add mover="add_atoms"/>"""
    
    # Add design cycles
    for i in range(design_cycles):
        xml_script += f"""
            <Add mover="design_cycle"/>"""
    
    xml_script += """
        </ParsedProtocol>
    </MOVERS>
    
    <PROTOCOLS>
        <Add mover="design"/>
    </PROTOCOLS>
    
    <OUTPUT scorefxn="ref2015"/>
</ROSETTASCRIPTS>"""
    
    # Write XML script to file
    script_path = os.path.join(project.working_dir, "interface_design.xml")
    with open(script_path, "w") as f:
        f.write(xml_script)
    
    # Create options file
    target_pdb = project.input_files["target_pdb"]
    ligand_params = project.output_files["ligand_params"]
    
    options = {
        "database": ROSETTA_DB_PATH,
        "in:file:s": target_pdb,
        "in:file:extra_res_fa": ligand_params,
        "parser:protocol": script_path,
        "out:file:silent": os.path.join(project.working_dir, "designs.silent"),
        "out:suffix": "_design",
        "nstruct": num_designs,
        "out:file:scorefile": "scores.sc",
        "ignore_unrecognized_res": None,
        "ignore_zero_occupancy": None,
        "ex1": None,
        "ex2": None
    }
    
    options_file = os.path.join(project.working_dir, "design_options.txt")
    create_rosetta_options_file(options, options_file)
    
    # Execute RosettaScripts command
    cmd = [
        ROSETTA_SCRIPTS_PATH,
        "@" + options_file
    ]
    
    if ctx:
        ctx.info(f"Running interface design with command: {' '.join(cmd)}")
    
    try:
        stdout, stderr, return_code = await run_command(cmd, cwd=project.working_dir)
        
        # Log output
        log_file = os.path.join(project.working_dir, "interface_design.log")
        with open(log_file, "w") as f:
            f.write("STDOUT:\n")
            f.write(stdout)
            f.write("\nSTDERR:\n")
            f.write(stderr)
        
        if return_code != 0:
            project.current_step = None
            return f"Error: Interface design failed with return code {return_code}. See log file: {log_file}"
        
        # Extract PDB models from silent file
        silent_file = os.path.join(project.working_dir, "designs.silent")
        
        if not os.path.exists(silent_file):
            project.current_step = None
            return f"Error: Silent file not generated. See log file: {log_file}"
        
        # Create script to extract PDBs from silent file
        extract_script = os.path.join(project.working_dir, "extract_pdbs.py")
        with open(extract_script, "w") as f:
            f.write("""
import os
import sys
import subprocess

silent_file = sys.argv[1]
output_dir = sys.argv[2]

# Create output directory
os.makedirs(output_dir, exist_ok=True)

# Run extract_pdbs
cmd = [
    os.path.join(os.environ.get("ROSETTA_BIN_PATH", ""), "extract_pdbs.default.linuxgccrelease"),
    "-in:file:silent", silent_file,
    "-out:path:pdb", output_dir
]

subprocess.run(cmd, check=True)
print(f"Extracted PDBs to {output_dir}")
""")
        
        # Extract PDBs
        pdb_output_dir = os.path.join(project.working_dir, "pdb_designs")
        extract_cmd = [
            "python",
            extract_script,
            silent_file,
            pdb_output_dir
        ]
        
        if ctx:
            ctx.info(f"Extracting PDBs with command: {' '.join(extract_cmd)}")
        
        extract_stdout, extract_stderr, extract_return_code = await run_command(extract_cmd, cwd=project.working_dir)
        
        if extract_return_code != 0:
            project.current_step = None
            return f"Error: Failed to extract PDBs from silent file. Error: {extract_stderr}"
        
        # Find extracted PDB files
        pdb_files = glob.glob(os.path.join(pdb_output_dir, "*.pdb"))
        
        # Store output file paths
        for i, pdb_file in enumerate(pdb_files):
            base_name = os.path.basename(pdb_file)
            project.add_output_file(f"design_{i+1}", pdb_file)
        
        project.add_output_file("designs_silent", silent_file)
        
        # Complete step
        project.complete_step("interface_design")
        
        # Save project state
        save_project_state()
        
        return f"Interface design completed successfully. Generated {len(pdb_files)} design models.\nOutput directory: {pdb_output_dir}"
    
    except Exception as e:
        project.current_step = None
        error_msg = f"Error running interface design: {str(e)}"
        logger.error(error_msg)
        return error_msg

@mcp.tool()
async def analyze_interface_designs(
    project_id: str,
    ctx: Context = None
) -> str:
    """Analyze interface design models to rank and score them
    
    Args:
        project_id: ID of the project
        
    Returns:
        Analysis results with rankings and scores
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    if project.type != ProjectType.SMALL_MOLECULE:
        return f"Error: Project '{project_id}' is not a small molecule design project."
    
    # Check if design was completed
    if "interface_design" not in project.steps_completed:
        return "Error: Interface design step has not been completed. Run run_interface_design first."
    
    # Start step
    project.start_step("design_analysis")
    
    # Find design files
    design_files = []
    for key, path in project.output_files.items():
        if key.startswith("design_") and path.endswith(".pdb"):
            design_files.append(path)
    
    if not design_files:
        project.current_step = None
        return "Error: No design models found."
    
    # Create analysis script
    analysis_script = os.path.join(project.working_dir, "analyze_interfaces.xml")
    
    with open(analysis_script, "w") as f:
        f.write("""<ROSETTASCRIPTS>
    <SCOREFXNS>
        <ScoreFunction name="ref2015" weights="ref2015"/>
    </SCOREFXNS>
    
    <RESIDUE_SELECTORS>
        <ResidueSelector name="ligand" chains="X"/>
        <NotResidueSelector name="not_ligand" selector="ligand"/>
    </RESIDUE_SELECTORS>
    
    <FILTERS>
        <InterfaceAnalyzer name="interface" scorefxn="ref2015" interface_selector="ligand" packstat="true" pack_input="true" pack_separated="true" compute_packstat="true" ligand_mode="true"/>
    </FILTERS>
    
    <MOVERS>
        <RunSimpleFilters name="analyze" filters="interface"/>
    </MOVERS>
    
    <PROTOCOLS>
        <Add mover="analyze"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>""")
    
    # Create options file template
    options_template = {
        "database": ROSETTA_DB_PATH,
        "parser:protocol": analysis_script,
        "in:file:extra_res_fa": project.output_files["ligand_params"],
        "out:file:score_only": None,
        "ignore_unrecognized_res": None,
        "ignore_zero_occupancy": None
    }
    
    # Analyze each design
    results = []
    
    for design_file in design_files:
        design_name = os.path.basename(design_file)
        
        # Create specific options file for this design
        options = options_template.copy()
        options["in:file:s"] = design_file
        options["out:file:score_only"] = f"{design_name}.score"
        
        options_file = os.path.join(project.working_dir, f"{design_name}_analysis.txt")
        create_rosetta_options_file(options, options_file)
        
        # Execute Rosetta Scripts command for analysis
        cmd = [
            ROSETTA_SCRIPTS_PATH,
            "@" + options_file
        ]
        
        if ctx:
            ctx.info(f"Analyzing design {design_name} with command: {' '.join(cmd)}")
        
        try:
            stdout, stderr, return_code = await run_command(cmd, cwd=project.working_dir)
            
            # Parse analysis results from stdout
            dG_bind = "N/A"
            packstat = "N/A"
            total_score = "N/A"
            
            for line in stdout.split("\n"):
                if "dG_bind" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "dG_bind":
                            dG_bind = parts[i+1]
                elif "packstat" in line:
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part == "packstat":
                            packstat = parts[i+1]
                elif "total_score" in line and design_name in line:
                    parts = line.split()
                    total_score = parts[1]  # Assuming total_score is the second column
            
            results.append({
                "design": design_name,
                "total_score": total_score,
                "dG_bind": dG_bind,
                "packstat": packstat,
                "path": design_file
            })
            
        except Exception as e:
            return f"Error analyzing design {design_name}: {str(e)}"
    
    # Sort results by binding energy (lower/more negative is better)
    try:
        sorted_results = sorted(results, 
                               key=lambda x: float(x["dG_bind"]) if x["dG_bind"] != "N/A" else float('inf'))
    except Exception:
        # Fall back to total score if dG_bind can't be parsed
        sorted_results = sorted(results, 
                               key=lambda x: float(x["total_score"]) if x["total_score"] != "N/A" else float('inf'))
    
    # Generate report
    analysis_file = os.path.join(project.working_dir, "interface_analysis.txt")
    
    with open(analysis_file, "w") as f:
        f.write("Small Molecule Interface Design Analysis\n")
        f.write("======================================\n\n")
        f.write("Ranked by Binding Energy (dG_bind, lower/more negative is better):\n\n")
        
        for i, result in enumerate(sorted_results, 1):
            f.write(f"Rank {i}: {result['design']}\n")
            f.write(f"  Total Score: {result['total_score']}\n")
            f.write(f"  Binding Energy (dG_bind): {result['dG_bind']}\n")
            f.write(f"  Packing Statistic: {result['packstat']}\n")
            f.write(f"  Path: {result['path']}\n")
            f.write("\n")
    
    # Save best design
    best_design = sorted_results[0]["design"]
    best_design_path = sorted_results[0]["path"]
    best_design_output = os.path.join(project.working_dir, "best_design.pdb")
    shutil.copy(best_design_path, best_design_output)
    
    project.add_output_file("best_design", best_design_output)
    project.add_output_file("analysis", analysis_file)
    
    # Complete step
    project.complete_step("design_analysis")
    
    # Save project state
    save_project_state()
    
    # Format results for return
    result_str = "Small Molecule Interface Design Analysis:\n\n"
    result_str += f"Total designs analyzed: {len(results)}\n"
    result_str += f"Best design: {best_design}\n"
    result_str += f"  Binding Energy (dG_bind): {sorted_results[0]['dG_bind']}\n"
    result_str += f"  Packing Statistic: {sorted_results[0]['packstat']}\n\n"
    
    result_str += "Top 5 Designs:\n"
    for i, result in enumerate(sorted_results[:5], 1):
        result_str += f"{i}. {result['design']} - dG_bind: {result['dG_bind']}, packstat: {result['packstat']}\n"
    
    result_str += f"\nDetailed analysis saved to: {analysis_file}"
    result_str += f"\nBest design saved to: {best_design_output}"
    
    return result_str

# =============================================================================
# Common Tools for Both Design Types
# =============================================================================

@mcp.tool()
async def export_pdb(
    project_id: str,
    output_name: str,
    model_type: str = "best",
    ctx: Context = None
) -> str:
    """Export a PDB file from a project
    
    Args:
        project_id: ID of the project
        output_name: Name to save the exported file as
        model_type: Type of model to export ('best', 'all', or specific model ID)
        
    Returns:
        Path to the exported PDB file
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    # Determine which file(s) to export
    export_files = []
    
    if model_type == "best":
        if "best_model" in project.output_files:
            export_files.append((project.output_files["best_model"], output_name))
        elif "best_design" in project.output_files:
            export_files.append((project.output_files["best_design"], output_name))
        else:
            return "Error: No best model or design found in project."
    
    elif model_type == "all":
        # Find all model files
        model_files = []
        for key, path in project.output_files.items():
            if (key.startswith("model_") or key.startswith("design_")) and path.endswith(".pdb"):
                model_files.append((path, f"{output_name}_{key}"))
        
        if not model_files:
            return "Error: No model or design files found in project."
        
        export_files = model_files
    
    else:
        # Look for specific model
        found = False
        for key, path in project.output_files.items():
            if key == model_type and path.endswith(".pdb"):
                export_files.append((path, output_name))
                found = True
                break
        
        if not found:
            return f"Error: Model '{model_type}' not found in project."
    
    # Create export directory
    export_dir = os.path.join(OUTPUT_DIR, project_id)
    os.makedirs(export_dir, exist_ok=True)
    
    # Export files
    exported_paths = []
    
    for source_path, file_name in export_files:
        if not file_name.endswith(".pdb"):
            file_name += ".pdb"
        
        export_path = os.path.join(export_dir, file_name)
        shutil.copy(source_path, export_path)
        exported_paths.append(export_path)
    
    # Format result
    if len(exported_paths) == 1:
        return f"PDB exported to: {exported_paths[0]}"
    else:
        result = f"Exported {len(exported_paths)} PDB files:\n"
        for path in exported_paths:
            result += f"- {path}\n"
        return result

@mcp.tool()
async def list_rosetta_executables(ctx: Context = None) -> str:
    """List all available Rosetta executable files in the ROSETTA_BIN_PATH directory
    
    Returns:
        List of available Rosetta executables
    """
    if not ROSETTA_BIN_PATH or not os.path.exists(ROSETTA_BIN_PATH):
        return f"Error: ROSETTA_BIN_PATH ({ROSETTA_BIN_PATH}) does not exist or is not set."
    
    # Find all executable files in the bin directory
    executables = []
    for root, dirs, files in os.walk(ROSETTA_BIN_PATH):
        for file in files:
            if file.endswith("macosclangrelease") or file.endswith("linuxgccrelease"):
                executables.append(os.path.join(root, file))
    
    if not executables:
        return "No Rosetta executables found in ROSETTA_BIN_PATH."
    
    # Format result
    result = f"Available Rosetta executables in {ROSETTA_BIN_PATH}:\n\n"
    
    # Group by type
    by_type = {}
    for exe_path in executables:
        exe_name = os.path.basename(exe_path)
        prefix = exe_name.split('.')[0]
        if prefix not in by_type:
            by_type[prefix] = []
        by_type[prefix].append(exe_name)
    
    for prefix, exes in sorted(by_type.items()):
        result += f"{prefix}:\n"
        for exe in sorted(exes):
            result += f"  - {exe}\n"
        result += "\n"
    
    result += "\nTo use a specific executable, provide its full path when needed."
    
    return result

@mcp.tool()
async def list_available_pdbs(
    project_id: str,
    ctx: Context = None
) -> str:
    """List all available PDB files in a project
    
    Args:
        project_id: ID of the project
        
    Returns:
        List of available PDB files with descriptions
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    # Find all PDB files
    pdb_files = []
    
    # Check input files
    for key, path in project.input_files.items():
        if path.endswith(".pdb"):
            pdb_files.append((key, path, "Input"))
    
    # Check output files
    for key, path in project.output_files.items():
        if path.endswith(".pdb"):
            pdb_files.append((key, path, "Output"))
    
    if not pdb_files:
        return "No PDB files found in project."
    
    # Format result
    result = f"Available PDB files for project {project_id} ({project.name}):\n\n"
    
    for file_id, file_path, file_type in pdb_files:
        file_size = os.path.getsize(file_path) / 1024  # KB
        result += f"ID: {file_id}\n"
        result += f"Path: {file_path}\n"
        result += f"Type: {file_type}\n"
        result += f"Size: {file_size:.2f} KB\n\n"
    
    return result

@mcp.tool()
async def visualize_pdb(
    project_id: str,
    model_id: str,
    ctx: Context = None
) -> str:
    """Generate visualization commands for a PDB file
    
    Args:
        project_id: ID of the project
        model_id: ID of the model to visualize
        
    Returns:
        PyMOL or Chimera visualization commands
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    # Find the specified model
    model_path = None
    
    # Check input files
    if model_id in project.input_files:
        model_path = project.input_files[model_id]
    
    # Check output files
    if model_id in project.output_files:
        model_path = project.output_files[model_id]
    
    if not model_path or not model_path.endswith(".pdb"):
        return f"Error: PDB file with ID '{model_id}' not found in project."
    
    # Generate visualization commands
    result = f"Visualization commands for {model_id} ({model_path}):\n\n"
    
    # PyMOL commands
    result += "PyMOL Commands:\n"
    result += f"load {model_path}\n"
    
    if project.type == ProjectType.ANTIBODY:
        result += "show cartoon\n"
        result += "color cyan, chain H\n"
        result += "color pink, chain L\n"
        result += "select cdr_h1, chain H and resi 26-32\n"
        result += "select cdr_h2, chain H and resi 52-56\n"
        result += "select cdr_h3, chain H and resi 95-102\n"
        result += "select cdr_l1, chain L and resi 24-34\n"
        result += "select cdr_l2, chain L and resi 50-56\n"
        result += "select cdr_l3, chain L and resi 89-97\n"
        result += "color red, cdr_h3\n"
        result += "color orange, cdr_h1 or cdr_h2 or cdr_l1 or cdr_l2 or cdr_l3\n"
    
    elif project.type == ProjectType.SMALL_MOLECULE:
        result += "show cartoon\n"
        result += "show sticks, chain X\n"
        result += "color green, chain X\n"
        result += "select binding_site, byres protein within 5 of chain X\n"
        result += "show sticks, binding_site\n"
        result += "color cyan, binding_site\n"
    
    # Chimera commands
    result += "\nChimera Commands:\n"
    result += f"open {model_path}\n"
    result += "style stick #0:.X\n"  # For ligand, if present
    
    if project.type == ProjectType.ANTIBODY:
        result += "color cyan #0:H@\n"
        result += "color pink #0:L@\n"
        result += "style ribbon #0\n"
    
    elif project.type == ProjectType.SMALL_MOLECULE:
        result += "style ribbon #0\n"
        result += "color green #0:.X@\n"
        result += "select #0:5.0@.X\n"
        result += "style stick sel\n"
        result += "color cyan sel\n"
    
    return result

@mcp.tool()
async def cleanup_project(
    project_id: str,
    remove_intermediates: bool = False,
    ctx: Context = None
) -> str:
    """Clean up a project by removing temporary files
    
    Args:
        project_id: ID of the project
        remove_intermediates: Whether to remove intermediate files (default: False)
        
    Returns:
        Status of the cleanup
    """
    if project_id not in active_projects:
        return f"Error: Project with ID '{project_id}' not found."
    
    project = active_projects[project_id]
    
    # Get list of files to keep
    files_to_keep = []
    
    # Always keep input files
    for key, path in project.input_files.items():
        if os.path.exists(path):
            files_to_keep.append(os.path.abspath(path))
    
    # Always keep best model/design and analysis files
    if "best_model" in project.output_files:
        if os.path.exists(project.output_files["best_model"]):
            files_to_keep.append(os.path.abspath(project.output_files["best_model"]))
    
    if "best_design" in project.output_files:
        if os.path.exists(project.output_files["best_design"]):
            files_to_keep.append(os.path.abspath(project.output_files["best_design"]))
    
    if "analysis" in project.output_files:
        if os.path.exists(project.output_files["analysis"]):
            files_to_keep.append(os.path.abspath(project.output_files["analysis"]))
    
    # Keep or remove intermediate files based on flag
    if not remove_intermediates:
        for key, path in project.output_files.items():
            if os.path.exists(path):
                files_to_keep.append(os.path.abspath(path))
    
    # Get all files in project directory
    all_files = []
    for root, dirs, files in os.walk(project.working_dir):
        for file in files:
            file_path = os.path.abspath(os.path.join(root, file))
            all_files.append(file_path)
    
    # Determine files to remove
    files_to_remove = [f for f in all_files if f not in files_to_keep]
    
    # Remove files
    removed_count = 0
    for file_path in files_to_remove:
        try:
            os.remove(file_path)
            removed_count += 1
        except OSError:
            continue
    
    return f"Cleanup completed. Removed {removed_count} files from project {project_id}."

# Load project state on startup
load_project_state()

# Run the server if this script is executed directly
if __name__ == "__main__":
    mcp.run()