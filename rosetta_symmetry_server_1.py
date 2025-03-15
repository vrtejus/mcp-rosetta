from mcp.server.fastmcp import FastMCP, Context
import os
import tempfile
import asyncio
import shutil
from typing import Optional, List, Dict
import glob

# Create an MCP server
mcp = FastMCP("RosettaSymmServer")

# Path to the Rosetta script - should be configurable through environment variables
ROSETTA_SCRIPT_PATH = os.environ.get(
    "ROSETTA_SCRIPT_PATH", 
    "/Users/ammachi/Developer/GitHub/MCP/mcp-rosetta/make_symmdef_file.pl"
)

# Default directories for input and output files
PDB_DIR = os.environ.get(
    "PDB_DIR",
    "/Users/ammachi/Developer/GitHub/MCP/mcp-rosetta/pdb_files"
)
OUTPUT_DIR = os.environ.get(
    "OUTPUT_DIR",
    "/Users/ammachi/Developer/GitHub/MCP/mcp-rosetta/output_files"
)

# Create directories if they don't exist
os.makedirs(PDB_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

@mcp.tool()
async def generate_symmetry_definition(
    pdb_file: str,
    mode: str = "NCS",
    chain_a: str = "A",
    chain_i: str = "B",
    ctx: Context = None
) -> str:
    """Generate a symmetry definition file for a protein structure using Rosetta's make_symmdef_file.pl.
    
    Args:
        pdb_file: The content of the PDB file as a string
        mode: Symmetry mode (NCS, CRYST, etc.)
        chain_a: The chain to keep
        chain_i: The chain to base symmetry off
        ctx: MCP context
        
    Returns:
        The symmetry definition file content
    """
    # Check if the Rosetta script exists
    if not os.path.exists(ROSETTA_SCRIPT_PATH):
        return f"Error: Rosetta script not found at {ROSETTA_SCRIPT_PATH}"
    
    # Validate input parameters
    if not pdb_file or not pdb_file.strip():
        return "Error: PDB file content is required"
    
    # Create a temporary directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Log progress
        if ctx:
            ctx.info(f"Created temporary directory for processing: {temp_dir}")
        
        # Write the PDB content to a temporary file
        input_pdb_path = os.path.join(temp_dir, "input.pdb")
        with open(input_pdb_path, "w") as f:
            f.write(pdb_file)
        
        if ctx:
            ctx.info(f"Saved PDB content to {input_pdb_path}")
        
        # Build the command
        cmd = [
            ROSETTA_SCRIPT_PATH,
            "-m", mode,
            "-a", chain_a,
            "-i", chain_i,
            "-p", input_pdb_path
        ]
        
        if ctx:
            ctx.info(f"Running command: {' '.join(cmd)}")
        
        try:
            # Run the process
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            # Wait for the process to complete
            stdout_bytes, stderr_bytes = await process.communicate()
            stdout = stdout_bytes.decode('utf-8')
            stderr = stderr_bytes.decode('utf-8')
            
            if process.returncode != 0:
                if ctx:
                    ctx.info(f"Error running Rosetta script: {stderr}")
                return f"Error running Rosetta script: {stderr}"
            
            # Check if stderr contains warning messages
            if stderr and stderr.strip():
                if ctx:
                    ctx.info(f"Warning from Rosetta script: {stderr}")
            
            # The symmetry definition is in stdout
            symmetry_def = stdout
            
            # List all files in the temp directory to see what was created
            files_created = [f for f in os.listdir(temp_dir) if f != os.path.basename(input_pdb_path)]
            
            if ctx:
                ctx.info("Successfully generated symmetry definition")
                ctx.info(f"Files created: {', '.join(files_created)}")
            
            # Create a formatted response
            response = f"Successfully generated symmetry definition for {os.path.basename(input_pdb_path)}.\n\n"
            response += f"Symmetry Definition:\n{symmetry_def}\n\n"
            
            if files_created:
                response += f"Additional files generated by Rosetta (not included in this response):\n"
                response += "\n".join(files_created)
            
            return response
        except Exception as e:
            if ctx:
                ctx.info(f"Error: {str(e)}")
            return f"Error generating symmetry definition: {str(e)}"

def is_text_file(file_path: str) -> bool:
    """Check if a file is a text file by reading a chunk and checking for null bytes."""
    try:
        with open(file_path, 'rb') as f:
            chunk = f.read(8192)
            return b'\0' not in chunk
    except:
        return False

@mcp.tool()
async def generate_symmetry_save_files(
    filename: str,
    mode: str = "NCS",
    chain_a: str = "A",
    chain_i: str = "B",
    output_prefix: str = None,
    ctx: Context = None
) -> str:
    """Generate symmetry definition and save all output files to the output directory.
    
    Args:
        filename: Name of the PDB file in the PDB_DIR
        mode: Symmetry mode (NCS, CRYST, etc.)
        chain_a: The chain to keep
        chain_i: The chain to base symmetry off
        output_prefix: Optional prefix for output files (defaults to PDB filename without extension)
        ctx: MCP context
        
    Returns:
        Summary of files created and their locations
    """
    # Check if the Rosetta script exists
    if not os.path.exists(ROSETTA_SCRIPT_PATH):
        return f"Error: Rosetta script not found at {ROSETTA_SCRIPT_PATH}"
    
    # Validate input file
    input_pdb_path = os.path.join(PDB_DIR, filename)
    if not os.path.exists(input_pdb_path):
        return f"Error: PDB file not found at {input_pdb_path}"
    
    # Determine output prefix
    if not output_prefix:
        output_prefix = os.path.splitext(filename)[0]
    
    # Create a temporary directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Log progress
        if ctx:
            ctx.info(f"Created temporary directory for processing: {temp_dir}")
        
        # Copy the PDB file to temp directory
        temp_pdb_path = os.path.join(temp_dir, "input.pdb")
        shutil.copy(input_pdb_path, temp_pdb_path)
        
        if ctx:
            ctx.info(f"Copied PDB file to {temp_pdb_path}")
        
        # Build the command
        cmd = [
            ROSETTA_SCRIPT_PATH,
            "-m", mode,
            "-a", chain_a,
            "-i", chain_i,
            "-p", temp_pdb_path
        ]
        
        if ctx:
            ctx.info(f"Running command: {' '.join(cmd)}")
        
        try:
            # Run the process
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            # Wait for the process to complete
            stdout_bytes, stderr_bytes = await process.communicate()
            stdout = stdout_bytes.decode('utf-8')
            stderr = stderr_bytes.decode('utf-8')
            
            if process.returncode != 0:
                if ctx:
                    ctx.info(f"Error running Rosetta script: {stderr}")
                return f"Error running Rosetta script: {stderr}"
            
            # Check if stderr contains warning messages
            if stderr and stderr.strip():
                if ctx:
                    ctx.info(f"Warning from Rosetta script: {stderr}")
            
            # The symmetry definition is in stdout
            symmetry_def = stdout
            
            # Save symmetry definition to file
            symm_file_path = os.path.join(OUTPUT_DIR, f"{output_prefix}.symm")
            with open(symm_file_path, 'w') as f:
                f.write(symmetry_def)
            
            # Copy all other generated files to output directory
            copied_files = [f"{output_prefix}.symm"]
            for file_path in glob.glob(os.path.join(temp_dir, "input*")):
                if os.path.basename(file_path) != "input.pdb":  # Skip the input file
                    file_name = os.path.basename(file_path)
                    new_name = file_name.replace("input", output_prefix)
                    dest_path = os.path.join(OUTPUT_DIR, new_name)
                    shutil.copy(file_path, dest_path)
                    copied_files.append(new_name)
            
            if ctx:
                ctx.info(f"Successfully generated symmetry files and saved to {OUTPUT_DIR}")
                ctx.info(f"Files created: {', '.join(copied_files)}")
            
            # Create a formatted response
            response = f"Successfully generated symmetry definition and saved output files:\n\n"
            
            for file_name in copied_files:
                response += f"- {file_name} → {os.path.join(OUTPUT_DIR, file_name)}\n"
            
            response += f"\nSymmetry definition was saved to: {symm_file_path}\n\n"
            response += f"First few lines of the symmetry definition:\n"
            response += "\n".join(symmetry_def.split("\n")[:10])
            response += "\n...\n"
            
            return response
        except Exception as e:
            if ctx:
                ctx.info(f"Error: {str(e)}")
            return f"Error generating symmetry definition: {str(e)}"

@mcp.tool()
async def get_specific_symmetry_file(
    filename: str,
    file_type: str,
    mode: str = "NCS",
    chain_a: str = "A",
    chain_i: str = "B",
    ctx: Context = None
) -> str:
    """Generate symmetry definition and return the content of a specific output file.
    
    Args:
        filename: Name of the PDB file in the PDB_DIR
        file_type: Type of file to retrieve (e.g., "symm", "INPUT", "model_AB", "symm.pdb")
        mode: Symmetry mode (NCS, CRYST, etc.)
        chain_a: The chain to keep
        chain_i: The chain to base symmetry off
        ctx: MCP context
        
    Returns:
        The content of the requested file
    """
    # Check if the Rosetta script exists
    if not os.path.exists(ROSETTA_SCRIPT_PATH):
        return f"Error: Rosetta script not found at {ROSETTA_SCRIPT_PATH}"
    
    # Validate input file
    input_pdb_path = os.path.join(PDB_DIR, filename)
    if not os.path.exists(input_pdb_path):
        return f"Error: PDB file not found at {input_pdb_path}"
    
    # Validate file type
    if not file_type or not file_type.strip():
        return "Error: file_type is required"
    
    # Create a temporary directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy the PDB file to temp directory
        temp_pdb_path = os.path.join(temp_dir, "input.pdb")
        shutil.copy(input_pdb_path, temp_pdb_path)
        
        # Build the command
        cmd = [
            ROSETTA_SCRIPT_PATH,
            "-m", mode,
            "-a", chain_a,
            "-i", chain_i,
            "-p", temp_pdb_path
        ]
        
        try:
            # Run the process
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            # Wait for the process to complete
            stdout_bytes, stderr_bytes = await process.communicate()
            stdout = stdout_bytes.decode('utf-8')
            
            if process.returncode != 0:
                return f"Error running Rosetta script: {stderr_bytes.decode('utf-8')}"
            
            # If the file type is "symm", return the symmetry definition from stdout
            if file_type.lower() == "symm":
                return f"Symmetry Definition File:\n{stdout}"
            
            # Look for the requested file
            file_patterns = [
                file_type,
                f"input_{file_type}",
                f"input{file_type}",
                f"{file_type}_input",
                f"{file_type}.input"
            ]
            
            # Check if any of the patterns match a file in the directory
            found_file = None
            for pattern in file_patterns:
                for fname in os.listdir(temp_dir):
                    if fname == pattern or pattern in fname:
                        found_file = os.path.join(temp_dir, fname)
                        break
                if found_file:
                    break
            
            # If the file was found, return its content
            if found_file and os.path.exists(found_file):
                if is_text_file(found_file):
                    with open(found_file, 'r') as f:
                        content = f.read()
                    
                    # Also save to output directory for future reference
                    output_name = os.path.basename(found_file).replace("input", os.path.splitext(filename)[0])
                    output_path = os.path.join(OUTPUT_DIR, output_name)
                    shutil.copy(found_file, output_path)
                    
                    return f"Content of {os.path.basename(found_file)} (saved to {output_path}):\n\n{content}"
                else:
                    # Also save to output directory for future reference
                    output_name = os.path.basename(found_file).replace("input", os.path.splitext(filename)[0])
                    output_path = os.path.join(OUTPUT_DIR, output_name)
                    shutil.copy(found_file, output_path)
                    
                    return f"File {os.path.basename(found_file)} is a binary file and cannot be displayed as text. It has been saved to {output_path}."
            
            # If file was not found, return a list of available files
            files_created = [f for f in os.listdir(temp_dir) if f != os.path.basename(temp_pdb_path)]
            return f"File '{file_type}' not found. Available files are:\n{', '.join(files_created)}"
        except Exception as e:
            if ctx:
                ctx.info(f"Error: {str(e)}")
            return f"Error retrieving symmetry file: {str(e)}"

@mcp.tool()
async def list_pdb_files(ctx: Context = None) -> str:
    """List all available PDB files in the input directory.
    
    Returns:
        A list of PDB files available for processing
    """
    try:
        # Get all PDB files in the directory
        pdb_files = [f for f in os.listdir(PDB_DIR) if f.endswith(('.pdb', '.PDB'))]
        
        if not pdb_files:
            return f"No PDB files found in {PDB_DIR}"
        
        # Create a formatted response
        response = f"Found {len(pdb_files)} PDB files in {PDB_DIR}:\n\n"
        
        for i, file_name in enumerate(sorted(pdb_files), 1):
            file_path = os.path.join(PDB_DIR, file_name)
            file_size = os.path.getsize(file_path)
            file_time = os.path.getmtime(file_path)
            
            # Format file size
            if file_size < 1024:
                size_str = f"{file_size} B"
            elif file_size < 1024 * 1024:
                size_str = f"{file_size / 1024:.1f} KB"
            else:
                size_str = f"{file_size / (1024 * 1024):.1f} MB"
            
            # Format file time
            time_str = os.path.getmtime(file_path)
            
            response += f"{i}. {file_name} ({size_str})\n"
        
        return response
    
    except Exception as e:
        if ctx:
            ctx.info(f"Error listing PDB files: {str(e)}")
        return f"Error listing PDB files: {str(e)}"

@mcp.tool()
async def list_output_files(ctx: Context = None) -> str:
    """List all available output files in the output directory.
    
    Returns:
        A list of output files that have been generated
    """
    try:
        # Get all files in the output directory
        output_files = os.listdir(OUTPUT_DIR)
        
        if not output_files:
            return f"No output files found in {OUTPUT_DIR}"
        
        # Group files by their prefix
        file_groups = {}
        for file_name in output_files:
            # Try to extract the prefix (part before the first underscore or dot)
            prefix = file_name.split('_')[0].split('.')[0]
            
            if prefix not in file_groups:
                file_groups[prefix] = []
            
            file_groups[prefix].append(file_name)
        
        # Create a formatted response
        response = f"Found {len(output_files)} output files in {OUTPUT_DIR}:\n\n"
        
        for prefix, files in sorted(file_groups.items()):
            response += f"Files for {prefix}:\n"
            
            for file_name in sorted(files):
                file_path = os.path.join(OUTPUT_DIR, file_name)
                file_size = os.path.getsize(file_path)
                
                # Format file size
                if file_size < 1024:
                    size_str = f"{file_size} B"
                elif file_size < 1024 * 1024:
                    size_str = f"{file_size / 1024:.1f} KB"
                else:
                    size_str = f"{file_size / (1024 * 1024):.1f} MB"
                
                response += f"  - {file_name} ({size_str})\n"
            
            response += "\n"
        
        return response
    
    except Exception as e:
        if ctx:
            ctx.info(f"Error listing output files: {str(e)}")
        return f"Error listing output files: {str(e)}"

if __name__ == "__main__":
    mcp.run()