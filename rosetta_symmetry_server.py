from mcp.server.fastmcp import FastMCP, Context
import os
import tempfile
import asyncio
from typing import Optional

# Create an MCP server
mcp = FastMCP("RosettaSymmServer")

# Path to the Rosetta script - should be configurable through environment variables
ROSETTA_SCRIPT_PATH = os.environ.get(
    "ROSETTA_SCRIPT_PATH", 
    "/Users/ammachi/Developer/GitHub/MCP/mcp-rosetta/make_symmdef_file.pl"
)

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
async def get_specific_symmetry_file(
    pdb_file: str,
    file_name: str,
    mode: str = "NCS",
    chain_a: str = "A",
    chain_i: str = "B",
    ctx: Context = None
) -> str:
    """Generate and retrieve a specific file from the symmetry definition process.
    
    Args:
        pdb_file: The content of the PDB file as a string
        file_name: The name of the file to retrieve (e.g., "input_symm.pdb")
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
    
    # Validate input parameters
    if not pdb_file or not pdb_file.strip():
        return "Error: PDB file content is required"
    
    if not file_name or not file_name.strip():
        return "Error: file_name is required"
    
    # Create a temporary directory for processing
    with tempfile.TemporaryDirectory() as temp_dir:
        # Write the PDB content to a temporary file
        input_pdb_path = os.path.join(temp_dir, "input.pdb")
        with open(input_pdb_path, "w") as f:
            f.write(pdb_file)
        
        # Build the command
        cmd = [
            ROSETTA_SCRIPT_PATH,
            "-m", mode,
            "-a", chain_a,
            "-i", chain_i,
            "-p", input_pdb_path
        ]
        
        try:
            # Run the process
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            # Wait for the process to complete
            stdout_bytes, _ = await process.communicate()
            
            if process.returncode != 0:
                return "Error running Rosetta script"
            
            # The base name of the PDB file
            base_name = os.path.basename(input_pdb_path).split('.')[0]
            
            # List of possible file patterns to check
            file_patterns = [
                file_name,
                f"{base_name}_{file_name}",
                f"{base_name}{file_name}",
                f"{file_name}_{base_name}",
                f"{file_name}.{base_name}"
            ]
            
            # Check if any of the patterns match a file in the directory
            found_file = None
            for pattern in file_patterns:
                for filename in os.listdir(temp_dir):
                    if filename == pattern or pattern in filename:
                        found_file = os.path.join(temp_dir, filename)
                        break
                if found_file:
                    break
            
            # If the file is the symm file (from stdout), return that
            if file_name.endswith('.symm') or file_name == 'symm':
                return f"Symmetry Definition File ({file_name}):\n{stdout_bytes.decode('utf-8')}"
            
            # If the file was found, return its content
            if found_file and os.path.exists(found_file):
                if is_text_file(found_file):
                    with open(found_file, 'r') as f:
                        content = f.read()
                    return f"Content of {os.path.basename(found_file)}:\n{content}"
                else:
                    return f"File {os.path.basename(found_file)} is a binary file and cannot be displayed as text"
            
            # If file was not found, return a list of available files
            files_created = [f for f in os.listdir(temp_dir) if f != os.path.basename(input_pdb_path)]
            return f"File '{file_name}' not found. Available files are:\n{', '.join(files_created)}"
        except Exception as e:
            if ctx:
                ctx.info(f"Error: {str(e)}")
            return f"Error retrieving symmetry file: {str(e)}"

if __name__ == "__main__":
    mcp.run()