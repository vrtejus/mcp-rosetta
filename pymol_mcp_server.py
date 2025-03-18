#!/usr/bin/env python3
"""
PyMOL MCP Server

This server implements the Model Context Protocol (MCP) for PyMOL integration.
It connects to the PyMOL socket plugin and allows sending PyMOL commands.
"""

from mcp.server.fastmcp import FastMCP, Context
import socket
import json
import logging
import asyncio
from contextlib import asynccontextmanager
from typing import AsyncIterator, Dict, Any, Optional

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("PyMOLMCPServer")

class PyMOLConnection:
    """Class to manage connection to the PyMOL socket plugin"""
    
    def __init__(self, host: str = 'localhost', port: int = 9876):
        self.host = host
        self.port = port
        self.sock: Optional[socket.socket] = None
        
    def connect(self) -> bool:
        """Connect to the PyMOL socket plugin"""
        if self.sock:
            return True
            
        try:
            self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            self.sock.connect((self.host, self.port))
            logger.info(f"Connected to PyMOL at {self.host}:{self.port}")
            return True
        except Exception as e:
            logger.error(f"Failed to connect to PyMOL: {str(e)}")
            self.sock = None
            return False
    
    def disconnect(self):
        """Disconnect from the PyMOL socket plugin"""
        if self.sock:
            try:
                self.sock.close()
            except Exception as e:
                logger.error(f"Error disconnecting from PyMOL: {str(e)}")
            finally:
                self.sock = None

    def send_command(self, code: str) -> Dict[str, Any]:
        """Send a PyMOL command to the socket plugin"""
        if not self.sock and not self.connect():
            raise ConnectionError("Not connected to PyMOL")
        
        command = {
            "type": "pymol_command",
            "code": code
        }
        
        try:
            logger.info(f"Sending command to PyMOL")
            self.sock.sendall(json.dumps(command).encode('utf-8'))
            
            # Wait for response
            self.sock.settimeout(10.0)  # 10 second timeout
            response_data = self.sock.recv(4096)
            
            if not response_data:
                logger.error("No response received from PyMOL")
                self.sock = None
                raise ConnectionError("Connection closed by PyMOL")
            
            response = json.loads(response_data.decode('utf-8'))
            return response
        except socket.timeout:
            logger.error("Socket timeout while waiting for response from PyMOL")
            self.sock = None
            raise Exception("Timeout waiting for PyMOL response")
        except ConnectionError as e:
            logger.error(f"Socket connection error: {str(e)}")
            self.sock = None
            raise
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON response from PyMOL: {str(e)}")
            raise Exception(f"Invalid response from PyMOL: {str(e)}")
        except Exception as e:
            logger.error(f"Error sending command to PyMOL: {str(e)}")
            self.sock = None
            raise

# Global connection for resources and tools
_pymol_connection = None

def get_pymol_connection() -> PyMOLConnection:
    """Get or create a persistent PyMOL connection"""
    global _pymol_connection
    
    # If we have an existing connection, check if it's still valid
    if _pymol_connection is not None:
        try:
            # Check connection by sending a ping command
            _pymol_connection.send_command("pass")  # Using a no-op command
            return _pymol_connection
        except Exception as e:
            # Connection is dead, close it and create a new one
            logger.warning(f"Existing connection is no longer valid: {str(e)}")
            try:
                _pymol_connection.disconnect()
            except:
                pass
            _pymol_connection = None
    
    # Create a new connection if needed
    if _pymol_connection is None:
        _pymol_connection = PyMOLConnection()
        if not _pymol_connection.connect():
            logger.error("Failed to connect to PyMOL")
            _pymol_connection = None
            raise Exception("Could not connect to PyMOL. Make sure the PyMOL socket plugin is running.")
        logger.info("Created new persistent connection to PyMOL")
    
    return _pymol_connection

@asynccontextmanager
async def server_lifespan(server: FastMCP) -> AsyncIterator[Dict[str, Any]]:
    """Manage server startup and shutdown lifecycle"""
    try:
        logger.info("PyMOL MCP server starting up")
        
        # Try to connect to PyMOL on startup to verify it's available
        try:
            pymol = get_pymol_connection()
            logger.info("Successfully connected to PyMOL on startup")
        except Exception as e:
            logger.warning(f"Could not connect to PyMOL on startup: {str(e)}")
            logger.warning("Make sure the PyMOL plugin is running and socket server is started")
        
        # Return an empty context - we're using the global connection
        yield {}
    finally:
        # Clean up the global connection on shutdown
        global _pymol_connection
        if _pymol_connection:
            logger.info("Disconnecting from PyMOL on shutdown")
            _pymol_connection.disconnect()
            _pymol_connection = None
        logger.info("PyMOL MCP server shut down")

# Create the MCP server with lifespan support
mcp = FastMCP(
    "PyMOLMCP",
    description="PyMOL integration through the Model Context Protocol",
    lifespan=server_lifespan
)

# Tool: Execute a PyMOL command
@mcp.tool()
def execute_pymol(ctx: Context, code: str) -> str:
    """
    Execute a PyMOL command.
    
    Parameters:
    - code: The PyMOL Python code to execute
    """
    try:
        pymol = get_pymol_connection()
        response = pymol.send_command(code)
        
        if response.get("status") == "success":
            return f"PyMOL command executed successfully"
        else:
            return f"Error executing PyMOL command: {response.get('message', 'Unknown error')}"
    except Exception as e:
        logger.error(f"Error in execute_pymol tool: {str(e)}")
        return f"Error executing PyMOL command: {str(e)}"

# Tool: Get PyMOL version
@mcp.tool()
def get_pymol_version(ctx: Context) -> str:
    """Get the version of PyMOL currently running"""
    try:
        pymol = get_pymol_connection()
        code = "import pymol; print(pymol.__version__)"
        response = pymol.send_command(code)
        
        if response.get("status") == "success":
            return f"Command sent to get PyMOL version. Check PyMOL console for output."
        else:
            return f"Error getting PyMOL version: {response.get('message', 'Unknown error')}"
    except Exception as e:
        logger.error(f"Error in get_pymol_version tool: {str(e)}")
        return f"Error getting PyMOL version: {str(e)}"

# Tool: Load a PDB file
@mcp.tool()
def load_pdb(ctx: Context, pdb_id: str) -> str:
    """
    Load a PDB file from the PDB database.
    
    Parameters:
    - pdb_id: The 4-letter PDB ID to load
    """
    try:
        pymol = get_pymol_connection()
        # Create PyMOL code to fetch the PDB file
        code = f"cmd.fetch('{pdb_id}', async_=0)"
        response = pymol.send_command(code)
        
        if response.get("status") == "success":
            return f"PDB {pdb_id} loaded successfully"
        else:
            return f"Error loading PDB {pdb_id}: {response.get('message', 'Unknown error')}"
    except Exception as e:
        logger.error(f"Error in load_pdb tool: {str(e)}")
        return f"Error loading PDB {pdb_id}: {str(e)}"
# Tool: Apply a visualization preset
@mcp.tool()
def apply_visualization(
    ctx: Context,
    style: str = "cartoon",
    color_scheme: str = "spectrum"
) -> str:
    """
    Apply a visualization preset to the current PyMOL session.
    
    Parameters:
    - style: Visualization style (cartoon, ribbon, surface, spheres, lines, sticks)
    - color_scheme: Color scheme (spectrum, chainbows, rainbow, bfactors)
    """
    try:
        pymol = get_pymol_connection()
        
        # Create code for the style
        style_code = ""
        if style == "cartoon":
            style_code = "cmd.hide('everything'); cmd.show('cartoon')"
        elif style == "ribbon":
            style_code = "cmd.hide('everything'); cmd.show('ribbon')"
        elif style == "surface":
            style_code = "cmd.hide('everything'); cmd.show('surface')"
        elif style == "spheres":
            style_code = "cmd.hide('everything'); cmd.show('spheres')"
        elif style == "lines":
            style_code = "cmd.hide('everything'); cmd.show('lines')"
        elif style == "sticks":
            style_code = "cmd.hide('everything'); cmd.show('sticks')"
        else:
            return f"Unknown style: {style}"
        
        # Create code for the color scheme
        color_code = ""
        if color_scheme == "spectrum":
            color_code = "cmd.spectrum('count', selection='all')"
        elif color_scheme == "chainbows":
            color_code = "cmd.util.chainbow('all')"
        elif color_scheme == "rainbow":
            color_code = "cmd.spectrum('count', 'rainbow', selection='all')"
        elif color_scheme == "bfactors":
            color_code = "cmd.spectrum('b', 'blue_white_red', selection='all')"
        else:
            return f"Unknown color scheme: {color_scheme}"
        
        # Combine and send the code
        code = f"{style_code}; {color_code}"
        response = pymol.send_command(code)
        
        if response.get("status") == "success":
            return f"Applied {style} style with {color_scheme} coloring"
        else:
            return f"Error applying visualization: {response.get('message', 'Unknown error')}"
    except Exception as e:
        logger.error(f"Error in apply_visualization tool: {str(e)}")
        return f"Error applying visualization: {str(e)}"

# Prompt templates
@mcp.prompt()
def pymol_tips() -> str:
    """Guidance for executing PyMOL commands via MCP"""
    return """When working with the PyMOL MCP integration, keep these tips in mind:

1. Always check if the PyMOL socket plugin is running before sending commands.

2. The execute_pymol tool allows you to run any valid PyMOL Python code. Here are some examples:
   - cmd.load('path/to/file.pdb')
   - cmd.show('cartoon')
   - cmd.color('red', 'resi 10')

3. For common operations, use the specialized tools like:
   - load_pdb: Load PDB files directly from the PDB database
   - apply_visualization: Apply predefined visualization styles

4. PyMOL commands are executed in the PyMOL environment and should use the cmd module.

5. Check the PyMOL console for any output produced by your commands.

6. When generating PyMOL code, ensure it's valid Python that imports necessary modules and uses proper syntax.
"""

def main():
    """Run the MCP server"""
    mcp.run()

if __name__ == "__main__":
    main()