from pymol import cmd
import os
import time
import threading

def monitor_pdb(pdb_path, object_name="structure", refresh_rate=1.0):
    # Convert refresh_rate to float (PyMOL passes it as string)
    try:
        refresh_rate = float(refresh_rate)
    except (ValueError, TypeError):
        print(f"Warning: Invalid refresh rate '{refresh_rate}', using 1.0 seconds")
        refresh_rate = 1.0
    
    # Expand tilde to full home directory path
    full_path = os.path.expanduser(pdb_path)
    
    print(f"Path: {full_path}")
    print(f"Object: {object_name}")
    print(f"Refresh: {refresh_rate} seconds")
    
    # First load the structure
    cmd.delete(object_name)
    cmd.load(full_path, object_name)
    print(f"Loaded {full_path} as {object_name}")
    print(f"Now monitoring for changes every {refresh_rate} seconds...")
    
    # Start monitoring in a background thread
    monitor_thread = threading.Thread(
        target=_monitor_file_thread,
        args=(full_path, object_name, refresh_rate),
        daemon=True
    )
    monitor_thread.start()
    return monitor_thread

def _monitor_file_thread(file_path, object_name, refresh_rate):
    last_modified = os.path.getmtime(file_path)
    print(f"Initial timestamp: {last_modified}")
    
    while True:
        time.sleep(refresh_rate)
        try:
            current_modified = os.path.getmtime(file_path)
            if current_modified > last_modified:
                print(f"File changed! Old: {last_modified}, New: {current_modified}")
                # File has changed - reload it
                cmd.delete(object_name)
                cmd.load(file_path, object_name)
                print(f"Reloaded {object_name} from {file_path}")
                last_modified = current_modified
        except Exception as e:
            print(f"Error checking file: {e}")

cmd.extend("monitor_pdb", monitor_pdb)