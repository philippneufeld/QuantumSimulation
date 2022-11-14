
import os
import socket

def get_miccell_folder() -> str:

    host = socket.gethostname()
    home = os.getenv("HOME")

    if host.startswith("calc") and len(host) == 5:
        host = "calcmaster"

    paths = {
        "": os.path.join(home, "remote_groups/MicCells"),
        "ludwigsburg": os.path.join(home, "MicCells"),
        "sost": os.path.join(home, "MicCells"),
        "PI5-calcmaster": "/mnt/groups/MicCells",
    }

    if host not in paths:
        host = ""
    
    return paths[host]

def get_default_app_folder(app_name: str) -> str:
    path = os.path.join(get_miccell_folder(), f"TraceGasSensing/Simulation/QNOSE_QuantumSimulation/{app_name}")
    os.makedirs(path, exist_ok=True)
    return path
