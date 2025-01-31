import subprocess
import psutil
import time
from datetime import datetime
import csv
import signal
import sys
import GPUtil

def get_gpu_info():
    try:
        gpus = GPUtil.getGPUs()
        if gpus:
            gpu = gpus[0]  # Using first GPU
            return {
                'gpu_memory_used': gpu.memoryUsed / 1024,  # Convert to GB
                'gpu_memory_total': gpu.memoryTotal / 1024,
                'gpu_utilization': gpu.load * 100
            }
    except:
        return None

def monitor_container(container_process, output_file):
    with open(output_file, 'w', newline='') as f:
        fieldnames = ['timestamp', 'cpu_percent', 'ram_used_gb', 'ram_total_gb',
                     'gpu_memory_used_gb', 'gpu_memory_total_gb', 'gpu_utilization']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        while container_process.poll() is None:
            try:
                cpu_percent = psutil.Process(container_process.pid).cpu_percent()
                ram = psutil.virtual_memory()
                gpu_info = get_gpu_info()
                
                metrics = {
                    'timestamp': datetime.now().isoformat(),
                    'cpu_percent': cpu_percent,
                    'ram_used_gb': ram.used / (1024**3),
                    'ram_total_gb': ram.total / (1024**3)
                }
                
                if gpu_info:
                    metrics.update(gpu_info)
                
                writer.writerow(metrics)
                f.flush()  # Ensure writing to file
                time.sleep(1)
                
            except psutil.NoSuchProcess:
                break
            except KeyboardInterrupt:
                break

def main():
    # Launch singularity container
    cmd = [
        "singularity",
        "exec",
        "--nv",
        "./boltz-cuda124-latest.sif",
        "boltz",
        "predict",
        "--cache",
        "/home/ariel/runtime/databases/boltz",
        "--output_format",
        "pdb",
        "--use_msa_server",
        "--diffusion_samples",
        "5",
        "./examples/ccr1.fasta"
    ]
    container = subprocess.Popen(cmd)
    
    try:
        monitor_container(container, 'container_metrics.csv')
    finally:
        # Cleanup
        if container.poll() is None:
            container.terminate()
            container.wait(timeout=5)

if __name__ == "__main__":
    main()