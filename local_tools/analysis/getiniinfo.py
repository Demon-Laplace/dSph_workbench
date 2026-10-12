import os
import configparser

def read_model_info(num):

    current_dir = os.getcwd()
    modelname = os.path.basename(current_dir)

    folder_name = f"{modelname}{num}"
    ini_filename = f"IC_{modelname}{num}.ini"
    ini_path = os.path.join('./modelfile', ini_filename)
    orbit_path = os.path.join(folder_name, f'elinfo_{modelname}{num}.csv')

    if not os.path.exists(ini_path):
        print(f"File not found: {ini_path}")
        return
    if not os.path.exists(orbit_path):
        print(f"找不到文件: {orbit_path}")
        return None
    
    with open(orbit_path, 'r') as f:

        first_line = f.readline().strip()
        data = first_line.replace("# Orbit Data:", "").strip()
        numbers = [float(x) for x in data.split(",")]
     
        second_line = f.readline().strip()
        fd_value = float(second_line.replace("# Timescale_fd:", "").strip())

    config = configparser.ConfigParser()
    config.optionxform = str 
    config.read(ini_path)

    gas_mass = float(config['Density_gas']['mass'])
    gas_scaleRadius = float(config['Density_gas']['scaleRadius'])
    gas_scaleHeight = float(config['Density_gas']['scaleHeight'])
    gas_nbody = int(config['Component_gas']['nbody'])

    try:
        bulge_mass = float(config['Density_bulge']['mass'])
        halo_mass = float(config['Density_halo']['mass'])
        bulge_nbody = int(config['Component_bulge']['nbody'])
        halo_nbody = int(config['Component_halo']['nbody'])

        stellar_mass = bulge_mass + halo_mass
        stellar_nbody = bulge_nbody + halo_nbody
    
    except:
        stellar_mass = float(config['Density_stellar']['mass'])
        stellar_nbody = int(config['Component_stellar']['nbody'])


    print(f"=== {modelname}{num} ===\n")

    print("=== Component Gas Info ===")
    print(f"  Gas Mass: {gas_mass * 1e10:.2E}")
    print(f"  Gas ScaleHeight: {gas_scaleHeight} kpc")
    print(f"  Gas ScaleRadius: {gas_scaleRadius} kpc")
    print(f"  Gas Nbody: {gas_nbody}\n")

    print("=== Component Stellar Info ===")
    print(f"  Stellar Mass: {stellar_mass * 1e10:.2E}")
    print(f"  Stellar Nbody: {stellar_nbody}\n")

    print("=== Orbit Info ===")
    print(f"  X: {numbers[0]:6.1f}    Vx = {numbers[3]:6.2f}")
    print(f"  Y: {numbers[1]:6.1f}    Vy = {numbers[4]:6.2f}")
    print(f"  Z: {numbers[2]:6.1f}    Vz = {numbers[5]:6.2f}\n")

    print("=== Simulation Info ===")
    print(f"  Feedback Value: {fd_value:.2f}\n")

if __name__ == "__main__":
    num = input("Input model number：")
    read_model_info(num)
