import numpy as np
import os

BOHR_TO_ANGSTROM = 0.529177210903

def parse_cube_file(filename):
    """Parse the CUBE file to extract voxel density data and atomic positions."""
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    atom_info_start = 2
    
    # Read the number of atoms and the origin of the grid
    natoms, origin_x, origin_y, origin_z = map(float, lines[atom_info_start].split())
    natoms = int(natoms)
    
    # Read the grid information
    nx, *vx = map(float, lines[atom_info_start + 1].split())
    ny, *vy = map(float, lines[atom_info_start + 2].split())
    nz, *vz = map(float, lines[atom_info_start + 3].split())
    nx, ny, nz = int(nx), int(ny), int(nz)
    
    # Parse atom positions and skip to voxel data
    atom_lines = lines[atom_info_start + 4:atom_info_start + 4 + natoms]
    atom_positions = []
    atomic_numbers = []
    for line in atom_lines:
        atomic_number, _, x, y, z = map(float, line.split())
        atomic_numbers.append(int(atomic_number))
        atom_positions.append((x, y, z))
    
    # Read voxel densities
    density_data = []
    density_lines = lines[atom_info_start + 4 + natoms:]
    for line in density_lines:
        density_data.extend(map(float, line.split()))
    density_data = np.array(density_data).reshape((nx, ny, nz))
    
    return {
        'natoms': natoms,
        'origin': np.array([origin_x, origin_y, origin_z]) * BOHR_TO_ANGSTROM,
        'grid': np.array([nx, ny, nz]),
        'voxel_vectors': np.array([vx, vy, vz]) * BOHR_TO_ANGSTROM,
        'atom_positions': np.array(atom_positions) * BOHR_TO_ANGSTROM,
        'atomic_numbers': atomic_numbers,
        'density': density_data,
    }

def shift_density_to_center(density, grid_shape):
    """Shift the maximum density to the center of the grid."""
    nx, ny, nz = grid_shape
    max_index = np.unravel_index(np.argmax(density), density.shape)
    center_index = np.array([nx // 2, ny // 2, nz // 2])
    shift = center_index - np.array(max_index)
    density = np.roll(density, shift, axis=(0, 1, 2))
    return density, shift

def shift_gyration_center_back(gyration_center, shift, voxel_vectors):
    """Shift the gyration center back to its original position."""
    shift_vector = np.dot(shift, voxel_vectors)
    return gyration_center - shift_vector

def calculate_gyration_properties_with_shift(parsed_cube):
    """Calculate gyration center and radius with boundary adjustment."""
    density = parsed_cube['density']
    grid_shape = parsed_cube['grid']
    
    # Shift maximum density to the center
    shifted_density, shift = shift_density_to_center(density, grid_shape)
    
    # Calculate gyration properties
    origin = parsed_cube['origin']
    voxel_vectors = parsed_cube['voxel_vectors']
    nx, ny, nz = grid_shape
    
    # Compute grid points in Cartesian coordinates
    x = origin[0] + np.arange(nx)[:, None, None] * voxel_vectors[0][0]
    y = origin[1] + np.arange(ny)[None, :, None] * voxel_vectors[1][1]
    z = origin[2] + np.arange(nz)[None, None, :] * voxel_vectors[2][2]
    
    # Create meshgrid for grid points
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Compute the density-weighted center of gyration
    total_density = np.sum(shifted_density)
    gx = np.sum(shifted_density * X) / total_density
    gy = np.sum(shifted_density * Y) / total_density
    gz = np.sum(shifted_density * Z) / total_density
    gyration_center = np.array([gx, gy, gz])
    
    # Compute the radius of gyration
    squared_distances = (X - gx)**2 + (Y - gy)**2 + (Z - gz)**2
    gyration_radius = np.sqrt(np.sum(shifted_density * squared_distances) / total_density)
    
    # Shift the gyration center back
    shifted_back_center = shift_gyration_center_back(gyration_center, shift, voxel_vectors)
    
    return shifted_back_center, gyration_radius

def process_trajectory(file_prefix, start, end, step, output_file):
    """Process cube files and write a combined XYZ trajectory file."""
    element_mapping = {1: "H", 8: "O", 19: "K"}
    
    with open(output_file, 'w') as xyz_file:
        for i in range(start, end + 1, step):
            filename = f"{file_prefix}{i}.cube"
            if not os.path.exists(filename):
                print(f"File {filename} not found. Skipping.")
                continue
            
            parsed_cube = parse_cube_file(filename)
            gyration_center, gyration_radius = calculate_gyration_properties_with_shift(parsed_cube)
            atom_positions = parsed_cube['atom_positions']
            atomic_numbers = parsed_cube['atomic_numbers']
            
            total_atoms = len(atom_positions) + 1
            xyz_file.write(f"{total_atoms}\n")
            xyz_file.write(f"Gyration Radius: {gyration_radius:.4f} Angstrom\n")
            
            for index, pos in enumerate(atom_positions):
                element = element_mapping.get(atomic_numbers[index], "C")
                xyz_file.write(f"{element} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
            
            xyz_file.write(f"X {gyration_center[0]:.6f} {gyration_center[1]:.6f} {gyration_center[2]:.6f}\n")

# Parameters
file_prefix = "wat-MD-SPIN_DENSITY-1_"
start = 0 
end = 5400 
step = 40
output_xyz_path = "trajectory_with_gyration_center.xyz"

# Process trajectory
process_trajectory(file_prefix, start, end, step, output_xyz_path)

print(f"Combined trajectory with gyration center written to {output_xyz_path}")
