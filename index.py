import numpy as np
import math
import matplotlib.pyplot as plt

# Global variable declaration
dict_input = {}  # Dictionary of parameters from input.txt


# This function read input.txt file and assign parameters to dict_input
def read_input():
    f_input = open("input.txt", "r")
    data = f_input.readlines()
    f_input.close()
    # Finding parameters by split, and convert to integer/float
    dict_input['polymer_lengths'] = data[0][data[0].find('[') + 1:data[0].find(']')].split(',')
    for i in range(len(dict_input['polymer_lengths'])):
        dict_input['polymer_lengths'][i] = int(dict_input['polymer_lengths'][i])
    dict_input['l'] = float(data[1].strip().split('=')[1])
    dict_input['diameter'] = float(data[2].strip().split('=')[1])
    dict_input['dim'] = int(data[3].strip().split('=')[1])
    dict_input['angle_interval'] = float(data[4].strip().split('=')[1])
    dict_input['mean_err'] = float(data[5].strip().split('=')[1])
    dict_input['max_tries'] = int(data[6].strip().split('=')[1])


# This function simulate a single polymer
# Input-polymer length, dimension and the file name to enter the coordinates
# Output-radii of the polymer
def single_polymer_sim(polymer_length: int, dim: int, file_name: str):
    radii = 0
    theta_list_deg = []
    phi_list_deg = []

    # List of phi in degrees
    for i in np.arange(0, 360, dict_input['angle_interval']):
        phi_list_deg.append(i)

    # 3D model
    if dim == 3:
        # List of theta in degrees
        for i in np.arange(0, 180, dict_input['angle_interval']):
            theta_list_deg.append(i)

        polymer = [[0.0, 0.0, 0.0]]  # List of monomers coordinates(starting at origin)
        i = 0
        # Build the polymer, one monomer at a time
        while i < polymer_length:
            theta_deg = theta_list_deg[np.random.randint(0, len(theta_list_deg))]  # Chose random theta from theta list
            phi_deg = phi_list_deg[np.random.randint(0, len(phi_list_deg))]  # Chose random phi from phi list
            # Convert angles to radian
            theta_rad = (theta_deg * np.pi) / 180.0
            phi_rad = (phi_deg * np.pi) / 180.0
            last_x = polymer[len(polymer) - 1][0]
            last_y = polymer[len(polymer) - 1][1]
            last_z = polymer[len(polymer) - 1][2]
            # Build the next monomer by spherical coordinates
            x = last_x + dict_input['l'] * np.sin(theta_rad) * np.cos(phi_rad)
            y = last_y + dict_input['l'] * np.sin(theta_rad) * np.sin(phi_rad)
            z = last_z + dict_input['l'] * np.cos(theta_rad)
            if check_monomer_distance(polymer, dim, x, y, z):
                # If the new monomer doesn't collide the others monomers it will append to the list
                polymer.append([x, y, z])
                i += 1

        # Open new file or overwrite any existing file with the same name to export the coordinates
        f_coordinates = open(file_name, 'w')
        # Headline
        f_coordinates.write('X COOR  \t  Y COOR  \t  Z COOR \n')
        # Write coordinates with round it into five important digits after the dot
        for point in polymer:
            f_coordinates.write('{0:+.5f}\t{1:+.5f}\t{2:+.5f}\n'.format(point[0], point[1], point[2]))
        f_coordinates.close()

        # Calculate radii
        radii = math.sqrt(
            polymer[len(polymer) - 1][0] ** 2 + polymer[len(polymer) - 1][1] ** 2 + polymer[len(polymer) - 1][2] ** 2)


    # 2D model
    elif dim == 2:
        polymer = [[0.0, 0.0]]  # List of monomers coordinates(starting at origin)
        i = 0
        # Build the polymer, one monomer at a time
        while i < polymer_length:
            phi_deg = phi_list_deg[np.random.randint(0, len(phi_list_deg))]  # Chose random phi from phi list
            phi_rad = (phi_deg * np.pi) / 180.0  # Convert angle to radian
            last_x = polymer[len(polymer) - 1][0]
            last_y = polymer[len(polymer) - 1][1]
            # Build the next monomer by polar coordinates
            x = last_x + dict_input['l'] * np.cos(phi_rad)
            y = last_y + dict_input['l'] * np.sin(phi_rad)
            if check_monomer_distance(polymer, dim, x, y, 0):
                # If the new monomer doesn't collide others monomers it will append to the list
                polymer.append([x, y])
                i += 1

        # Open new file_name or overwrite any existing file with the same name to export the coordinates
        f_coordinates = open(file_name, 'w')
        # Headline
        f_coordinates.write('X COOR	 \t  Y COOR \n')
        # Write coordinates with round it into five important digits after the dot
        for point in polymer:
            f_coordinates.write('{0:+.5f}\t{1:+.5f}\n'.format(point[0], point[1]))
        f_coordinates.close()

        # Calculate radii
        radii = math.sqrt(polymer[len(polymer) - 1][0] ** 2 + polymer[len(polymer) - 1][1] ** 2)
    return radii


# This function check if new monomer collide the others monomers
# Input-list of existing monomers, dimension and the new monomer coordinates
# Output-return False if the new monomer is too close to other,and True if it does not
def check_monomer_distance(polymer: [], dim: int, x_i: float, y_i: float, z_i: float):
    if dim == 3:
        for monomer in polymer:
            if math.sqrt((x_i - monomer[0]) ** 2 + (y_i - monomer[1]) ** 2 + (z_i - monomer[2]) ** 2) < dict_input[
                'diameter']:
                return False

    elif dim == 2:
        for monomer in polymer:
            if math.sqrt((x_i - monomer[0]) ** 2 + (y_i - monomer[1]) ** 2) < dict_input['diameter']:
                return False
    return True


# This function append radius of specific polymer to a suitable text file
# Input-list of radius and polymer length
def write_radius(radius: list, polymer_length: int):
    f_name_output = "radii_" + str(dict_input['dim']) + "d_N" + str(polymer_length) + "_l" + str(
        dict_input['l']) + ".txt"
    # Opens a file for both appending and reading
    f_radii = open(f_name_output, "+a")
    f_radii.seek(0)
    # If file is empty-add headline
    if f_radii.readline(5) != 'Radii':
        f_radii.write('Radii\n')
    # Add radius to the file
    for radii in radius:
        f_radii.write('{:.5f}\n'.format(radii))
    f_radii.close()


# This function read radius in text file of specific polymer
# Input-polymer length
# Output-return list of radius
def read_radius(polymer_length: int):
    radius = []
    f_name_output = "radii_" + str(dict_input['dim']) + "d_N" + str(polymer_length) + "_l" + str(
        dict_input['l']) + ".txt"
    f_input = open(f_name_output, "r")
    data = f_input.readlines()
    f_input.close()
    for i in data[1:]:
        radius.append(float(i.strip()))
    return radius


######################################################### Main #########################################################
read_input()
mean_radius_list = []  # List of all the mean radius
# This loop find mean radii and build Histogram graph for each polymer length
for length in dict_input['polymer_lengths']:
    radius_list = []  # List of all radius
    mean_radius = 0
    stopping_criteria = False
    # Repeat until stopping criteria is True/number of maximum polymers is reached
    while not stopping_criteria:
        # Simulate single polymer and append the radii to the list
        radius_list.append(single_polymer_sim(length, dict_input['dim'], 'coordinates.txt'))
        new_mean_radius = sum(radius_list) / len(radius_list)  # Calculate new mean radius

        # If stopping criteria is True the loop will end
        if (abs((new_mean_radius - mean_radius) / new_mean_radius)) < dict_input['mean_err']:
            stopping_criteria = True
        # If number of maximum polymers is reached warning will be print and the loop will end
        elif len(radius_list) == dict_input['max_tries']:
            print('Warning: The Radius for polymer length {0} did not converge in {1} tries'.format(length,
                                                                                                    len(radius_list)))
            stopping_criteria = True
        else:
            mean_radius = new_mean_radius

    # Append radius to file
    write_radius(radius_list, length)
    # Read all radius in the file
    radius_list = read_radius(length)
    mean_radius_list.append(sum(radius_list) / len(radius_list))  # Calculate mean radius

    # Histogram graph of all the radii from the file
    plt.figure()
    plt.hist(radius_list)
    plt.xlabel('Radii')
    plt.ylabel('Probability')
    plt.title('Histogram for: dim-{} N={} l={}'.format(dict_input['dim'], length, dict_input['l']))

# Graph Mean Radius Vs. Polymers Length
if len(dict_input['polymer_lengths']) > 1:
    plt.figure()
    plt.plot(dict_input['polymer_lengths'], mean_radius_list)
elif len(dict_input['polymer_lengths']) == 1:
    plt.figure()
    plt.plot(dict_input['polymer_lengths'], mean_radius_list, 'bo')
plt.xlabel('Polymer Length')
plt.ylabel('Mean Radii')
plt.title('Mean Radius Vs. Polymers Length for: dim-{} l={}'.format(dict_input['dim'], dict_input['l']))
plt.grid()
plt.show()

# Print when the program finished
print('program end')
