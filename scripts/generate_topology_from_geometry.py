import math
import sys

def get_distance(c1, c2):
    return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)

def main(pdbqt_file, output_itp):
    atoms = []
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = name[0]
                if name[0].isdigit(): element = name[1]
                mass = 12.01
                if element == 'H': mass = 1.008
                elif element == 'N': mass = 14.01
                elif element == 'O': mass = 16.00
                elif element == 'S': mass = 32.06
                elif element == 'F': mass = 19.00
                elif element == 'C' and len(name) > 1 and name[1] == 'l': mass = 35.45 # Cl
                
                atoms.append({'name': name, 'mass': mass, 'coord': (x/10.0, y/10.0, z/10.0)})

    bonds = []
    THRESHOLD = 0.17 
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            dist = get_distance(atoms[i]['coord'], atoms[j]['coord'])
            if dist < THRESHOLD:
                bonds.append((i+1, j+1, dist))

    with open(output_itp, 'w') as f:
        f.write('[ moleculetype ]\nUNL 3\n\n[ atoms ]\n')
        for i, a in enumerate(atoms):
            atype = 'HA' if a['mass'] < 2.0 else 'CA'
            f.write('{:4d} {:4s} 1 UNL {:4s} {:4d} 0.00 {:.3f}\n'.format(i+1, atype, a['name'], i+1, a['mass']))
        f.write('\n[ bonds ]\n')
        for b in bonds:
            f.write('{:4d} {:4d} 1 {:.5f} 150000.0\n'.format(b[0], b[1], b[2]))
    print('Done. {} atoms, {} bonds.'.format(len(atoms), len(bonds)))

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
