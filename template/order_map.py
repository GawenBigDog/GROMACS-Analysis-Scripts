#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2017.03.07 - A. Lorkowski - Changed script to provide order parameters for individual lipids.
#                           - Updated code with current Martini standard.

# 2011.11.27 - Helgi I. Ingolfsson - Fix POPC and POPE (tail order is flipped in regular Martini 2.1 itp)

from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess
import numpy as np

### SOME TOOLS

# parse a .gro file
# return a list of coordinates
def read_gro(file, atoms):
  line_counter = 0
  number_of_particles = 0
  first, second = [], []
  for line in open(file):
    if line_counter == 1:
      number_of_particles = int(line)
    elif line_counter > 1 and line_counter < number_of_particles + 2:
      if line[10:15].strip() == atoms[0]:
        first.append([float(line[21:30]), float(line[31:40]), float(line[41:50])])
      elif line[10:15].strip() == atoms[1]:
        second.append([float(line[21:30]), float(line[31:40]), float(line[41:50])])
    line_counter += 1
  return [first, second]

def get_coord(file, atom):
  line_counter = 0
  number_of_particles = 0
  position = []
  for line in open(file):
    if line_counter == 1:
      number_of_particles = int(line)
    elif line_counter > 1 and line_counter < number_of_particles + 2:
      if line[10:15].strip() == atom:
        position.append([float(line[21:30]), float(line[31:40]), float(line[41:50])])
    line_counter += 1
  return position

### REAL STUFF

if len(argv) != 10:
  # coments/usage
  print '''
  Compute (second rank) order parameter, defined as:

    P2 = 0.5*(3*<cos²(theta)> - 1)

  where "theta" is the angle between the bond and the bilayer normal.
  P2 = 1      perfect alignement with the bilayer normal
  P2 = -0.5   anti-alignement
  P2 = 0      random orientation

  All lipids defined in the "martini_v2.0_lipids.itp" file can be analyzed
  with this script.
  Usage: %s <traj file> <tpr file> <initial time> <final time> <skip frames> <bilayer normal - xyz> <lipid type>

    > %s traj.xtc topol.tpr 0 10000 5 0 0 1 DSPC

  will for example read a 10ns trajectory of 64 DSPC lipids, calculating the order parameter for 
  every 5th frame and averaging the results. P2 will be calculated relative to the y-axis.

  WARNING script will output all frames in one go, into files called frame_dump_XXX.gro and 
  then remove them so don't have any other files with this name in the current directory.
  ''' % (argv[0], argv[0])
  exit(0)

else:

  # snapshots
  trajfile = argv[1]
  tprfile = argv[2]
  initial_time = int(argv[3])
  final_time = int(argv[4])
  traj_skip = int(argv[5])
  # (normalized) orientation of bilayer normal
  orientation_of_bilayer_normal = [float(argv[6]), float(argv[7]), float(argv[8])]
  norm = sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
  for i in range(3):
    orientation_of_bilayer_normal[i] /= norm
  stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
    orientation_of_bilayer_normal[0], \
    orientation_of_bilayer_normal[1], \
    orientation_of_bilayer_normal[2]  \
  ))
  # number of lipids
  # number_of_lipids = int(argv[8])
  # lipid type
  lipid_type = argv[9]

  # output legend
  phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
  phosphatidicacid_bond_names = " PO4-GL1 GL1-GL2 "
  phosphatidylserine_bond_names = " CNO-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylinositol_bond_names = " C1-C2 C1-C3 C2-C3 C1-PO4 PO4-GL1 GL1-GL2 "
  lysophosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
  sphingomyelin_bond_names = " NC3-PO4 PO4-AM1 AM1-AM2 "
  ceramide_bond_names = " AM1-AM2 "

  # PCs
  if   lipid_type == "DAPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-D1B C1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"; hbCount = 3
  elif lipid_type == "DEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"; hbCount = 3
  elif lipid_type == "DHPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"; hbCount = 3
  elif lipid_type == "DLPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DOPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"; hbCount = 3
  elif lipid_type == "DPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "DSPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"; hbCount = 3
  elif lipid_type == "PAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PGPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PIPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "POPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PUPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PRPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  
  # PEs
  elif lipid_type == "DHPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C1B-C2B\n"; hbCount = 3
  elif lipid_type == "DLPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"; hbCount = 3
  elif lipid_type == "DOPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-D3A D3A-C4A C4A-C5A C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"; hbCount = 3
  elif lipid_type == "DSPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C4A-C5A C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"; hbCount = 3
  elif lipid_type == "DPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PAPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PGPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PIPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "POPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PQPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PRPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PUPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  
  # PPCS
  elif lipid_type == "PPCS": bond_names = " NC3-PO4 PO4-AM1 AM1-AM2 AM1-C1A GL2-D1B C1A-C2A C2A-C3A C3A-C4A D1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3

  # PAs
  elif lipid_type == "POPA": bond_names = phosphatidicacid_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 2
  elif lipid_type == "PGPA": bond_names = phosphatidicacid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 2

  # PSs
  elif lipid_type == "PGPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PRPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3

  # PIs
  elif lipid_type == "PAPI": bond_names = phosphatidylinositol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 6
  elif lipid_type == "POPI": bond_names = phosphatidylinositol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 6
  elif lipid_type == "PVPI": bond_names = phosphatidylinositol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 6

  # PPC
  elif lipid_type == "PPC": bond_names =  lysophosphatidylcholine_bond_names + "GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3

  # SMs
  elif lipid_type == "DPSM": bond_names =  sphingomyelin_bond_names + "AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PGSM": bond_names =  sphingomyelin_bond_names + "AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"; hbCount = 3
  elif lipid_type == "PNSM": bond_names =  sphingomyelin_bond_names + "AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"; hbCount = 3
  elif lipid_type == "POSM": bond_names =  sphingomyelin_bond_names + "AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-D2B D2B-C3B C3B-C4B\n"; hbCount = 3
  elif lipid_type == "PVSM": bond_names =  sphingomyelin_bond_names + "AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-D3B D3B-C4B\n"; hbCount = 3

  # CEs
  elif lipid_type == "DPCE": bond_names =  ceramide_bond_names + "AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-C3B C3B-C4B\n"; hbCount = 1
  
  # output legend
  output_legend = "  Frame" + bond_names 

  # write the stuff
  stdout.write("\n " + output_legend)
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n")
  #output = open('order.dat', 'w')
  #output.write(output_legend)
  #output.write(("-"*(len(output_legend) - 1)) + "\n")

  # Output all frame using trjconv 
  stdout.write("Output all coordinate files \n")
  command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o frame_dump_.gro > /dev/null" % (lipid_type, trajfile, tprfile, initial_time, final_time, traj_skip)
  print command
  subprocess.call(command, shell=True)

  command = "grep %s frame_dump_0.gro | grep PO4 | wc -l" % (lipid_type)

  number_of_lipids = int(subprocess.check_output(command, shell=True))

  # For each dumped frame
  stdout.write("Starting P2 calculation")
  byFrame_order_parameters = []
  byFrame_position = []
  file_count = 0
  bonds = []
  
  while True:
    filename = "frame_dump_" + str(file_count) + ".gro"
    if not path.isfile(filename) or path.getsize(filename) == 0:
        break
    
    stdout.write("Taking care of snapshot %s \n" % filename)

    # compute order parameter for each bond, for each snapshot
    #current_order_parameters = []
    
    byRes_order_parameters = []

    # Store each frame's order parameters in an array
    #byFrame_order_parameters = []
    
    # bonds respectively involved in the head,
    #                             in the junction head-tail,
    #                             in each tail
    bonds = []

    for bond_name in bond_names.split():
      bonds.append(bond_name.split("-"))

    for i in range(number_of_lipids):
      current_order_parameters = []
      for bond in bonds:
        # parse .gro file, grep bead coordinates
        first, second = read_gro(filename, bond)

        # compute order parameter for each lipid
        order_parameter = 0.0
        # vector between the two previous beads (orientation doesn't matter)
        vector = [0.0, 0.0, 0.0]
        for j in range(3):
          vector[j] = first[i][j] - second[i][j]
        norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
        # compute projection on the bilayer normal
        projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
        # order parameter
        order_parameter = projection**2/norm2
        current_order_parameters.append(0.5*(3.0*(order_parameter) - 1.0))

      # compute final averaged order parameter
      # store everything in lists
      # current_order_parameters.append(0.5*(3.0*(order_parameter/number_of_lipids) - 1.0))
      byRes_order_parameters.append(current_order_parameters)

    atomRef = "PO4"
    position = get_coord(filename, atomRef)

    byFrame_position.append(position)
    byFrame_order_parameters.append(byRes_order_parameters)
      
    # Don't need to write results for individual bond
    '''
    results = "%7i" % file_count
    for order_parameter in current_order_parameters:
      results += "%8.3f" % order_parameter
    stdout.write(" " + results + "\n")
    output.write(results + "\n")
    '''
    
    remove(filename)
    file_count += 1
  # End while loop

  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
  stdout.write("Snapshots analysis done.%s\n" % (" "*56))
  stdout.write("Computing averages...\n")

  # We don't reall need the time average of the order parameter
  '''
  # average order parameter
  averaged_order_parameters = []
  for i in range(len(bonds)):
    sum = 0.0
    for j in range(len(order_parameters)):
      sum += order_parameters[j][i]
    averaged_order_parameters.append(sum/len(order_parameters))
 
  # write results
  stdout.write("\n         " + bond_names)
  stdout.write(("-"*(len(output_legend) - 1)) + "\n")
  output.write(("-"*(len(output_legend) - 1)) + "\n")
  results = "average"
  for order_parameter in averaged_order_parameters:
    results += "%8.3f" % order_parameter
  stdout.write(" " + results + "\n")
  output.write(results + "\n")
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
  '''

  # Write abs average order parameters <Sn> (for carbon chains only)
  # WARNING this works with currenct lipids (all have defined x5 none carbon bonds) but for manually added lipids this might not be true

  # Here is the major modification A. Lorkowski made.
  # Instead of the absolute average order parameter of the time average,
  # the absolute average order parameter per time frame is calculated.

  output = open('order-'+str(lipid_type)+'.dat', 'w')
  format_txt = '{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('{LIPID}','{FRAME}','{RESNUM}','{X}','{Y}','{Z}','{ORDERPARAM}')

  output.write(format_txt)

  for i in range(len(byFrame_order_parameters)):
    #output = open('order-'+str(i)+'.dat', 'w')
    for j in range(number_of_lipids):
      ave_chain_s = 0
      for k in byFrame_order_parameters[i][j][hbCount:]:
        ave_chain_s += k
      average_txt = '{:>10} {:>10} {:>10} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f}\n'.format(lipid_type, i, j, byFrame_position[i][j][0], byFrame_position[i][j][1], byFrame_position[i][j][2], ( ave_chain_s / (len(byFrame_order_parameters[i][j])-hbCount)))
      #stdout.write(average_txt)
      output.write(average_txt)
    #output.close()
  stdout.write("Results written in \"order-"+str(lipid_type)+".dat\".\n")
  stdout.write("File in the following format: {LIPID} {FRAME} {RESNUM} {X} {Y} {Z} {ORDERPARAM}\n")
  output.close()

  # Below calculates the time average per residue.
  # The PBC effects must first be taken into account before this can be used.
  '''
  avg_output = open('avg_order-'+str(lipid_type)+'.dat', 'w')
  avg_format_txt = '{:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}\n'.format('{LIPID}','{-----}','{RESNUM}','{X}','{Y}','{Z}','{ORDERPARAM}')

  avg_output.write(avg_format_txt)

  for j in range(number_of_lipids):
    ave_chain_s = 0
    avg_x = 0
    avg_y = 0
    avg_z = 0
    for i in range(len(byFrame_order_parameters)):
      for k in byFrame_order_parameters[i][j][hbCount:]:
        ave_chain_s += k
      avg_x += byFrame_position[i][j][0]
      avg_y += byFrame_position[i][j][1]
      avg_z += byFrame_position[i][j][2]
    average_txt = '{:>10} {:>10} {:>10} {:>10.3f} {:>10.3f} {:>10.3f} {:>10.3f}\n'.format(lipid_type, '-----', j, avg_x/len(byFrame_order_parameters), avg_y/len(byFrame_order_parameters), avg_z/len(byFrame_order_parameters), ( ave_chain_s / ((len(byFrame_order_parameters[i][j])-hbCount)*len(byFrame_order_parameters))))
    avg_output.write(average_txt)
      
  stdout.write("Results written in \"avg_order-"+str(lipid_type)+".dat\".\n")
  stdout.write("File in the following format: {LIPID} {FRAME} {RESNUM} {X} {Y} {Z} {ORDERPARAM}\n")
  avg_output.close()
  '''
