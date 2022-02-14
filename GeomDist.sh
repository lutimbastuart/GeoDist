#!/usr/bin/bash 
clear
echo -e "\e[34m \e[25m"
echo -e " \t                  Intertromics AB / AfriBioMol AB \n               Protein Structure Modelling / Drug Discovery Department \n        \t                  By Lutimba Stuart"

###-----------Documentation--------------------
: <<ScriptDocument_README
The name of the sample molecule is: 3-(2,4-Dichlorophenyl)-4-(1-methyl-1H-indol-3-yl)-1H-pyrrole-2,5-dione
Common Name: GSK-3 Inhibitor IV

Description;
This Script is written to:
1. Calculate and print the distance between atoms like C1 and C3 of the provide molecule 
2. Calculate and print the geometric center of the molecule.
3. Compering the results with the existing molecular Dynamics tools (VMD). To run with VMD uncomment the block please!
How to run the Script
    A. Locate / Path the script and the provided PDB file in one folder
    B. bash  GeomDist.sh 
ScriptDocument_README
echo -e "\e[97m"
##---------------END Documentation-------------
####-------Load the molecule for analysis------------##
echo -e "\tStructure File in .pdb format will be required:\n"
echo -e "\e[1m1. Calculating the Distance Between Atoms:\n"
echo -e "\e[0m \e[90m" 

###-----------User defined input --------------
read -e -p "Please in put a the PDB file_name of the Structure: 'Default is INH.pdb':   " PDB_file
read -e -p "Please select the resid number of the First Atom:   " atm_1
read -e -p "Please select the resid number of the Second Atom:   " atm_2

echo -e "Printing the Selected atoms\n"
echo "First Atom Selected" && echo | awk -v atmsel1="$atm_1" '{if ($2==atmsel1) print echo ($3 "\t" $4 "\t" "Chain  " $5)}' $PDB_file

echo "Second Atom Selected" && echo | awk -v atmsel2="$atm_2" '{if ($2==atmsel2) print echo ($3 "\t" $4 "\t" "Chain  " $5)}' $PDB_file

echo -e "\n"

#From First Principles 
#dist = √((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2) #dist=represents distance and (x,y,z) are Atomic Coordinates in 3D space
echo -e "\e[34m \e[3mExtracting Atomic Coordinates from the PDB file INH.pdb (GSK-3 Inhibitor IV MOlecule): ........."

#Atomic Coordinates of C1 atom in (x,y,z) dimension 
x_1=`echo | awk -v atmsel1="$atm_1" '{if ($2==atmsel1) print $7}' $PDB_file`
y_1=`echo | awk -v atmsel1="$atm_1" '{if ($2==atmsel1) print $8}' $PDB_file`
z_1=`echo | awk -v atmsel1="$atm_1" '{if ($2==atmsel1) print $9}' $PDB_file`

echo -e "\e[0mAtomic Coordinates of Atom $atm_1 (X Y Z):  \n $x_1\t$y_1 \t$z_1"

#Atomic Coordinates of C3 atom in (x,y,z) dimension 
x_2=`echo | awk -v atmsel2="$atm_2" '{if ($2==atmsel2) print $7}' $PDB_file`
y_2=`echo | awk -v atmsel2="$atm_2" '{if ($2==atmsel2) print $8}' $PDB_file`
z_2=`echo | awk -v atmsel2="$atm_2" '{if ($2==atmsel2) print $9}' $PDB_file`

echo -e "\nAtomic Coordinates of Atom $atm_2 (X Y Z):  \n $x_2 \t$y_2 \t$z_2"

echo -e "\n\e[1m \e[34m \e[3mCalculating the distance between atoms exploiting the above atomic coordinate information: ..."
echo -e "Implimenting the fomular Dist = √((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)"
echo -e "\e[0m"

Dist=$(echo "($x_2 - $x_1)^2 + ($y_2 - $y_1)^2 + ($z_2 - $z_1)^2" | bc -l) 
Dist_btn_atoms=`echo "scale=10; sqrt($Dist)" | bc -l`

echo -e "\e[1m \e[3mAns.1 The distance between atom $atm_1 and atom $atm_2 is \e[104m$Dist_btn_atoms Å\e[0m"

echo -e "\n\e[1m2.  Calculate and print the geometric center of the molecule | $PDB_file"

########################################################################################
##To obtain the Geometric Center of the molecule:                                      #
#1. Get the Coordinates of the X Y Z axis                                              #
#2. Calculate the average/mean of each axis                                            # 
#3. Report and Print the Geometric center as the obtained avearge in X, and Y dimension#
########################################################################################

echo -e "\n\e[34m \e[3mExtracting Atomic Coordinate Columns in each dimension from the PDB file | $PDB_file:......... " 

#Selecting Coordinates in the respective Dimension
X_axis=`echo | awk '{print $6}' $PDB_file | wc -w`  # In the X dimension 
Y_axis=`echo | awk '{print $7}' $PDB_file | wc -w`  # In the Y dimension 
Z_axis=`echo | awk '{print $8}' $PDB_file | wc -w`  # In the Z dimension 

# Getting the average/Mean  of each axis (X, Y and Z)   
x_sum=`echo | awk '{counts+=$6}; END{print counts}' $PDB_file`
y_sum=`echo | awk '{counts+=$7}; END{print counts}' $PDB_file`
z_sum=`echo | awk '{counts+=$8}; END{print counts}' $PDB_file`

echo -e "\n\e[95m \e[3mAns.2  \t         Geometric Center of the molecule (Å)" 
echo -e "\e[0m"
x_geocent=`echo "scale=5; ($x_sum/$X_axis)" | bc -l`
y_geocent=`echo "scale=5; ($y_sum/$Y_axis)" | bc -l`
z_geocent=`echo "scale=5; ($z_sum/$Z_axis)" | bc -l`

echo -e "\e[1m \e[3m" 
echo -e  "  X     \t\t    Y                      Z"
printf "%f\t\t" "$x_geocent" "$y_geocent"  "$z_geocent" 

echo -e "\n\n................................END.............................."
echo -e "\e[2m\t\t           Intetromics AB"

#To compare using the existing molecular dynamics analysis tools (Visual Molecular Dynamics -VMD)
<< 'Compare_Results'
vmd -dispdev text <<START
#Geting the distance between atoms 
mol new \$PDB_file
set atom_1 [atomselect top "name C1"]
set atom_2 [atomselect top "name C3"]

measure rmsd \$atom_2 \$atom_1

#Geting the Geometry center
set mol_all [atomselect top all]
measure center \$mol_all
START
Compare_Results
