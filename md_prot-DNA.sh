#!/bin/bash  -v
#SBATCH --partition=GPUSP4
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH -J ProtDNA					#job_name
#SBATCH --time=80:00:00
#SBATCH --gres=gpu:tesla:2

cd /scratch/cinthia/md_prot-dna

###GROMACS_2019###

module load Gromacs/2019.3-cuda
export OMP_NUM_THREADS=8

#Molecular Dynamics of protein-dna in water box

mkdir em npt nvt md

# Input = pdb file  
# First, preparing pdb file. For this, exclude waters molecules, from complex pdb.
# DON'T FORGET # to edit DNA pdb file in a text editor: Exclude waters and, for the 5'nucleotides, exclude the P, OP1, OP2 (5' fosfate group it's not necessary) according to dna.rtp of amber ff!! 

#topology of protein:
# Converting protein.pdb to gromacs input. Force field : Amber-99SB-ff with bsc1-ff modifications for nucleic acids, water type: spce
# hydrogen bonds in the pdb file are ignored (-ighn).

gmx pdb2gmx -f complex.pdb -o complex.gro -water spce -ff amber99sb_bsc1mod -ignh 

# Adding simulation box (cubic, min dist from wall 1.0 nm)
gmx editconf -f complex.gro -o complex_newbox.gro -c -d 1.4 -bt cubic

# Adding water mols (water type: spc216)
gmx solvate -cp complex_newbox.gro -cs spc216.gro -o complex_solv.gro -p topol.top

# Adding ions
# Add couter ions to the system (any mdp can be used for this step, the goal is only update the topology file);
gmx grompp -f ions.mdp -c complex_solv.gro -p topol.top -o ions.tpr

gmx editconf -f complex_solv.gro -o complex_waterbox.pdb

# Select SOL (you do not want to replace part of your molecule with ions)
echo "14" | gmx genion -s ions.tpr -o complex_solv_ions.gro -p topol.top -conc 0.15 -pname K -nname CL -neutral

# Energy minimization of protein + water (min.mdp -> boundary = xyz)

cd em

gmx grompp -f ../em_prot-dna.mdp -c ../complex_solv_ions.gro -p ../topol.top -o em_prot-dna.tpr

gmx mdrun -v -deffnm em_prot-dna -pin on -ntmpi 2 -ntomp 8

echo "10" "0" | gmx energy -f em_prot-dna.edr -o potential.xvg

# NVT Temperature control:
# Proper control of temperature coupling is a sensitive issue.
# Coupling every "moleculetype" (At top.top file) to its own thermostatting group is a bad idea. 
# For instance, if you do the following in mdp file":
# "nvt.mdp> tc-grps = Protein DNA SOL K", Your system will probably blow up, since the temperature coupling
# algorithms are not stable enough to control the fluctuations in kinetic energy that groups with a few atoms 
# (i.e., DNA and NA) will produce. Do not couple every single species in your system separately!!! NEVER!!
# In this case, is better to use in nvt file npt option tc-grps: DNA_Protein Water_and_Ions (Always uses only two groups for couple)!!!
# We need to create an DNA_Protein to these groups before use Grompp, then, type: 
# (1: DNA........... 3: Protein.........q: exit) with the command below:

echo -e "1 | 12 | \n q" | gmx make_ndx -f em_prot-dna.gro -o ../DNA_Protein.ndx

cd ../nvt

gmx grompp -f ../nvt_prot-dna.mdp -c ../em/em_prot-dna.gro -r ../em/em_prot-dna.gro -p ../topol.top -o nvt_prot-dna.tpr -n ../DNA_Protein.ndx

gmx mdrun -v -deffnm nvt_prot-dna -pin on -ntmpi 2 -ntomp 8

echo "15" "0" | gmx energy -f nvt_prot-dna.edr -o temperature.xvg

# NPT Pressure control

cd ../npt

gmx grompp -f ../npt_prot-dna.mdp -c ../nvt/nvt_prot-dna.gro -r ../nvt/nvt_prot-dna.gro -t ../nvt/nvt_prot-dna.cpt -p ../topol.top -o npt_prot-dna.tpr -n ../DNA_Protein.ndx

gmx mdrun -v -deffnm npt_prot-dna -pin on -ntmpi 2 -ntomp 8

echo "18" "0" | gmx energy -f npt_prot-dna.edr -o pressure.xvg

echo "24" "0" | gmx energy -f npt_prot-dna.edr -o density.xvg

# Real MD 

cd ../md
 
gmx grompp -f ../md_prot-dna.mdp -c ../npt/npt_prot-dna.gro -t ../npt/npt_prot-dna.cpt -p ../topol.top -o md_prot-dna.tpr -n ../DNA_Protein.ndx

gmx mdrun -v -deffnm md_prot-dna -pin on -ntmpi 2 -ntomp 8

#Trajectory Analysis

#Check trajectory

cd md
gmx check -f md_prot-dna.xtc

# Analysis generation of RMSD files and Gyration files
# This command line below convert trr in xtc to extract pdb. The option -pbc corrects the trajectory file.

echo "1" "0" |  gmx trjconv -s md_prot-dna.tpr -f md_prot-dna.xtc -o md_prot-dna_center.xtc -pbc mol -center -tu ns -ur compact

#Smoother trajectory visualization (to make a movie after)

echo "1" "0" | gmx trjconv -s md_prot-dna.tpr -f md_prot-dna_center.xtc -o md_prot-dna_fit.xtc -fit rot+trans

#Save a pdb file with pdb coordenates saved at every 10000 ps (10ns) (Multiple coordenates in a pdb file)
#dt = ps

echo "22" "0"  | gmx trjconv -f md_prot-dna.xtc -s md_prot-dna.tpr -n ../DNA_Protein.ndx -o md_prot-dna_.pdb -pbc nojump -dt 10000 -sep

# Protein coordinates Trajectory Extraction each 50 ps

echo "22" "0" | gmx trjconv -s md_prot-dna.tpr -f md_prot-dna_fit.xtc -n ../DNA_Protein.ndx -o movie.pdb -pbc nojump -dt 50

#Analysis generation of RMSD files and Gyration files
# Analyze is whether the protein is stable and close to the experimental structure by RMSD (em.tpr files vs md_prot-dna_center.tpr)
# RMSD analysis of protein backbone

echo "4" "4" | gmx rms -s ../em/em_prot-dna.tpr -f md_prot-dna_center.xtc -o rmsd_PROTbound_crystal.xvg -tu ns

echo "4" "4" | gmx rms -s md_prot-dna.tpr -f md_prot-dna_center.xtc -o rmsd_PROTbound_equilibrated.xvg -tu ns

# RMSD analysis of the DNA

echo "12" "12" | gmx rms -s ../em/em_prot-dna.tpr -f md_prot-dna_center.xtc -tu ns -o rmsd_DNAbound_crystal.xvg

echo "12" "12" | gmx rms -s md_prot-dna.tpr -f md_prot-dna_center.xtc -tu ns -o rmsd_DNAbound_equilibrated.xvg

# Radii of gyration analysis of protein

echo "1" | gmx gyrate -s md_prot-dna.tpr -f md_prot-dna_center.xtc -o gyrate_PROTbound.xvg

# Raddi of gyration analysis of DNA

echo "12" | gmx gyrate -s md_prot-dna.tpr -f md_prot-dna_center.xtc -o gyrate_DNAbound.xvg

# RMSF analysis of protein aa.

echo "4" | gmx rmsf -s md_prot-dna.tpr -f md_prot-dna_center.xtc -o rmsf_PROTbound.xvg -oq rmsf_bfac_PROTbound.pdb -od rmsdev_PROTbound.xvg -oc correl_PROTbound.xvg -res

# RMSF analysis of DNA

echo "12" | gmx rmsf -s md_prot-dna.tpr -f md_prot-dna_center.xtc -o rmsf_DNAbound.xvg -oq rmsf_bfacDNA_bound.pdb -od rmsdev_DNAbound.xvg -oc correl_DNAbound.xvg -res
