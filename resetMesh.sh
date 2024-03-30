#! /bin/sh.

source $HOME/foam/foam-extend-5.0/etc/bashrc
rm -r 0.*
rm -r constant/polyMesh
rm CellCenterValues.txt.npy

python3 SimpleBlkMesh.py
blockMesh
python3 setFullCurrent.py
mv 0/Jcoil 0/JCoil
