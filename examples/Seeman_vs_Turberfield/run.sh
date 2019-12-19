cd Seeman && 
staplehedra ../../../ply_files/tetrahedron.ply --spacers 0 --bp 21 --coarse-steps 1e5 --fine-steps 1e5 --oxdna-steps 1e6 -d out --nicks 1 &  
cd ..

cd Turberfield && 
staplehedra ../../../ply_files/tetrahedron.ply --spacers 0 --bp 21 --coarse-steps 1e5 --fine-steps 1e5 --oxdna-steps 1e6 -d out --nicks 0 --ldmoverhangfile ../../../overhang_files/ldm_tet_seq.overhang & 

