mkdir sim_tet_0 && cd sim_tet_0 && ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py --coarse-steps 1e6 --fine-steps 1e6 --oxdna-steps 1e6 -d out --bp 11 --spacers 0 --nicks 1  ../plyfiles/tetrahedron.ply &
mkdir sim_cube_0 && cd sim_cube_0 && ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py --coarse-steps 1e6 --fine-steps 1e6 --oxdna-steps 1e6 -d out --bp 11 --spacers 0 --nicks 1  ../plyfiles/cube.ply &
mkdir sim_oct_0 && cd sim_oct_0 && ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py --coarse-steps 1e6 --fine-steps 1e6 --oxdna-steps 1e6 -d out --bp 11 --spacers 0 --nicks 1  ../plyfiles/octahedron.ply & 

mkdir sim_tet_4 && cd sim_tet_4 && ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py --coarse-steps 1e6 --fine-steps 1e6 --oxdna-steps 1e6 -d out --bp 11 --spacers 4 --nicks 1  ../plyfiles/tetrahedron.ply &
mkdir sim_cube_4 && cd sim_cube_4 && ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py --coarse-steps 1e6 --fine-steps 1e6 --oxdna-steps 1e6 -d out --bp 11 --spacers 4 --nicks 1  ../plyfiles/cube.ply &
mkdir sim_oct_4 && cd sim_oct_4 && ~/anaconda3/bin/python ../../../bin6/polymebabyonemoretime.py --coarse-steps 1e6 --fine-steps 1e6 --oxdna-steps 1e6 -d out --bp 11 --spacers 4 --nicks 1  ../plyfiles/octahedron.ply & 
