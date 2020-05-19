#mpiexec -n 4 --machinefile hosts.txt --map-by node ./proj3 6 53
time mpiexec --mca btl self --mca btl_openib_cpc_include rdmacm --machinefile ./hosts.txt -n 64 --map-by node ./LU_par 1000 1
# time mpiexec --machinefile ./hosts.txt -n 9 --map-by node ./proj3 9 10
