#! /bin/bash
#
# Based on the run_wrapper for JUWELS Booster by Luis Kornblueh
# Adapted to Levante (only 2 NICs instead of 4 as on JUWELS Booster)
#
# nvidia-smi topo -mp
#       GPU0    GPU1 GPU2      GPU3   mlx5_0      mlx5_1 CPU Affinity        NUMA Affinity
# GPU0   X     SYS  SYS        SYS      SYS    SYS  48-63,176-191      3
# GPU1  SYS         X   SYS        SYS      PIX    SYS  16-31,144-159      1
# GPU2  SYS        SYS       X         SYS      SYS    PIX  112-127,240-255    7
# GPU3        SYS      SYS    SYS   X     SYS  SYS        80-95,208-223    5
# mlx5_0      SYS    PIX  SYS        SYS       X         SYS              
# mlx5_1 SYS       SYS     PIX   SYS SYS        X          
#
# Legend:
#
#  X    = Self
#  SYS  = Connection traversing PCIe as well as the SMP interconnect between NUMA nodes (e.g., QPI/UPI)
#  NODE = Connection traversing PCIe as well as the interconnect between PCIe Host Bridges within a NUMA node
#  PHB  = Connection traversing PCIe as well as a PCIe Host Bridge (typically the CPU)
#  PXB  = Connection traversing multiple PCIe bridges (without traversing the PCIe Host Bridge)
#  PIX  = Connection traversing at most a single PCIe bridge
#
#____________________________________________________________________________________________________


io_tasks=0

while getopts n:o:e: argv
do
    case "${argv}" in
        n) mpi_total_procs=${OPTARG};;
        o) io_tasks=${OPTARG};;
        e) executable=${OPTARG};;
    esac
done

set -eu

lrank=$SLURM_LOCALID%4

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^vader,tcp,openib,smcuda"

# need to check in run script that the variables make sense and are
# exported!

(( compute_tasks = mpi_total_procs - io_tasks ))

if (( SLURM_PROCID < compute_tasks ))
then

    echo Compute process $SLURM_LOCALID on $(hostname)

    numanode=(2-3 0-1 6-7 4-5)
    gpus=(0 1 2 3)
    nics=(mlx5_0:1 mlx5_0:1 mlx5_1:1 mlx5_1:1)
    reorder=(0 1 2 3)

    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]})
    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]})

    export UCX_NET_DEVICES=${nic_reorder[lrank]}
    export CUDA_VISIBLE_DEVICES=${gpus[${reorder[lrank]}]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=16384

    export UCX_IB_GPU_DIRECT_RDMA=yes

    export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy
    export UCX_MEMTYPE_CACHE=n

    export JULIA_CUDA_MEMORY_POOL=none

    #export UCX_WARN_UNUSED_ENV_VARS=n

    #export UCX_ERROR_SIGNALS="SIGILL,SIGSEGV,SIGBUS,SIGFPE"

    #export OMPI_MCA_pml=ucx
    #export OMPI_MCA_osc=ucx

else

    echo IO process $SLURM_LOCALID on $(hostname)

    numanode=(2-3 0-1 6-7 4-5)
    nics=(mlx5_0:1 mlx5_0:1 mlx5_1:1 mlx5_1:1)    
    reorder=(0 1 2 3)
 
    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]}) 

    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]}) 
    
    export UCX_NET_DEVICES=${nic_reorder[lrank]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=16384

    export UCX_IB_GPU_DIRECT_RDMA=yes

    export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy    
    export UCX_MEMTYPE_CACHE=n

    export JULIA_CUDA_MEMORY_POOL=none

    export UCX_WARN_UNUSED_ENV_VARS=n

fi

numactl --cpunodebind=${numanode_reorder[$lrank]} --membind=${numanode_reorder[$lrank]} $executable


