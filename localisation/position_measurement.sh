#!/bin/bash -l
#SBATCH --account=pawsey0272
#SBATCH --partition=work
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=8
#SBATCH --mem=50G
#SBATCH --nodes=1

source /software/projects/pawsey0272/nhurleywalker/.login

#export SINGULARITY_BINDPATH=

set -x

image=J1755-25_1729341386_scan18_deeper-MFS-image.fits
crop=${image%.fits}_crop.fits
singularity exec -B $PWD ${GXCONTAINER} getfits -o $crop $image 17:55:34.9 -25:27:49.1 500 500

singularity exec -B $PWD ${GXCONTAINER} BANE --cores=1 ${crop}
singularity exec -B $PWD ${GXCONTAINER} aegean --autoload --seedclip=6 --floodclip=3 --table=${crop} ${crop}
