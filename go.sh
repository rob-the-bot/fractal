#!/bin/bash
#SBATCH --job-name=qfract
#SBATCH --ntasks=3
#SBATCH --nodes=3
#SBATCH --cpus-per-task=16
#SBATCH --mem=2G
#SBATCH --time=2-0:00
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err
#SBATCH --partition=cosc

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

DATE=$(date +"%Y%m%d%H%M")
echo "time started  "$DATE

mkdir "${SLURM_ARRAY_TASK_ID}"
cd "${SLURM_ARRAY_TASK_ID}"

module load mpi
time mpirun ../fract ../"${SLURM_ARRAY_TASK_ID}"

(for f in ./*.png;
	do ffmpeg -loglevel panic -framerate 10 -i "$f" -vf "zoompan=z=zoom+0.1:x='iw/2-(iw/zoom/2)':y='ih/2-(ih/zoom/2)':d='10':s=1920x1080,scale=hd1080" -c:v libx264 -crf 16 -y "$f.mp4";
	done)

rm *.png

ffmpeg -loglevel panic -f concat -safe 0 -i <(for f in ./*.png.mp4; do echo "file '$PWD/$f'"; done) -c copy ../"${SLURM_ARRAY_TASK_ID}".mp4

rm *.png.mp4

DATE=$(date +"%Y%m%d%H%M")
echo "time finished "$DATE
