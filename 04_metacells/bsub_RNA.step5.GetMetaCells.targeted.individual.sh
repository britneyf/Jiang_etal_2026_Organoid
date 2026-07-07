project=HTA.multiome_plasticity.combined.042023
log_dir=/data/peer/chanj3/$project/logs/RNA_step5.GetMetaCells.targeted.individual.100123
mkdir -p $log_dir

cd $log_dir


mdir=/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.individual.042023


for sample in `ls $mdir`; do
    num_bc=`wc -l $mdir/$sample/obs.$sample.042023.csv | awk '{print $1}'`

    for target_mc in 10 15 20 25 30; do 
	num_mc=$((num_bc/target_mc))

	odir=~/data/HTA.multiome_plasticity.combined.042023/out.RNA.individual.042023/$sample/SEACells/$target_mc
	mkdir -p $odir
	if `ls $odir | grep -q .$num_mc.h5ad`; then
	    :
	else    
	    bsub -J MC.$sample.$num_mc -n 10 -R span[ptile=10] -R rusage[mem=6] -W 24:00 -o MC.$sample.$num_mc.%J.stdout -eo MC.$sample.$num_mc.%J.stderr "conda activate multiome_env; /home/chanj3/anaconda3/envs/multiome_env/bin/python /home/chanj3/scripts/$project/GetMetaCells.individual.py $sample $num_mc"
	fi
    done
done

cd -
