date=020124
project=CRC_ZFP36L2.092023/Organoid/
log_dir=/data/chanjlab/$project/logs_new/metacells/metacells.$date/
mkdir -p $log_dir

cd $log_dir

in_dir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/input/metacells/

for sample in $(ls $in_dir); do
    num_bc=$(wc -l $in_dir/$sample/obs.$sample.csv | awk '{print $1}')

    for target_mc in 10 15 20 25 30; do 
        num_mc=$((num_bc/target_mc))

        out_dir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output_new/metacells/metacells.$date/$sample/$target_mc/
        mkdir -p $out_dir
        if [ ! -f $out_dir/metacells_adata.$sample.$target_mc.h5ad ]; then
            echo "Processing sample: $sample"
            bsub -J MC.$sample.$target_mc -n 10 -R span[ptile=10] -R rusage[mem=6] -W 24:00 -o MC.$sample.$target_mc.%J.stdout -eo MC.$sample.$target_mc.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; /home/forsythb/anaconda3/envs/scrna/bin/python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/9_metacells.py $sample $num_mc $target_mc"
        else
            echo "Output file already exists for sample: $sample, skipping."
        fi
    done
done

cd -


# project=CRC_ZFP36L2.092023/Organoid/
# log_dir=/data/chanjlab/$project/logs/metacells
# #log_dir=/data/peer/chanj3/$project/logs/RNA_step5.GetMetaCells.targeted.individual.100123
# mkdir -p $log_dir

# cd $log_dir

# mdir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/input/metacells/
# #mdir=/data/peer/chanj3/HTA.multiome_plasticity.combined.042023/out.RNA.individual.042023

# for sample in `ls $mdir`; do
#     num_bc=`wc -l $mdir/$sample/obs.$sample.csv | awk '{print $1}'`

#     for target_mc in 10 15 20 25 30; do 
# 	num_mc=$((num_bc/target_mc))

# 	odir=/data/chanjlab/CRC_ZFP36L2.092023/Organoid/output/metacells/$sample/SEACells/$target_mc
# 	mkdir -p $odir
# 	if `ls $odir | grep -q .$num_mc.h5ad`; then
# 	    :
# 	else    
# 	    bsub -J MC.$sample.$num_mc -n 10 -R span[ptile=10] -R rusage[mem=6] -W 24:00 -o MC.$sample.$num_mc.%J.stdout -eo MC.$sample.$num_mc.%J.stderr "source ~/anaconda3/etc/profile.d/conda.sh; conda activate scrna; /home/forsythb/anaconda3/envs/scrna/bin/python /data/chanjlab/forsythb/organoid_analysis_pipeline_scripts/9_metacells.py $sample $num_mc"
        
# 	fi
#     done
# done

cd -
