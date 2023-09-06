cd /home/zzhu/metilene_vscode

SECONDS=0;./metilene -s 1 -t 1000 ./metilene_DMR_simulation_scripts/bg/beta_40_3_1.2_0.8_ratio.1.metilene > ./metilene_output/N5x10_chr10_digits2.s1t1000c100.hpc.res;duration=$SECONDS;echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 1000 threads on HPC."

SECONDS=0;./metilene -s 1 -t 1000 ./metilene_DMR_simulation_scripts/bg/beta_40_3_1.2_0.8_ratio.0.87.metilene > ./metilene_output/N5x10_chr10_digits2.s1t1000c87.hpc.res;duration=$SECONDS;echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 1000 threads on HPC."

SECONDS=0;./metilene -s 1 -t 1000 ./metilene_DMR_simulation_scripts/bg/beta_40_3_1.2_0.8_ratio.0.73.metilene > ./metilene_output/N5x10_chr10_digits2.s1t1000c73.hpc.res;duration=$SECONDS;echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 1000 threads on HPC."

SECONDS=0;./metilene -s 1 -t 1000 ./metilene_DMR_simulation_scripts/bg/beta_40_3_1.2_0.8_ratio.0.6.metilene > ./metilene_output/N5x10_chr10_digits2.s1t1000c60.hpc.res;duration=$SECONDS;echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 1000 threads on HPC."
