#!/usr/bin/env bash

SECONDS=0
# do some work
./metilene -s 1 beta.40.3_ratio.1.metilene > tmp.res
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 1 thread." >> elapsed.out


SECONDS=0
# do some work
./metilene -s 1 -t 10 beta.40.3_ratio.1.metilene > tmp.res
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 10 threads." >> elapsed.out


SECONDS=0
# do some work
./metilene -s 1 -t 128 beta.40.3_ratio.1.metilene > tmp.res
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 128 threads." >> elapsed.out


SECONDS=0
# do some work
./metilene -s 1 -t 1000 beta.40.3_ratio.1.metilene > tmp.res
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed, using 1000 threads." >> elapsed.out


rm tmp.res
