 mxqsub --processors=1 --memory=1G --runtime=24h --tmpdir=1G --stdout ./out_t1.txt --stderr ./err_t1.txt ./mxq_N5x10_chr10_digits2_t1.sh
 mxqsub --processors=10 --memory=10G --runtime=1h --tmpdir=1G --stdout ./out_t10.txt --stderr ./err_t10.txt ./mxq_N5x10_chr10_digits2_t10.sh
 mxqsub --processors=100 --memory=100G --runtime=1h --tmpdir=1G --stdout ./out_t100.txt --stderr ./err_t100.txt ./mxq_N5x10_chr10_digits2_t100.sh
 mxqsub --processors=100 --memory=100G --runtime=1h --tmpdir=1G --stdout ./out_t1000.txt --stderr ./err_t1000.txt ./mxq_N5x10_chr10_digits2_t1000.sh
