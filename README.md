# swSpTRSV
Sparse triangular solve for Sunway architecture

1.Login to Sunway TaihuLight

- Land the homepage of Sunway TaihuLight Supercomputer: http://www.nsccwx.cn/wxcyw/

- Select one VPN service: 'Telecom', 'Unicom' or 'China Mobile' on the top of the website. Please choose the best one for better connection.

- Login to Sunway TaihuLight Supercomputer: ssh 41.0.0.188

Please email to clarencewxl@gmail.com for further questions.

2.Experiments for the Temperature Map

- Build the binary:

cd ./SpTRSV_map

makeall.sh

- Run a single benchmark

./test.sh -I benchmark_name

Note: The benchmarks' names are listed in the file: benchmarks\_list in the same directory

- Use nohup to run all the benchmarks

nohup ./total_run.sh &

Note: Use bjobs to check whether the jobs finish or not. Each job will product a new file in the subdirectory 'results' in the current directory.

- The referenced performance results are listed in the file 'ref_map' in the current directory.

Note: There are totally 4x6x2057=49368 works and will cost about 36 hours. As works can be currently done using different processors, depending on the task-schedule system, you do not need to wait for finishing all the works and can compare any work's result with the referenced result as long as it is finished.

3.Experiments for the Scatter Plot

- Build the binary:

cd ./SpTRSV_tuning

make

- Run a single benchmark

./run.sh -I benchmark_name

Note: The benchmarks' names are listed in the file: benchmarks\_list in the same directory

- Use nohup to run all the benchmarks

nohup ./total_run.sh &

Note: Use bjobs to check whether the jobs finish or not. Each job will product a new file in the subdirectory 'results' in the current directory.

- The referenced performance results are listed in the file 'ref_SwSpTRSV', 'ref_serial' and 'ref_levelsets' in the current directory. 

Note: There are totally $2057$ works and will cost about 2 hours. As works can be currently done using different processors, depending on the task-schedule system, you do not need to wait for finishing all the works and can compare any work's result with the referenced result as long as it is finished.
