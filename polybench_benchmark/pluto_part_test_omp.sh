#!/bin/bash
workdir=$(cd `dirname $0`; pwd)
# echo "当前的工作路径是:" $workdir
cd $workdir



# 配置PLC路径
PLC=/home/sheen/pluto/polycc

# 设置编译脚本
sed -n "s/OMPSET=\" \"/OMPSET=\"-fopenmp\"/p" taffo_compiler.sh
sed -i "s/OMPSET=\" \"/OMPSET=\"-fopenmp\"/g" taffo_compiler.sh



array=( 2 4 8 )
for element in ${array[@]}
do
    # 先将run的脚本设置成对应的核数
    sed -n "s/OMP_NUM_THREADS=1/OMP_NUM_THREADS=${element}/p" taffo_run.sh
    sed -i "s/OMP_NUM_THREADS=1/OMP_NUM_THREADS=${element}/g" taffo_run.sh

    all_benchs=$(cat ./utilities/benchmark_list)
    for bench in $all_benchs;
    do
        benchdir=$(dirname $bench)
        benchname=$(basename $benchdir)
        # echo $benchdir " " $benchname
        # 进入测试用例的目录进行测试
        cd $benchdir
        # echo -e "32\n32\n32\n32" > tile.sizes
        $PLC ${benchname}.c --tile --noprevector --cloogsh --parallel -o ${benchname}.pluto.c

        # 返回测试脚本目录
        cd $workdir
    done

    # 返回测试脚本目录
    cd $workdir
    # 准备
    rm -rf build
    rm -rf result-out

    ./taffo_collect-fe-stats.sh pluto_test_result_${element}

    echo "pluto omp ${element}, over!"

    # 进行清理
    cd $workdir
    cd utilities/
    perl clean.pl ../
    perl makefile-gen.pl ../  -cfg

    # 返回测试脚本目录
    cd $workdir


    # 重置run的脚本
    sed -n "s/OMP_NUM_THREADS=${element}/OMP_NUM_THREADS=1/p" taffo_run.sh
    sed -i "s/OMP_NUM_THREADS=${element}/OMP_NUM_THREADS=1/g" taffo_run.sh
done



# 复原编译脚本
sed -n "s/OMPSET=\"-fopenmp\"/OMPSET=\" \"/p" taffo_compiler.sh
sed -i "s/OMPSET=\"-fopenmp\"/OMPSET=\" \"/g" taffo_compiler.sh