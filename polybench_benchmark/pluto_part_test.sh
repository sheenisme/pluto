#!/bin/bash
workdir=$(cd `dirname $0`; pwd)
# echo "当前的工作路径是:" $workdir
cd $workdir



# 配置PLC路径
PLC=/media/sheen/文件/ubuntu_code_data/pluto/polycc


all_benchs=$(cat ./utilities/benchmark_list)
for bench in $all_benchs;
do
    benchdir=$(dirname $bench)
    benchname=$(basename $benchdir)
    # echo $benchdir " " $benchname
    # 进入测试用例的目录进行测试
    cd $benchdir
    echo -e "32\n32\n32\n32" > tile.sizes
    $PLC ${benchname}.c --tile --noprevector -o ${benchname}.pluto.c

    # 返回测试脚本目录
    cd $workdir
done



# 返回测试脚本目录
cd $workdir
# 准备
rm -rf build
rm -rf result-out


./taffo_collect-fe-stats.sh pluto.vs.origion_0305

echo "pluto vs origion, over!"
