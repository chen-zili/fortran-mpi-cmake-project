# fortran-mpi-cmake-project

为演示在Fortran中使用MPI，该仓库提供下面四个简单的例子：

1. 阻塞点对点通信
2. 非阻塞点对点通信
3. 非阻塞集体通信
4. 单边通信

## TAU

为了分析代码性能，这里简单演示了TAU的使用。

### 下载

首先从[TAU官方网站](https://www.cs.uoregon.edu/research/tau/downloads.php)下载安装所需PDT文件以及TAU目录

### PDT安装

随后需要安装PDT文件，解压后在目录中运行：
```shell
./configure -ICPC -prefix=/home/chenzili/pdt
make install -j8
```

### TAU安装

接下来配置、编译和安装 TAU (需要先切换到Python 2环境)
```shell
./configure -c++=mpiicpc -cc=mpiicc -fortran=mpiifort -mpi -ompt -iowrapper -bfd=download -dwarf=download -otf=download -unwind=download -pdt=/home/chenzili/pdt -prefix=/home/chenzili/tau
make install -j8
```

详细过程见 PDT 和 TAU 目录中的 Readme 和 INSTALL

### TAU使用

这里仅仅使用基于命令行的tau_exec工具

```shell
mpirun -np 4 tau_exec -T mpi,openmp,ompt,pdt -ompt ./test_01
pprof
```

## 参考链接

- [Message Passing Interface (MPI)] (https://hpc-tutorials.llnl.gov/mpi/)
- [TAU: Tuning and Analysis Utilities](https://hpc.llnl.gov/software/development-environment-software/tau-tuning-and-analysis-utilities)
- [TAU Installation](https://lsi2.ugr.es/jmantas/ppr/ayuda/datos/instalaciones/Install_TAU_en.pdf)
- [INTERMEDIATE MPI](https://enccs.github.io/intermediate-mpi/#)
- [ADVANCED MPI](https://hpc.llnl.gov/sites/default/files/DavidCronkSlides.pdf)
