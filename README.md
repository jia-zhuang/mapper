# mapper: 用于 NIPT 和 CNV 的短序列比对工具

## 0. 简介

mapper 是一个短序列比对工具，它可以高效地将固定长度的短序列比对到人基因组上。非常适合 NIPT/CNV 等低覆盖度二代测序产生的短序列的比对。它具有以下特点：

- 固定长度的短序列
- 精准匹配（不支持错配和 gap）
- 速度快 10M条/分钟，且并行效率不降低

从原理上讲，mapper 是基于索引的序列比对。在比对之前需要对人基因组进行索引，这个过程大约耗时 20~30 分钟，生成的索引文件大小约 38G。这个索引文件将被加载并常驻于共享内存中。之后的比对会直接检索存放于共享内存中的索引，速度非常快，并允许多个进程访问，彼此不影响，并行效率高。在实际应用中，对于一张测序芯片（96样本，每个样本数据量大约为 10M 对 PE reads），32核并行处理，比对时间在 5 分钟左右，因此，非常适合临床应用。

## 1. 系统要求

mapper 适用于服务器，不适合个人电脑。

- 硬盘大于 50G
- 内存大于 50G
- Linux 系统

## 2. 安装

```
git clone https://github.com/jia-zhuang/mapper.git
cd mapper
make
```

产生可执行文件 mapper，可以拷贝到任意位置使用。

## 3. 使用

- 建立索引

下载人参考基因组文件（hg19 或 hg38），预处理，只保留 1~22，X，Y，M 共 25 条染色体，且染色体命名规则为 chrN（N = 1~22，X，Y，M）。

```bash
mapper index hg19.fa
```

- 将索引加载到共享内存

```bash
mapper shm -l hg19.fa
```

可通过如下命令检出索引是否已存在于共享内存中：

```bash
mapper shm -c hg19.fa
```

或者打印共享内存里的文件名

```bash
mapper shm -p
```

销毁共享内存中的索引
```bash
mapper shm -d hg19.fa
```

- 比对

测序文件格式为 Fastq，支持压缩格式 .gz

```bash
mapper align hg19.fa read.fq.gz    # SE
mapper align hg19.fa read1.fq.gz read2.fq.gz    # PE
```

- 直接允许命令获取帮助

```bash
mapper
mapper index
mapper shm
mapper align
```

## 4. 注意事项和局限性

- 不适用于常规的序列比对，对于长序列、序列长度可变、允许错配和 gap的情况，请使用 [bwa](https://github.com/lh3/bwa)。
- 如果某个 reads 在基因组上存在多位点，该 reads 会被丢弃。
- 程序默认 reads 读长为 36bp，如果你的测序读长为其他值，编译前请修改 `utils.h` 文件中 `MER_LEN` 参数，但要确保是 4 的倍数。

