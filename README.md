# Gravitational Lens

## 简介

一个计算黑洞引力透镜效应的demo。

不解释，上图。

![Icon](https://github.com/yuchenxi2000/gravitational-lens/blob/master/pics/res6000-r800-d50-small.jpg)

## 光子轨迹具体是怎么算的？

### 旧版（糊，别看了）

见 数学推导.pdf

纠错（2019-11-26）：

1. 第5页开头theta_m少了一个常数1/w/sqrt(1+m)
2. 轨道常数D的计算过于复杂，其实D=1/sqrt(cc)。（第1页第5行方程直接取u=0，注意D=d{theta}/du=1/sqrt(cc)）
3. pdf里使用的图片采样方式是平铺，新版c++程序里使用天空盒。

### 讲稿

讲稿.pdf

## 我也想跑程序渲染图片。

src里c++代码自行编译。其中GNU Scientific Library需要自行下载、安装。

