# Gravitational Lens

## 简介

一个计算黑洞引力透镜效应的demo。

不解释，上图。

![Icon](https://github.com/yuchenxi2000/gravitational-lens/blob/master/pics/photon88.bmp)

![Icon](https://github.com/yuchenxi2000/gravitational-lens/blob/master/pics/photon0.bmp)

## 光子轨迹具体是怎么算的？

见"数学推导.pdf"。

自己算的。如果有错欢迎指正。

## 我也想跑程序渲染图片。

你需要如下library:

1. glew（OpenGL相关的库）
2. glfw（OpenGL相关的库）
3. opengl（显卡要支持）
4. SOIL（图像加载、存储的库。比较简陋。）
5. GSL（GNU科学计算库）

我觉得123不需要，但是SOIL没有这些运行不起来。我用123成功编译，但可能有其他的解决办法。
