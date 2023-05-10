### introduction：

​	伴随汽车产业与物联网的深度融合，汽车的功能和设计趋于智能化，越来越多的雷达、传感器及人工智能模块加装在汽车上。其中，雷达可以帮助汽车获取周边物体的位置信息和速度信息，起到了不可或缺的作用。汽车防撞雷达不仅可以帮助驾驶员更清楚地了解周围障碍物的情况，提高驾驶的安全性，而且对周边人物的探测信息也可以作为自动驾驶领域的信息来源，是边缘计算的重要一环。汽车防撞雷达的设计需要考虑到雷达探测的距离范围、速度范围、测距精度、测速精度、距离分辨率等因素。

​	本项目是笔者《雷达原理与系统》课程设计的仿真部分，主要讨论了汽车防撞雷达中雷达波形的参数设计、计算和性能仿真。考虑到雷达计算量和精度的问题，我们选择了基于MSLFSK波形体制的雷达，可以利用该波形解算出目标耦合的距离和速度信息。雷达战术指标确定后，通过计算仿真，得到了我们需要构建的波形参数。最后通过测距测速仿真观察到目标的距离速度精度可以满足战术指标要求。

### explain：

![屏幕截图 2023-05-09 203901](.\pics\屏幕截图 2023-05-09 203901.jpg)

![屏幕截图 2023-05-09 203927](.\pics\屏幕截图 2023-05-09 203927.jpg)

​	根据汽车驾驶的一般情况，我们的雷达战术指标要求雷达最大可以探测到前方120m的目标，又根据探测目标例如：静物、行人、单车、汽车甚至火车的移动速度，确定雷达探测速度范围为-240~240km/h，测距和测速精度分别为0.5m和1.8km/h，考虑到目标距离和速度的可区分性，确定距离分辨率为1m，速度分辨率为1km/h。

​	计算方法和公式不再详述，直接run cal.m，根据战术指标得到调频带宽、相参处理间隔、LFSK信号个数、LFSK信号中频率键控次数、LFSK信号中频率键控时宽、LFSK信号中键控频差、不同LFSK信号频差 等参数。

​	接着利用MSLFSK波形对目标进行测距测速仿真。MSLFSK_noisefree.m中雷达信号未加噪，MSLFSK_noise.m中信号加噪，这里只分析MSLFSK_noise.m的处理过程：

· 设置目标距离和速度

· 产生发射信号MSLFSK波形时频图、收发时频图

· 生成发射信号、接收信号（信道加噪）和混频信号（接收机加噪）

· 将混频信号数据重排，分成M路LFSK，分别进行傅里叶变换，1通道LFSK信号频谱

· 进行频谱周期图积累后，M路LFSK频谱

· 提取出fmslfsk(p)后，计算M路LFSK频谱在p处的相位差，计算各模糊距离，并解模糊，得到目标测量结果

· 通过仿真测量结果，在信道噪声和接收机噪声的共同影响下，目标测量的误差仍保持在一个较低水平，可以满足我们上文设定的测距精度和测速精度。



### conclusion：

