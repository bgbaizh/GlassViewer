# GlassViewer

金属玻璃序参量计算包

Program Process：

    1，读取POSCAR
    2，分析POSCAR序参量
    3，返回数据和图像

    先形成相应模块
    而后集成到ofdfttest package 中实现job模式

Function:

    1.	Pair distribution function 	    Pyscal
    2.	Structure factor             	Undeveloped(Fourier to PDF)
    3.	Coordination number         	Pyscal
    4.	Chemical short-range order  	Undeveloped(Use CN)
    5.	Bond angle distribution      	Undeveloped
    6.	Bond orientational order    	Pyscal
    7.	Common neighbor analysis    	Pyscal
    8.	Voronoi tessellation        	Pyscal