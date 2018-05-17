倾角道集叠前时间偏移备份说明

CreateDangle.c 由叠前共偏移距数据输入道常孔径叠前时间偏移产生倾角道集
XCreateDangle CreateDangle.c的可执行壳文件

jswxdcdpd.c 对CreateDangle.c生成的倾角道集进行单CDP扫描

jswxdcdpq.c 对CreateDangle.c生成的倾角道集进行全CDP扫描

DangleSmALL.c 由叠前共偏移距数据输入道常孔径叠前时间偏移产生倾角道集，并进行倾角扫描
XDangleSmALL DangleSmALL.c的可执行壳文件

srdmigD.c 采用jswxdcdpq.c或DangleSmALL.c生成的倾角信息，进行输入道单偏移距偏移
XsrdmigD srdmigD.c的可执行壳文件

srdmigA.c 采用jswxdcdpq.c或DangleSmALL.c生成的倾角信息，进行输入道全偏移距偏移
XsrdmigA srdmigA.c的可执行壳文件

srdmigAcmp.c 采用jswxdcdpq.c或DangleSmALL.c生成的倾角信息，进行输入道全偏移距偏移，并输出共crp道集
XsrdmigAcmp srdmigAcmp.c的可执行壳文件

scdmigD.c 采用jswxdcdpq.c或DangleSmALL.c生成的倾角信息，进行输出道单偏移距偏移
XscdmigD scdmigD.c的可执行壳文件

scdmigA.c 采用jswxdcdpq.c或DangleSmALL.c生成的倾角信息，进行输出道全偏移距偏移
XscdmigA scdmigA.c的可执行壳文件

srdrecreateDAngle.c 采用jswxdcdpq.c或DangleSmALL.c生成的倾角信息，进行输入道全偏移距偏移，并再次生成倾角道集，验证算法
XsrdrecreateDAngle srdrecreateDAngle.c的可执行壳文件
