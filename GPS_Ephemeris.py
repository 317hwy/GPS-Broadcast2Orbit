ephemeris={
    #记录头（第1行前半）
    "PRN":None,
    "YEAR":None,
    "MONTH":None,
    "DAY":None,
    "HOUR":None,
    "MINUTE":None,
    "SECOND":None,

    #卫星钟参数（第1行后半）
    "af0":None,     #卫星钟偏差（s）
    "af1":None,     #钟漂移（s/s）
    "af2":None,     #钟漂移率（s/s^2）

    #轨道参数（第2行）
    "iode":None,        #星历数据期号
    "crs":None,         #轨道半径正弦调和改正项（m）
    "delta_n":None,     #平均角速度改正值（rad/s）
    "m0":None,          #参考时刻平近点角（rad）

    #轨道参数（第3行）
    "cuc":None,                 #纬度参数余弦调和改正项（rad）
    "eccentricity":None,        #轨道偏心率
    "cus":None,                 #纬度参数正弦调和改正项（rad）
    "sqrt_a":None,              #轨道半长轴的平方根（m^0.5）

    #轨道参数（第4行）
    "toe":None,         #参考时刻（s）
    "cic":None,         #升交点赤纬参数余弦调和改正项（rad）
    "omega0":None,      #升交点赤经（rad）
    "cis":None,         #升交点赤纬参数正弦调和改正项（rad）

    #轨道参数（第5行）
    "i0":None,          #轨道倾角（rad）
    "crc":None,         #轨道半径余弦调和改正项（m）
    "omega":None,       #近地点角距（rad）
    "omega_dot":None,   #升交点赤经率（rad/s）

    #状态参数（第6行）
    "idot":None,            #轨道倾角率（rad/s）
    "codes_on_l2":None,     #L2频道码标志
    "gps_week":None,        #GPS周数
    "l2_p_data_flag":None,  #L2频道数据标志

    #精度/健康/延迟（第7行）
    "sv_accuracy":None,     #卫星精度指标（URA）
    "sv_health":None,       #卫星健康状况
    "tgd":None,             #群时延差改正（s）
    "iodc":None,            #时钟数据期号

    #尾部参数（第8行）
    "transmission_time_of_message":None,    #消息发送时刻（s）
    "fit_interval":None                     #拟合区间（s）
}