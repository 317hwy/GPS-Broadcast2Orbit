from GPS_Ephemeris import ephemeris
from copy import deepcopy
from itertools import islice

def _to_float(s:str):
    s=s.strip().replace("D","E").replace("d","E")
    return float(s) if s else None

def _to_int(s:str):
    s=s.strip().replace("D","E").replace("d","E")
    return int(float(s)) if s else None

#按行从广播星历的"END OF HEADER"后开始读取文件，得到仅含有星历的文件
def save_nav_body_file(input_path,output_path):
    with open(input_path,'r',encoding='utf-8',errors='ignore') as fin , open(output_path,'w',encoding='utf-8',errors='ignore') as fout:
        is_header=False
        for line in fin:
            if "END OF HEADER" in line:
                is_header=True
            elif is_header:
                fout.write(line)
    if not is_header:
        raise ValueError("未找到'END OF HEADER'，请检查输入文件格式")

#从仅含有星历的文件中读取数据，返回一个包含所有卫星星历的列表
def read_nav_body_file(input_path):
    with open(input_path,'r',encoding='utf-8',errors='ignore') as f:
        #建立一个空列表来存储所有卫星的星历数据
        ephemeris_list=[]
        #每次读取八行，分别对应一个卫星的星历数据
        while True:
            eph=deepcopy(ephemeris)   #每次读取前先复制一个空的星历字典
            block=list(islice(f,8))     #一次读八
            if not block:
                break
            if len(block)<8:
                raise ValueError(f"星历块不完整：期望8行，实际{len(block)}行")
            block=[line.rstrip("\n") for line in block] #去掉每行末尾的换行符
            #第1行
            eph['PRN']=block[0][0:3].strip()
            eph['YEAR']=int(block[0][4:8].strip())
            eph['MONTH']=int(block[0][9:11].strip())
            eph['DAY']=int(block[0][12:14].strip())
            eph['HOUR']=int(block[0][15:17].strip())
            eph['MINUTE']=int(block[0][18:20].strip())
            eph['SECOND']=int(block[0][21:23].strip())

            eph['af0']=_to_float(block[0][23:42])
            eph['af1']=_to_float(block[0][42:61])
            eph['af2']=_to_float(block[0][61:80])
            #第2行
            eph['iode']=_to_float(block[1][4:23])
            eph['crs']=_to_float(block[1][23:42])     
            eph['delta_n']=_to_float(block[1][42:61])
            eph['m0']=_to_float(block[1][61:80])
            #第3行
            eph['cuc']=_to_float(block[2][4:23])
            eph['eccentricity']=_to_float(block[2][23:42])     
            eph['cus']=_to_float(block[2][42:61])
            eph['sqrt_a']=_to_float(block[2][61:80])
            #第4行
            eph['toe']=_to_float(block[3][4:23])
            eph['cic']=_to_float(block[3][23:42])
            eph['omega0']=_to_float(block[3][42:61])
            eph['cis']=_to_float(block[3][61:80])
            #第5行
            eph['i0']=_to_float(block[4][4:23])
            eph['crc']=_to_float(block[4][23:42]) 
            eph['omega']=_to_float(block[4][42:61])
            eph['omega_dot']=_to_float(block[4][61:80])
            #第6行
            eph['idot']=_to_float(block[5][4:23])
            eph['codes_on_l2']=_to_float(block[5][23:42])
            eph['gps_week']=_to_int(block[5][42:61])
            eph['l2_p_data_flag']=_to_float(block[5][61:80])
            #第7行
            eph['sv_accuracy']=_to_float(block[6][4:23])
            eph['sv_health']=_to_float(block[6][23:42])
            eph['tgd']=_to_float(block[6][42:61])
            eph['iodc']=_to_float(block[6][61:80])
            #第8行
            eph['transmission_time_of_message']=_to_float(block[7][4:23])
            eph['fit_interval']=_to_float(block[7][23:42])
            ephemeris_list.append(eph)
    return ephemeris_list