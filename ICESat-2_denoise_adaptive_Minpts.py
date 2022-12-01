import pandas as pd
import os
import math
import numpy as np

Radius = 2
PI = 3.141593
h2 = 50


#该函数用于寻找递增的数列list中值小于value_end的最大索引index_start
def index_find(list, value_end,index_start):
    for i in range(index_start,len(list)):
        if(i == len(list)-1):
            return i
        else:
            if(list[i]>value_end):
                return i-1


#该函数用于将递增的数列list按照步长step和切片长度slice_dis进行切片
#每个切片结果都以（index_start,index_end)的形式存储在slice_index数列中
def list_slice(list,step,slice_dis):
    slice_num = int((max(list)-min(list))/step)+1
    slice_index = [[0,0] for i in range(slice_num)]
    for i in range(slice_num):
        if(i==0):
            slice_index[i][0] = 0
            slice_index[i][1] = index_find(list,list[slice_index[i][0]]+slice_dis,slice_index[i][0])
        else:
            slice_index[i][0] = index_find(list,list[slice_index[i-1][0]]+step,slice_index[i-1][0])
            slice_index[i][1] = index_find(list, list[slice_index[i][0]] + slice_dis, slice_index[i][0])
            if(slice_index[i][0] == slice_index[i][1] and slice_index[i][0]!=len(list)-1):
                slice_index[i][0] +=1
                # slice_index[i][1] +=1
    return slice_index


#该函数用于统计递增的数列list中值最顶端和最低端长度为h2部分数据的数量
#因为光子数据中最顶部和最底部的数据一般都是噪声数据
def list_statis(list):
    num = 0
    for i in list:
        if (i < (min(list) + h2) or i>(max(list) - h2)):
            num += 1
    return num


#该函数用于根据邻域点数量区分信号点和噪声点
def noi_sig_class(along_track, height, MinPts, Radius, clss_res):
    for i in range(len(along_track)):
        num = 0
        for j in range(len(along_track)):
            if(j==1):
                continue
            #首先初筛一遍，以Radius为半边长，以(along_track[i],height[i])为中心构建矩形
            #只对矩形内的点进行距离的计算
            if (along_track[j] <= along_track[i] + Radius and along_track[j] >= along_track[i] - Radius and height[j] <=
                    height[i] + Radius and height[j] >= height[i] - Radius):
                #计算i点(along_track[i],height[i])和j点(along_track[j],height[j])的欧氏距离
                dis = math.hypot(along_track[i] - along_track[j], height[i] - height[j])
                #如果欧式距离小于设定的阈值Radius，则将j点视为i点的邻域点
                if (dis < Radius):
                    num += 1
        #如果i点的邻域点数量num比设定的阈值MinPts小，则将i点视为噪声点，反之视为信号点
        if (num < MinPts):
            clss_res.append(0)
        else:
            clss_res.append(1)


def DRAGANN(input_csv, outputpath):
    #读取csv中的属性
    formmer_path, csv_filename = os.path.split(input_csv)
    df = pd.read_csv(input_csv)
    # df = df.drop(df[df['Signal Confidence'] <= 0].index)
    df = df.sort_values(by='Along-Track (m)', ascending=True)

    lat = list(df['Latitude (deg)'])
    lon = list(df['Longitude (deg)'])
    along_track = list(df['Along-Track (m)'])
    along_track = [i/50 for i in along_track]
    h = list(df['Height (m MSL)'])
    
    #如果along_track里没有数据，就输出csv的名字，并返回0
    try:
        along_start = along_track[0]
    except:
        print(input_csv)
        return 0

    #初始化一些变量
    dis = 3#分段的步长
    clss_res = []
    cou = 0#用于记录迭代的次数
    along_track_tmp = []#用于记录每段数据中的沿轨道长度
    h_tmp = []#用于记录每段数据中的高程
    
    #下面的循环主要是对along_track进行分段，分段的距离为dis（注意这个是值的距离，而不是索引的距离)
    #然后对每段进行信噪区分
    #分段的意义是为了消除不同地方信噪的阈值不同的影响
    for i in range(len(along_track)):
    
        #如果along_track[i]到along_start的距离小于dis，且i点不是最后一个点
        #表明此时还没有遍历到目前这段数据的尽头
        #因此只需要将along_track[i]和h[i]分别记录到along_track_tmp和h_tmp
        #等到完成一段数据的遍历后再进行下一步的计算
        if (along_track[i] < along_start + dis and i != len(along_track) - 1):
            along_track_tmp.append(along_track[i])
            h_tmp.append(h[i])
        
        #如果along_track[i]到along_start的距离小于dis，且i点是最后一个点
        #表明此时虽然该段数据的长度还没有到达极限，但已经是最后一个数据了
        #那么首先要将这个末端点一起保存进along_track_tmp和h_tmp，然后对该段的数据进行信噪分离
        elif (along_track[i] < along_start + dis and i == len(along_track) - 1):
            
            #首先要将这个末端点一起保存进along_track_tmp和h_tmp
            along_track_tmp.append(along_track[i])
            h_tmp.append(h[i])
            cou += 1
            print('finish 1000*' + str(cou))
            
            #更新along_start（其实没啥必要，但是为了保持程序的对称）
            along_start = along_track[i]
            
            #根据论文《Satellite-derived bathymetry using the ICESat-2 lidar and Sentinel-2 imagery datasets》
            #动态计算每一段数据的最小邻域点数阈值MinPts
            if (len(h_tmp) >= 1):
                if((max(h_tmp) - min(h_tmp)) ==0):
                    SN1 = PI * Radius * Radius * len(h_tmp) / (dis * (0.0001))
                else:
                    SN1 = PI * Radius * Radius * len(h_tmp) / (dis * (max(h_tmp) - min(h_tmp)))
                SN2 = PI * Radius * Radius * list_statis(h_tmp) / (dis * h2*2)
                MinPts = (2 * SN1 - SN2) / math.log(2 * SN1 / SN2)
                #为了防止计算的阈值MinPts太小
                if (MinPts < 3):
                    MinPts = 3
                #根据计算的邻域点数量阈值MinPts和设定的邻域范围阈值Radius区分信噪点
                noi_sig_class(along_track_tmp, h_tmp, MinPts, Radius, clss_res)
                # print(MinPts)
            else:
                clss_res.append(0)
        
        #如果along_track[i]到along_start的距离大于dis
        #表明此时还遍历到目前这段数据的尽头了
        #那么就要对该段的数据进行信噪分离，并且更新参数
        else:
            cou += 1
            print('finish 1000*' + str(cou))
            along_start = along_track[i]
            if (len(h_tmp) >= 1):
                if((max(h_tmp) - min(h_tmp)) ==0):
                    SN1 = PI * Radius * Radius * len(h_tmp) / (dis * (0.0001))
                else:
                    SN1 = PI * Radius * Radius * len(h_tmp) / (dis * (max(h_tmp) - min(h_tmp)))
                SN2 = PI * Radius * Radius * list_statis(h_tmp) / (dis * h2*2)
                MinPts = (2 * SN1 - SN2) / math.log(2 * SN1 / SN2)
                if (MinPts < 3):
                    MinPts = 3
                noi_sig_class(along_track_tmp, h_tmp, MinPts, Radius, clss_res)
                # print(MinPts)
            else:
                clss_res.append(0)
            along_track_tmp = [along_track[i]]
            h_tmp = [h[i]]
    
    #将clss_res中保存的信噪结果作为新的属性数据导出到csv中
    df['clss_res'] = clss_res
    df = df.drop(df[df['clss_res'] == 0].index)
    df.to_csv(os.path.join(outputpath, 'pro' + csv_filename))


#该函数用来对上步去噪后的信号点进行进一步的简化
#采取的方法是分段取中位数
#这样可以防止上步去噪后依旧存在一些噪声点
def gather(inputcsv,outputpath):
    formmer_path, csv_filename = os.path.split(inputcsv)
    df = pd.read_csv(inputcsv)
    # df = df.drop(df[df['Signal Confidence'] <= 0].index)
    df = df.sort_values(by='Along-Track (m)', ascending=True)

    lat = list(df['Latitude (deg)'])
    lon = list(df['Longitude (deg)'])
    along_track = list(df['Along-Track (m)'])
    # along_track = [i / 50 for i in along_track]
    h = list(df['Height (m MSL)'])

    step = 10
    dis = 10
    slice_index = list_slice(along_track, step, dis)

    new_alongtrack = []
    new_lat = []
    new_lon = []
    new_h = []
    
    #分段取中位数
    for i in range(len(slice_index)):
        alongtrack_tmp = along_track[slice_index[i][0]:slice_index[i][1]]
        lat_tmp = lat[slice_index[i][0]:slice_index[i][1]]
        lon_tmp = lon[slice_index[i][0]:slice_index[i][1]]
        h_tmp = h[slice_index[i][0]:slice_index[i][1]]

        new_alongtrack.append(np.median(alongtrack_tmp))
        new_lat.append(np.median(lat_tmp))
        new_lon.append(np.median(lon_tmp))
        new_h.append(np.median(h_tmp))

    data_new = {'Latitude': new_lat, 'Longitude': new_lon, 'Along-Track (m)': new_alongtrack, 'Height (m MSL)': new_h}
    df_new = pd.DataFrame(data_new)
    df_new.to_csv(os.path.join(outputpath, 'pro_' + csv_filename), index=False)




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    denoise_path = r'F:\ATL03process\data_process\new_process\de_noise'
    gatherpath = r'F:\ATL03process\data_process\new_process\gather'
    inputpath = r'F:\ATL03process\data_process\test_data'

    for file in os.listdir(inputpath):
        print(file)
        DRAGANN(os.path.join(inputpath, file), denoise_path)
        gather(os.path.join(denoise_path, 'pro'+file), gatherpath)
