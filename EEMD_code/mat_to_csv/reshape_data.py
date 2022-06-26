import csv
import os
from os import path

import numpy as np
import scipy.io as scio

def reshape_data(url_in,key1,key_ref,ref_value,url_out,outFileName):

    file1 = os.listdir(url_in)
    #file1.sort(key=lambda x: int(x[:-5]))
    #print(file1)
    for f in file1:
        if os.path.isdir(url_in): #判断是否是文件夹
            url2 = path.join(url_in, f)
            #print(url2)
            file2 = os.listdir(url2)
            # 对每个文件名将扩展名前的字符串转化为数字，然后以数字为key来进行排序，这样便能按照我们的心意来排序了
            file2.sort(key=lambda x: int(x[:-4]))
            for f in file2:
                real_url_in = path.join(url2, f)
                #print(real_url_in)
                real_url_out = path.join(url_out, outFileName)
                #读取指定信号和reference
                data_in = scio.loadmat(real_url_in)
                specific_signal = data_in[key1]
                specific_signal_ref = data_in[key_ref]
                reference_value = data_in[ref_value]

                temp_signal = np.zeros(np.size(specific_signal_ref))
                #信号大于ref则将信号截取为ref长度
                if np.size(specific_signal) > np.size(specific_signal_ref):

                    for i in range(np.size(temp_signal)):
                        temp_signal[i] = specific_signal[0,i]
                #信号小于ref则将信号补均值至长度等于ref
                elif np.size(specific_signal) < np.size(specific_signal_ref):

                    for i in range(np.size(specific_signal)):
                        temp_signal[i] = specific_signal[0,i]
                    for j in range(np.size(specific_signal_ref))[np.size(specific_signal):]:
                        temp_signal[j] = np.mean(specific_signal)
                #信号长度等于ref
                else:
                    for i in range(np.size(temp_signal)):

                        temp_signal[i] = specific_signal[0,i]


                #写入文件
                csvfile = open(real_url_out, 'a', newline='')

                writer = csv.writer(csvfile)
                #writer.writerow(temp_signal)
                writer.writerow(reference_value)


            csvfile.close()

#数据文件夹地址
url_in = r'C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480s切片数据含EP'
#要转换数据名
key1 = 'co2'
#reference数据名称
key_ref = 'co2'
#reference value
ref_value = 'reference_value'
#输出数据文件夹地址
url_out = r'C:\Users\Administrator\Desktop\毕业设计\数据\初始数据\EEDM提取数据\480切片数据SQI优化csv'
#数据文件名称
outFileName = 'reference_value.csv'
reshape_data(url_in,key1,key_ref,ref_value,url_out,outFileName)









