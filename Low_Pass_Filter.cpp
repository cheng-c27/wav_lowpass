#include <iostream>
#include <fstream>
#include <string>
#include <math.h> 
using namespace std;
/*
    变量区
*/
#define M_PI (3.1415926535897932384)// Π值
double* pre_filter;                 // 滤波前的数据(范围[-1,1])
double* pre_filter_chan1;           // 滤波前的通道1数据
double* pre_filter_chan2;           // 滤波前的通道2数据
double* after_filter_chan1;         // 滤波后的通道1数据
double* after_filter_chan2;         // 滤波后的通道2数据
double* after_filter;               // 滤波处理后的完整数据
uint8_t* after_filter_demap_8bits;  // 8位滤波后数据反映射结果
uint16_t* after_filter_demap_16bits;// 16位滤波后数据反映射结果
double One_FF = 0.018;              // 单通道低通200Hz的ff = 200 / 11025
double Two_FF = 0.0045;             // 双通道低通200Hz的ff = 200 / 44100

// WAV文件数据结构体
struct wav_struct
{
    unsigned long  ChunkSize;       // 文件大小
    unsigned long  SubChunk1Size;   // FMT数据块的大小
    unsigned short AudioFormat;     // 音频格式说明
    unsigned short NumChannels;     // 通道数
    unsigned long  SampleRate;      // 采样率
    unsigned long  ByteRate;        // Byte率
    unsigned short BlockAlign;      // 数据块对齐单元
    unsigned short BitsPerSample;   // 采样时模数转换的分辨率
    unsigned long  SubChunk2Size;   // 数据块大小
    uint8_t*       data_8bits;      // 8位的数据
    uint16_t*      data_16bits;     // 16位的数据
}WAV;


/*
    WAV文件读取:
    filename：读入的WAV音频文件名
*/
void WavRead(string filename) {
    ifstream fs;                    // 输入文件流对象
    string path = filename + ".wav";// 文件路径（默认与cpp文件在同一目录下）

    // 以二进制的形式读入文件
    fs.open(path, ios::in | ios::binary);

    // 文件打开失败提示
    if (!fs.is_open()) {
        cout << "***文件打开失败,请检查是否有以下情况***" << endl;
        cout << "(1)文件名是否正确。" << endl;
        cout << "(2)文件是否在同一目录下。" << endl;
        return;
    }
    
    // 文件大小
    fs.seekg(0, ios::end);        // 文件移步到末尾
    WAV.ChunkSize = fs.tellg();   // 获取当前文件位置（即末尾处=>文件大小）

    // FMT块大小（16）
    fs.seekg(0x10);
    fs.read((char*)&WAV.SubChunk1Size, sizeof(WAV.SubChunk1Size));

    // 音频格式说明（20）
    fs.seekg(0x14);
    fs.read((char*)&WAV.AudioFormat, sizeof(WAV.AudioFormat));

    // 音频声道信息（22）
    fs.seekg(0x16);
    fs.read((char*)&WAV.NumChannels, sizeof(WAV.NumChannels));

    // 声道采样频率（24）
    fs.seekg(0x18);
    fs.read((char*)&WAV.SampleRate, sizeof(WAV.SampleRate));

    // 码率（28）
    fs.seekg(0x1c);
    fs.read((char*)&WAV.ByteRate, sizeof(WAV.ByteRate));

    // 数据块对齐单元（32）
    fs.seekg(0x20);
    fs.read((char*)&WAV.BlockAlign, sizeof(WAV.BlockAlign));

    // 样本位数（34）
    fs.seekg(0x22);
    fs.read((char*)&WAV.BitsPerSample, sizeof(WAV.BitsPerSample));

    // 数据块的大小（40）
    fs.seekg(0x28);
    fs.read((char*)&WAV.SubChunk2Size, sizeof(WAV.SubChunk2Size));

    // 根据不同的数据位数存入不同类型的数组中
    if (WAV.BitsPerSample == 8) {
        //cout << "数据为8位" << endl;
        WAV.data_8bits = new uint8_t[WAV.SubChunk2Size];

        //读取data块中的数据（44）
        fs.seekg(0x2c);
        fs.read((char*)WAV.data_8bits, WAV.SubChunk2Size * sizeof(uint8_t) );
    }
    else if (WAV.BitsPerSample == 16) {
        //cout << "数据为16位" << endl;
        WAV.data_16bits = new uint16_t[WAV.SubChunk2Size];

        //读取data块中的数据（44）
        fs.seekg(0x2c);
        fs.read((char*)WAV.data_16bits, WAV.SubChunk2Size * sizeof(uint16_t));
    }

    // 关闭文件
    fs.close();

    // 打印文件的基本信息
    cout << "********读入的音频数据********"        << endl;
    cout << "文件大小为    ：" << WAV.ChunkSize     << endl;
    cout << "FMT块大小为   ：" << WAV.SubChunk1Size << endl;
    cout << "音频格式说明  ：" << WAV.AudioFormat   << endl;
    cout << "音频通道数    ：" << WAV.NumChannels   << endl;
    cout << "采样频率      ：" << WAV.SampleRate    << endl;
    cout << "码率          ：" << WAV.ByteRate      << endl;
    cout << "数据块对齐单元：" << WAV.BlockAlign    << endl;
    cout << "样本位数      ：" << WAV.BitsPerSample << endl;
    cout << "音频数据大小  ：" << WAV.SubChunk2Size << endl;
    cout << "******************************"        << endl;

    printf("音频数据读取完毕！\n");
}

/*
    WAV文件输出
    filename： 输出文件的名字
    WAV_INPUT：包含输出文件格式的结构体
    data：     音频数据
*/
template<class T>
void WavWrite(string filename, wav_struct WAV_INPUT, T data[]) {
    ofstream out;// 输出文件流
    string path = filename + ".wav";
    out.open(path, ios::out | ios::binary);

    char RIFF[4] = { 'R','I','F','F' };
    char WAVE[4] = { 'W','A','V','E' };
    char FMT[4]  = { 'f','m','t',' ' };
    char DATA[4] = { 'd','a','t','a' };

    out.write((char*)RIFF, sizeof(RIFF));                                       // ChunkID:"RIFF"
    out.write((char*)&WAV_INPUT.ChunkSize, sizeof(WAV_INPUT.ChunkSize));        // 文件大小
    out.write((char*)WAVE, sizeof(WAVE));                                       // Format:"WAVE"
    out.write((char*)FMT, sizeof(FMT));                                         // SubChunk1ID:"FMT "
    out.write((char*)&WAV_INPUT.SubChunk1Size, sizeof(WAV_INPUT.SubChunk1Size));// fmt块大小
    out.write((char*)&WAV_INPUT.AudioFormat, sizeof(WAV_INPUT.AudioFormat));    // AudioFormat
    out.write((char*)&WAV_INPUT.NumChannels, sizeof(WAV_INPUT.NumChannels));    // 音频通道数
    out.write((char*)&WAV_INPUT.SampleRate, sizeof(WAV_INPUT.SampleRate));      // 音频采样率
    out.write((char*)&WAV_INPUT.ByteRate, sizeof(WAV_INPUT.ByteRate));          // 码率
    out.write((char*)&WAV_INPUT.BlockAlign, sizeof(WAV_INPUT.BlockAlign));      // 数据块对齐单元
    out.write((char*)&WAV_INPUT.BitsPerSample, sizeof(WAV_INPUT.BitsPerSample));// BPS
    out.write((char*)DATA, sizeof(DATA));                                       // SubChunk2ID:"DATA"
    out.write((char*)&WAV_INPUT.SubChunk2Size, sizeof(WAV_INPUT.SubChunk2Size));// 数据块大小
    
    if(WAV_INPUT.BitsPerSample == 8)
        out.write((char*)data, sizeof(uint8_t) * WAV_INPUT.SubChunk2Size);      // 滤波后的数据
    else if(WAV_INPUT.BitsPerSample == 16)
        out.write((char*)data, sizeof(uint16_t) * WAV_INPUT.SubChunk2Size);     // 写入滤波后的数据
    out.close();
    cout << "文件输出完毕!\n" << endl;
}

/*
    采样数据映射：:将数据映射到[-1,1]之间   => 便于处理
    input： 需要映射的数据
    output：映射后存储的数组
    size：  数据大小
    BPS：   样本数据位数
*/
template<class T>
void Map_Data(T* input,double* output, long size,unsigned short BPS) {
    unsigned long data_low;
    double float_data;  // 浮点数据

    if (BPS == 8) {
        for (unsigned long i = 0; i < size; i++)
        {

            if (input[i] > 128 || input[i] == 128) {// >=128时
                data_low = input[i] - 128;
                float_data = (double)(data_low / (double)128);
            }
            else {                                  // <128时
                data_low = 128 - input[i];
                float_data = -(double)(data_low / (double)128);
            }

            output[i] = float_data;// 存入数组

        }
    }
    else if (BPS == 16) {
        for (unsigned long i = 0; i < size; i++)
        {
            if (input[i] > 32768 || input[i] == 32768) {// >=32768时
                data_low = input[i] - 32768;
                float_data = (double)(data_low / (double)32768);
            }
            else {                                      // <32768时
                data_low = 32768 - input[i];
                float_data = -(double)(data_low / (double)32768);
            }

            output[i] = float_data;// 存入数组

        }
    }
    

    //printf("音频数据已映射到[-1,1]之间\n");
}

/*
    采样数据反映射:将[-1,1]之间的数据解映射 => 用于写WAV文件
    input： 需要映射的数据
    output：映射后存储的数组
    size：  数据大小
    BPS：   样本数据位数
*/
template<class T>
void DeMap_Data(double* input, T* output, long size, unsigned short BPS) {
    if (BPS == 8) {
        for (unsigned long i = 0; i < size; i = i + 1) {
            if (input[i] >= 0)
                output[i] = (input[i] * 128 + 127);
            else
                output[i] = (128 - (-input[i] / 1) * 128);
        }
    }
    else if (BPS == 16) {
        for (unsigned long i = 0; i < size; i = i + 1) {
            if (input[i] >= 0)
                output[i] = (input[i] * 32768 + 32767);
            else
                output[i] = (32768 - (-input[i] / 1) * 32768);
        }
    }
}

/*
    巴特沃斯低通滤波器（二阶）函数
    input： 需要处理的数据
    output：处理后的数据
    size：  数据长度
    ff：    截止频率 / 采样频率
*/
void Low_Pass_Filter(double* input, double* output, long size, double ff) {
    // 二阶滤波器参数
    const double ita = 1.0 / tan(M_PI * ff);
    const double q = sqrt(2.0);
    double b0 = 1.0 / (1.0 + q * ita + ita * ita);
    double b1 = 2 * b0;
    double b2 = b0;
    double a1 = 2.0 * (ita * ita - 1.0) * b0;
    double a2 = -(1.0 - q * ita + ita * ita) * b0;
    output[0] = 0.0;
    output[1] = 0.0;

    // 根据差分方程计算滤波后数据
    for (int i = 2; i < size; i++)
        output[i] = b0 * input[i] + b1 * input[i - 1] + b2 * input[i - 2] + a1 * output[i - 1] + a2 * output[i - 2];
}

/*
    操作测试函数
    File_Read： 读入的文件名
    File_Write：输出的文件名
*/
void Test_Process(string File_Read,string File_Write) {
    // 读WAV文件
    WavRead(File_Read);

    // 根据通道数和样本位数分类处理
    if (WAV.BitsPerSample == 8) {   // 8位数据
        // 数据映射
        pre_filter = new double[WAV.SubChunk2Size];
        Map_Data(WAV.data_8bits, pre_filter, WAV.SubChunk2Size,WAV.BitsPerSample);
        delete[]WAV.data_8bits;
        if (WAV.NumChannels == 2) { // 8位双通道数据
            long half_size = (long)(WAV.SubChunk2Size / 2);
            if (half_size % 2 == 1)
                half_size += 1;
            pre_filter_chan1 = new double[half_size];
            pre_filter_chan2 = new double[(long)(WAV.SubChunk2Size / 2)];
            // 双通道数据分开处理
            long c1 = 0, c2 = 0;
            for (int i = 0; i < WAV.SubChunk2Size; i++) {
                if (i % 2 == 0)
                    pre_filter_chan1[c1++] = pre_filter[i];
                else
                    pre_filter_chan2[c2++] = pre_filter[i];
            }
            delete[]pre_filter;
            // 双通道数据分别滤波
            after_filter_chan1 = new double[half_size];
            after_filter_chan2 = new double[(long)(WAV.SubChunk2Size / 2)];
            Low_Pass_Filter(pre_filter_chan1, after_filter_chan1, half_size, One_FF);
            Low_Pass_Filter(pre_filter_chan2, after_filter_chan2, (long)(WAV.SubChunk2Size / 2), One_FF);
            delete[]pre_filter_chan1;
            delete[]pre_filter_chan2;
            // 双通道数据合并输出
            after_filter = new double[WAV.SubChunk2Size];
            c1 = 0, c2 = 0;
            for (int i = 0; i < WAV.SubChunk2Size; i++) {
                if (i % 2 == 0)
                    after_filter[i] = after_filter_chan1[c1++];
                else
                    after_filter[i] = after_filter_chan2[c2++];
            }
            delete[]after_filter_chan1;
            delete[]after_filter_chan2;
        }
        else {// 8位单通道数据
            // 低通滤波
            after_filter = new double[WAV.SubChunk2Size];
            Low_Pass_Filter(pre_filter, after_filter, WAV.SubChunk2Size, One_FF);
            delete[]pre_filter;
        }
        // 反映射
        after_filter_demap_8bits = new unsigned char[WAV.SubChunk2Size];
        DeMap_Data(after_filter, after_filter_demap_8bits, WAV.SubChunk2Size,WAV.BitsPerSample);
        delete[]after_filter;
        // 输出wav文件
        WavWrite(File_Write, WAV, after_filter_demap_8bits);
        delete[]after_filter_demap_8bits;
    }
    else if (WAV.BitsPerSample == 16) { // 16位数据
        // 数据映射
        pre_filter = new double[WAV.SubChunk2Size];
        Map_Data(WAV.data_16bits, pre_filter, WAV.SubChunk2Size,WAV.BitsPerSample);
        delete[]WAV.data_16bits;
        if (WAV.NumChannels == 2) {     // 16位双通道数据
            long half_size = (long)(WAV.SubChunk2Size / 2);
            if (half_size % 2 == 1)
                half_size += 1;
            pre_filter_chan1 = new double[half_size];
            pre_filter_chan2 = new double[(long)(WAV.SubChunk2Size / 2)];
            // 双通道数据分开处理
            long c1 = 0, c2 = 0;
            for (int i = 0; i < WAV.SubChunk2Size; i++) {
                if (i % 2 == 0)
                    pre_filter_chan1[c1++] = pre_filter[i];
                else
                    pre_filter_chan2[c2++] = pre_filter[i];
            }
            delete[]pre_filter;
            // 双通道数据分别滤波
            after_filter_chan1 = new double[half_size];
            after_filter_chan2 = new double[(long)(WAV.SubChunk2Size / 2)];
            Low_Pass_Filter(pre_filter_chan1, after_filter_chan1, half_size, Two_FF);
            Low_Pass_Filter(pre_filter_chan2, after_filter_chan2, (long)(WAV.SubChunk2Size / 2), Two_FF);
            delete[]pre_filter_chan1;
            delete[]pre_filter_chan2;
            // 双通道数据合并输出
            after_filter = new double[WAV.SubChunk2Size];
            c1 = 0, c2 = 0;
            for (int i = 0; i < WAV.SubChunk2Size; i++) {
                if (i % 2 == 0)
                    after_filter[i] = after_filter_chan1[c1++];
                else
                    after_filter[i] = after_filter_chan2[c2++];
            }
            delete[]after_filter_chan1;
            delete[]after_filter_chan2;
        }
        else {// 16位单通道数据
            after_filter = new double[WAV.SubChunk2Size];
            Low_Pass_Filter(pre_filter, after_filter, WAV.SubChunk2Size, Two_FF);// 滤波
            delete[]pre_filter;
        }
        // 反映射
        after_filter_demap_16bits = new uint16_t[WAV.SubChunk2Size];
        DeMap_Data(after_filter, after_filter_demap_16bits, WAV.SubChunk2Size,WAV.BitsPerSample);
        delete[]after_filter;
        // 输出wav文件
        WavWrite(File_Write, WAV, after_filter_demap_16bits);
        delete[]after_filter_demap_16bits;
    }
}

/*           
    主函数        
*/
int main(int argc, char** argv)
{
    string F_R, F_W;

    cout << "请输入读入的文件名和输出的文件名：" << endl;

    while (cin >> F_R >> F_W) {
        Test_Process(F_R, F_W);
        cout << "可继续输入：" << endl;
    }

    return 0;
}