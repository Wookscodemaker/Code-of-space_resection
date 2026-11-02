#include <iostream>
#include<Eigen/Dense>
#include<cmath>
#include<vector>
#include<fstream>
#include<iomanip>//为了调用cout<<fixed<<setprecision(n);来设置输出精度
using namespace std;
class wk{//wellknow，已知点。
public:
    int n;
    double f;
    double m;
    vector<double> cp;//cp是controllpoint，控制点。
};
class ae{// angle_element角元素
public:
    double fi;//phi角的坐标，航偏角。
    double w;//omega角的坐标，旁偏角。
    double k;//kappa角的坐标，像转角。
    void init_ae(){
        fi=0;
        w=0;
        k=0;
    }
};
class ces{//coordinate_element of s,摄影中心的坐标
public:
    double xs;//s点的x坐标
    double ys;//s点的y坐标
    double zs;//s点的z坐标
    void init_ces(wk fc){//fc为这个函数专用的名称
        int j;
        int x_index=3;
        int y_index=4;
        int z_index=5;
        for(j=0;j<fc.n;j++){
            xs+=fc.cp[j*6+x_index];
            ys+=fc.cp[j*6+y_index];
            zs+=fc.cp[j*6+z_index];
        }
        xs/=fc.n;   
        ys/=fc.n;
        zs/=fc.n;
        int h=fc.f/1000*fc.m;
        zs=zs+h;
    }
};

int main()
{
    wk zx;//无任何意义，只是为了引用一个变量。
    string path;
    cout << "请输入文件路径(使用＼进行区分): ";
    cin >> path; // 用户手动在终端输入
    ifstream fin(path); // ifstream 可以直接接受 string
    if (!fin.is_open()) {
        cout << "打开失败！" << endl;
        return 1;
    }
    fin >>zx.n>>zx.f>>zx.m;
    zx.cp.resize(zx.n*6);
    for(int i=0;i<zx.n;i++){
        int l=i*6;
        fin >> zx.cp[l] >> zx.cp[l+1] >> zx.cp[l+2] >> zx.cp[l+3] >> zx.cp[l+4] >> zx.cp[l+5];
    }
    if(!fin){
        cout << "文件格式有误！" << endl;
        return 1;
    }else{
        cout << "输入成功!"<<endl;
    }
//以上完成数据的输入，下方将对所输入的值进行初始化。
    ces ps;//ps是摄影中心的坐标
    ae fa;//三个和角度相关的内方位元素
    ps.init_ces(zx);
    fa.init_ae();
    cout << "初始化成功！" << endl;

}