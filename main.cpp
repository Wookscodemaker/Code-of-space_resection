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
class Rmatrix{
public:
    double a1;
    double a2;
    double a3;
    double b1;
    double b2;
    double b3;
    double c1;
    double c2;
    double c3;
    void c_r_m(ae kt){//calculater of rotation matrix旋转矩阵计算器。
        a1=cos(kt.fi)*cos(kt.k)-sin(kt.fi)*sin(kt.w)*sin(kt.k);
        a2=(-cos(kt.fi)*sin(kt.k))-sin(kt.fi)*sin(kt.w)*cos(kt.k);
        a3=(-sin(kt.fi)*cos(kt.w));
        b1=cos(kt.w)*sin(kt.k);
        b2=cos(kt.w)*cos(kt.k);
        b3=(-sin(kt.w));
        c1=sin(kt.fi)*cos(kt.k)+cos(kt.fi)*sin(kt.w)*sin(kt.k);
        c2=(-sin(kt.fi)*sin(kt.k))+cos(kt.fi)*sin(kt.w)*cos(kt.k);
        c3=cos(kt.fi)*cos(kt.w);
    }
};
void x0y0z0(int n,vector<double>xb,vector<double> yb,vector<double> zb,Rmatrix rm,vector<double>x0,vector<double> y0,vector<double> z0){
    vector<double> xyz;
    double x,y,z;
    for(int j=0;j<n;j++){
        x=rm.a1*xb[j]+rm.a2*yb[j]+rm.a3*zb[j];
        x0.push_back(x);
        y=rm.b1*xb[j]+rm.b2*yb[j]+rm.b3*zb[j];
        y0.push_back(y);
        z=rm.c1*xb[j]+rm.c2*yb[j]+rm.c3*zb[j];
        z0.push_back(z);
    }
}
double f_x(double f,double x0,double z0){
    return -f*x0/z0;
}
double f_y(double f,double y0,double z0){
    return -f*y0/z0;
}
vector<double> cheif_calculator(wk xk,ces sp,ae av){//xk为读取的文件值，sp是是目前的坐标，av是角目前已知的值。
    vector<double> A;
    vector<double> L;
    vector<double> D;

    Rmatrix R;
    R.c_r_m(av);//计算a1、a2、a3、b1、b2、b3、c1、c2、c3的值，放入类R中
    vector<double> xb;
    vector<double> yb;
    vector<double> zb;
    for(int i=0;i < xk.n;i++){
        int l=i*6;
        xb.push_back(xk.cp[l+3]-sp.xs);
        yb.push_back(xk.cp[l+4]-sp.ys);
        zb.push_back(xk.cp[l+5]-sp.zs);
    }
    //前面是计算x-xs，y-ys，z-zs的值
    

}
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
    //下一步通过循环，进行最后的求解
   


}