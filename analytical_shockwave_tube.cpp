#include <iostream>
#include <cmath>
#include <fstream>

double solve_by_newton_method(double x, double gamma0, double gamma4, double a4, double a0, double p4, double p0, double u4, double u0);

int main() {
    //時間は自分で指定
    double t = 2;

    const double gamma0 = 1.4;
    const double gamma4 = 1.4;

    //格子点の大きさの指定
    const int mx = 101;
    const double xmin = -5.0;
    const double xmax = 5.0;
    const double x0 = 0;
    const double dx = (xmax - xmin) / (mx - 1);
    //座標に指定
    double x[mx];
    for (int i = 0; i < mx; i++) {
        x[i] = xmin + i * dx;
    }

    const double rho4 = 2.0;
    const double p4 = 0.8;
    const double u4 = 0.0;  //初期値の速度は０であることを固定
    const double a4 = sqrt(gamma4 * p4 / rho4);

    const double rho0 = 1;
    const double p0 = 0.4;
    const double u0 = 0.0;//初期値の速度は０であることを固定
    const double a0 = sqrt(gamma0 * p0 / rho0);

    double Ms;
    double y;
    double Msini;

    //各領域でのvariableを計算
    double rho_const[5];
    double u_const[5];
    double p_const[5];

    //まず、constな値を求める。
    //P
    Msini = 1.0; //初期値

    Ms = solve_by_newton_method(Msini, gamma0, gamma4, a4, a0, p4, p0, u4, u0);
    double Us = Ms*a0;
    std::cout << "Ms = " << Ms << std::endl;
    //領域０
    rho_const[0] = rho0;
    u_const[0] = u0;
    p_const[0] = p0;

    //領域１
    p_const[1] = p0*(2*gamma0*pow(Ms,2)-(gamma0-1))/(gamma0+1);
    u_const[1] = a0*2/(gamma0+1)*(Ms-1/Ms);
    rho_const[1] = rho0*((gamma0+1)*pow(Ms,2))/((gamma0-1)*pow(Ms,2)+2);

    //領域２
    p_const[2] = p_const[1];
    u_const[2] = u_const[1];
    double a2 = a4-(gamma4-1)*u_const[2]/2;
    rho_const[2] = gamma4*p_const[2]/pow(a2,2);

    //領域4
    //ここは膨張波領域なので値が定まらない

    //領域5
    rho_const[4] = rho4;
    u_const[4] = u4;
    p_const[4] = p4;

    //全ての格子点でのu,v,pを格納する
    double u[mx];
    double rho[mx];
    double p[mx];

    for (int i = 0; i<mx; i++){
        //領域０
        if(x[i] > x0 + Us*t){
            u[i] = u_const[0];
            rho[i] = rho_const[0];
            p[i] = p_const[0];
        }

        //領域１
        if(x[i] > x0 + u_const[1]*t and x[i] <= x0 + Us*t){
            u[i] = u_const[1];
            rho[i] = rho_const[1];
            p[i] = p_const[1];
        }

        //領域２
        if(x[i] > x0 + (u_const[2]-a2)*t and x[i] <= x0 + u_const[1]*t){
            u[i] = u_const[2];
            rho[i] = rho_const[2];
            p[i] = p_const[2];
        }

        //領域３
        if(x[i] > x0 - a4*t and x[i] <= x0 + (u_const[2]-a2)*t){
            double a = a4*(2/(gamma4+1) - (gamma4-1)/(gamma4+1)*(x[i]/(a4*t)));
            u[i] = 2*a4/(gamma4+1)*(1+(x[i]/(a4*t)));
            p[i] = p4 * pow(a/a4,2*gamma4/(gamma4-1));
            rho[i] = rho4*pow(a/a4,2/(gamma4-1));
        }

        //領域４
        if( x[i] <= x0 - a4*t){
            u[i] = u_const[4];
            rho[i] = rho_const[4];
            p[i] = p_const[4];
        }



    }

    //gnuplot用にファイル出力
    std::ofstream ofs("analytical_data.data");
    if (!ofs) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }

    for (int i = 0; i<mx ; i++){
        ofs << x[i] << " " << u[i] << " " << rho[i] << " " << p[i] << std::endl;
    }

    //csvファイルを出力
    std::ofstream ofscsv("analytical_data.csv");    //, std::ios::app
    if (!ofscsv) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }

    for (int i = 0; i<mx ; i++){
        ofscsv << p0 << "," << p4 << "," << rho0 << "," << rho4 << "," << x[i] << "," << t << "," << u[i] << "," << rho[i] << "," << p[i] << std::endl;
    }


    //test用ファイル出力
    std::ofstream ofstest("test2.csv");  //, std::ios::app
    if (!ofstest) {
        std::cerr << "ファイルオープンに失敗" << std::endl;
        std::exit(1);

    }

    for (int i = 0; i<mx ; i++){
        ofstest << p0 << "," << p4 << "," << rho0 << "," << rho4 << "," << x[i] << "," << t << std::endl;
    }

    return 0;

}


double solve_by_newton_method(double x, double gamma0, double gamma4, double a4, double a0, double p4, double p0, double u4, double u0) {

    auto func = [gamma0, gamma4, a4, a0, p4, p0, u4, u0](double Ms) {
        return (2*gamma0*pow(Ms,2)-(gamma0-1))/(gamma0+1)*pow((1-(gamma4-1)/(gamma0+1)*a0/a4*(Ms-1/Ms)),-2*gamma4/(gamma4-1))-p4/p0;
    };

    double y1;
    double y2;
    double x1;
    double x2;
    double diff_x;
    double error = 0.0001;
    double delta = 0.00001;
    double next_x;
    int cnt = 0;
    x1 = x;

    while (1) {

        x2 = x1 + delta;

        y1 = func(x1);
        y2 = func(x2);

        diff_x = (y2 - y1) / (x2 - x1);
        next_x = x1 - y1 / diff_x;

//        if (cnt <10){
//            std::cout << y1 << ',' << y2 << ',' << next_x << std::endl;
//        }

        if (fabs(y1 - 0) < error) {
            break;
        }
        x1 = next_x;
        cnt += 1;
    }

    return x1;

}