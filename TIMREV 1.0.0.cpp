#include <iostream>
#include <math.h>
#include<bits/stdc++.h>
#include<conio.h>
using namespace std;
class cross_section_properties {
public:
    double Area(double t, double d) {
        double a = t * d;
        return a;
    }
    double Ixx(double t, double d) {
        double Ix = (t * (pow(d, 3))) / 12;
        return Ix;
    }
    double Iyy(double t, double d) {
        double Iy = (d * (pow(t, 3))) / 12;
        return Iy;
    }
    double rxx(double t, double d) {
        double rx = sqrt(Ixx(t, d) / Area(t, d));
        return rx;
    }
    double ryy(double t, double d) {
        double ry = sqrt(Iyy(t, d) / Area(t, d));
        return ry;
    }
    double rbxx(double t, double d, double le,double k1) {
        double rbx = k1* ((le * 100)/d);
        return rbx;
    }
    double rbyy(double t, double d, double le,double k1) {
        double rby = k1 * ((le * 100)/t);
        return rby;
    }
};
class cross_section_modified_strength {
public :
    cross_section_properties csp1;
    double MXu(double t, double d, double fb, double kf, double phi)
    {
    double Ix = csp1.Ixx(t,d);
    double c = d/2;
    double Mmodx = (((fb * Ix)/c) * kf * phi) * (pow(10,-5));
    return Mmodx;
    }
    double MYu (double t, double d,double fb, double kf, double phi )
    {
        double Iy = csp1.Iyy(t,d);
        double c = t/2;
        double Mmody = (((fb * Iy)/c) * kf * phi) * (pow(10,-5));
        return Mmody;
    }
    double FVu (double t, double d, double fv, double kf, double phi)
    {
        double Vu = ((fv * d * t * 2 * kf * phi)/3) * (pow(10,-3));
        return Vu;
    }
    double FCu (double t, double d, double fc, double kf, double phi)
    {
        double a = csp1.Area(t,d);
        double fcu = ((a * fc * kf * phi)) * (pow(10,-3));
        return fcu;
    }
    double FTu (double t, double d, double ft, double kf, double phi)
    {
        double a = csp1.Area(t,d);
        double ftu = ((a * ft * kf * phi)) * (pow(10,-3));
        return ftu;
    }
    double Pcr(double t, double d,double Emin, double le, double k)
    {
        double Ix = csp1.Ixx(t,d);
        double Iy = csp1.Iyy (t,d);
        double pcr;
        if (Ix > Iy)
         pcr = (pow(M_PI,2) * Emin * Iy)/(pow((k*le*100),2)) * (pow(10,-3));
        else
             pcr = (pow(M_PI,2) * Emin * Ix)/(pow((k*le*100),2)) * (pow(10,-3));
        return pcr;
    }
};
class Design_check{
     public :
    double CompressiveForceCheck (double C, double cforce)
    {
        if (C > cforce || C == cforce)
        {
            cout<<"Axial load (compressive) capacity check result (design vs. analysis) :  O. K."<<endl<<endl;
            return 0;
        }
        else
        {
            cout<< "Axial load (compressive) capacity check result (design vs. analysis) :  \n"
                       "WARNING! Member failure eminent due to large compressive force!\n"
                       "Consider redesigning"<<endl<<endl;
        return 0;
        }
            }
            double TensileForceCheck (double T, double tforce)
            {
        if (T > tforce || T == tforce)
                {
                    cout<<"Axial load (tensile) capacity check result (design vs. analysis) :  O. K."<<endl<<endl;
                    return 0;
                }
                else
                {
                    cout<<"Axial load (tensile) capacity check result (design vs. analysis) :  \n"
                       "WARNING! Member failure eminent due to large tensile force!\n"
                       "Consider redesigning"<<endl<<endl;
                       return 0;
                }
            }
    double ShearForceCheck (double V, double vforce)
    {
        if (V > vforce || V== vforce)
        {
            cout<<"Shearing load capacity check result (design vs. analysis) :  O. K."<<endl<<endl;
            return 0;
        }
        else
        {
           cout<<"Shearing load capacity check result (design vs. analysis) :  \n"
               "WARNING! Member failure eminent due to large shearing force!\n"
               "Consider redesigning"<<endl<<endl;
               return 0;
        }
    }
    double MomentCheckX (double Mx, double mx)
    {
        if (Mx > mx || Mx == mx)
        {
            cout<<"Bending moment capacity (around x axis) check result (design vs. analysis) :  O. K."<<endl<<endl;
            return 0;
        }
        else
        {
            cout<<"Bending moment capacity around (x axis) check result (design vs. analysis) :  \n"
               "WARNING! Member failure eminent due to large bending moment!\n"
               "Consider redesigning"<<endl<<endl;
               return 0;
        }
    }
    double MomentCheckY (double My, double my)
    {
        if (My > my || My == my)
        {
            cout<<"Bending moment capacity (around y axis) check result (design vs. analysis) :  O. K."<<endl<<endl;
            return 0;
        }
        else
        {
            cout<<"Bending moment capacity around (y axis) check result (design vs. analysis) :  \n"
                  "WARNING! Member failure eminent due to large bending moment!\n"
                  "Consider redesigning"<<endl<<endl;
                  return 0;
        }
    }
};
int main() {
    double b,h,Fbn,Fvn,Fcn,Ftn,Pcr,Ksl,Rbx,Rby,Le,Emin,Ix,Iy,A,rx,ry,KFb,KFt,KFv,KFc,Phi,C,T,Mx,My,V;
    int i,j = 1;
    int runs,runcount;
    double cmax, tmax, vmax, mxmax, mymax;
    string str,strt,strn;
    Phi =0.75;
    KFb = 2.54;
    KFc = 2.40;
    KFt = 2.70;
    KFv = 2.88;
    double AlaskaPine[5] ={81,44,12,37,35700};
    double AlaskaHemlock [5] ={91,44,13,31,43400};
    double AlaskaNoelle [5] ={98,63,11,23,40600};
    double AlaskaYellowPine [5] ={95,56,16,36,38500};
    double Douglas [5] ={105,70,13,44,48300};
    double SiberiaHemlockPine [5] ={91,54,10,28,43400};
    cout<<"**********************************************WELCOME TO****************************************************************\n"
          "************************************************************************************************************************\n"
          "******                 *****   ****      ******      ****           *******            ****   *****************   ******\n"
          "*************   ************   ****   *   ****   *   ****   *******   *****   **************   ***************   *******\n"
          "*************   ************   ****   **   **   **   ****   ********   ****   ***************   *************   ******** \n"
          "*************   ************   ****   ***   *   **   ****   ******   ******   ****************   ***********   ********* \n"
          "*************   ************   ****   ****      **   ****   *****   *******         ***********   *********   ********** \n"
          "*************   ************   ****   ************   ****   ****   ********   ******************   *******   *********** \n"
          "*************   ************   ****   ************   ****   *****   *******   *******************   *****   ************ \n"
          "*************   ************   ****   ************   ****   ******   ******   ********************   ***   ************* \n"
          "*************   ************   ****   ************   ****   *******   *****            ************   *   ************** \n"
          "************************************************************************************************************************\n"
          "************************************************************************************************************************ \n"
          "************************************************************************************************************************ \n"
          "************************************************************************************************************************ \n"
          "*************A Revolutionary Timber Structures Design Software By Mohammad Amin Moghaddasi******************************"<<endl;
    cout<<"How many members you are planing to run design check on? :"<<endl;
    cin>>runs;
    for(runcount=0;runcount<runs;runcount ++)
    {
    cout<<"Enter the width (b), and height (h) of member cross section (cm) :"<<endl;
    cin>>b>>h;
    cout<<"Enter effective length of member (m):"<<endl;
    cin>>Le;
    cout<<"Enter the slenderness ratio (K) :"<<endl;
    cin>>Ksl;
    cout<<"Please select timber type, using corresponding number for each :"<<endl;
    cout<<"Alaska Pine (1) - Alaska Hemlock (2) - Alaska Noelle (3)\n"
          "Alaska Yellow Pine (4) - Douglas (5) - Siberia Hemlock Pine (Yolka) (6): "<<endl;
    cin>>i;
        while (j == 1) {
            if (i < 1 || i > 6) {
                cout<<endl<<"Your selection does not exist in the list. Please reselect timber type.\n "
                        "Alaska Pine (1) - Alaska Hemlock (2) - Alaska Noelle (3)\n"
                        "Alaska Yellow Pine (4) - Douglas (5) - Siberia Hemlock Pine (Yolka) (6) :"<<endl;
                cin>>i;
            } else
                j = 2;
        }
        if(i == 1){
            Fbn = AlaskaPine[0];
            Ftn = AlaskaPine[1];
            Fvn = AlaskaPine[2];
            Fcn = AlaskaPine[3];
            Emin = AlaskaPine[4];
        }else if(i == 2){
            Fbn = AlaskaHemlock[0];
            Ftn = AlaskaHemlock[1];
            Fvn = AlaskaHemlock[2];
            Fcn = AlaskaHemlock[3];
            Emin = AlaskaHemlock[4];
        }else if(i == 3){
            Fbn = AlaskaNoelle[0];
            Ftn = AlaskaNoelle[1];
            Fvn = AlaskaNoelle[2];
            Fcn = AlaskaNoelle[3];
            Emin = AlaskaNoelle[4];
        }else if(i == 4){
            Fbn = AlaskaYellowPine[0];
            Ftn = AlaskaYellowPine[1];
            Fvn = AlaskaYellowPine[2];
            Fcn = AlaskaYellowPine[3];
            Emin = AlaskaYellowPine[4];
        }else if(i == 5){
            Fbn = Douglas[0];
            Ftn = Douglas[1];
            Fvn = Douglas[2];
            Fcn = Douglas[3];
            Emin = Douglas[4];
        }else if(i == 6){
            Fbn = SiberiaHemlockPine[0];
            Ftn = SiberiaHemlockPine[1];
            Fvn = SiberiaHemlockPine[2];
            Fcn = SiberiaHemlockPine[3];
            Emin = SiberiaHemlockPine[4];
        }
        else
        return 0;
        cout<<"Enter member name and type :"<<endl;
        cin>>strn>>strt;
   cout<<"Enter maximum Compressive load value (tonf):"<<endl;
   cin>>cmax;
   cout<<endl;
    cout<<"Enter maximum Tensile load value (tonf):"<<endl;
    cin>>tmax;
    cout<<endl;
    cout<<"Enter maximum Shearing load value (tonf):"<<endl;
   cin>>vmax;
   cout<<endl;
    cout<<"Enter maximum Bending moment value(around x axis) (tonf.m) :"<<endl;
    cin>>mxmax;
    cout<<endl;
    cout<<"Enter maximum Bending moment value(around y axis) (tonf.m) :"<<endl;
    cin>>mymax;
    cout<<endl;
        cross_section_properties csp;
    A = csp.Area(b,h);
    Ix = csp.Ixx(b,h);
    Iy = csp.Iyy(b,h);
    rx = csp.rxx(b,h);
    ry = csp.ryy(b,h);
    Rbx = csp.rbxx(b,h,Le,Ksl);
    Rby = csp.rbyy(b,h,Le,Ksl);
     cout<<endl<<"Design control results for a "<<strt<<" named "<<strn<<" :"<<endl;
    cout<<"The cross section area is : "<<A<< " cm^2"<<endl;
    cout<<"The moment of inertia around x is : "<<Ix<<" cm^4"<<endl;
    cout<<"The moment of inertia around y is : "<<Iy<<" cm^4"<<endl;
    cout<<"The gyration radius around x is : "<<rx<<" cm"<<endl;
    cout<<"The gyration radius around y is : "<<ry<<" cm"<<endl;
    if (Rbx > Rby)
        cout<<"Bending will occur around x Axis"<<endl;
        else if (Rby > Rbx)
            cout<<"Bending will occur around y Axis"<<endl;
    else
        cout<<"Bending can occur around either both Axises (X and Y)"<<endl;
    if (50<Rbx && Rbx<75)
    cout<<"The slenderness value around x is : "<<Rbx<<" WARNING!"<<endl<<"The slenderness value is higher than 50!\n"
                                                                          " but it is O. K. for construction (Check Design according to A.I.T.C)"<<endl;
    else if (Rbx>75)
        cout<<"The slenderness value around x is : "<<Rbx<<" WARNING!"<<endl<<"The slenderness value is higher than 75!\n"
                                                                              " it is not suitable neither for design nor for construction (Redesign according to A.I.T.C )"<<endl;
    else
        cout<<"The slenderness value around x is : "<<Rbx<<" This value is O. K."<<endl;
    if (50<Rby && Rby<75)
        cout<<"The slenderness value around y is : "<<Rby<<" WARNING!"<<endl<<"The slenderness value is higher than 50!\n"
                                                                              "but it is O. K. for construction (Check Design)"<<endl;
    else if (Rby>75)
        cout<<"The slenderness value around y is : "<<Rby<<" WARNING!"<<endl<<"The slenderness value is higher than 75!\n"
                                                                              " it is not suitable neither for design nor for construction (Modify Design)"<<endl;
    else
        cout<<"The slenderness value around y is : "<<Rby<<" This value is O. K."<<endl;
    cross_section_modified_strength csms;
    if (csms.Pcr(b,h,Emin,Le,Ksl)== csms.FCu(b,h,Fcn,KFc,Phi)) {
        C = csms.FCu(b, h, Fcn, KFc, Phi);
         str = "Critical axial force is not dominant";
    }
    else {
        C = csms.Pcr(b, h, Emin, Le, Ksl);
        str = "Critical axial force is dominant";
    }
        T = csms.FTu(b,h,Ftn,KFt,Phi);
    V = csms.FVu(b,h,Fvn,KFv,Phi);
    Mx = csms.MXu(b,h,Fbn,KFb,Phi);
    My = csms.MYu(b,h,Fbn,KFb,Phi);
    Pcr = csms.Pcr(b,h,Emin,Le,Ksl);
    cout<<"The ultimate bending moment capacity of cross section around x is :"<<Mx<<" tonf.m"<<endl;
    cout<<"The ultimate bending moment capacity of cross section around y is :"<<My<<" tonf.m"<<endl;
    cout<<"The ultimate shear capacity of cross section is :"<<V<<" tonf"<<endl;
    cout<<"The ultimate tension capacity of cross section is :"<<T<<" tonf"<<endl;
    cout<<"The ultimate compression capacity of cross section is :"<<C<<" tonf"<<endl;
    cout<<str<<endl;
    cout<<"Critical axial load capacity of member is :"<<Pcr<<" tonf"<<endl<<endl;
    cout<<"****************************************************************"<<endl;
    cout<<"Design vs. Analysis Check Results :"<<endl;
    Design_check dcheck;
    dcheck.CompressiveForceCheck(C,cmax);
    dcheck.TensileForceCheck(T,tmax);
    dcheck.ShearForceCheck(V,vmax);
    dcheck.MomentCheckX(Mx,mxmax);
    dcheck.MomentCheckY(My,mymax);
    cout<<endl<<"**********************************************************"<<endl;
    }
    cout<<endl<<"Press Any Key to Continue ......."<<endl;
    _getch();
}
