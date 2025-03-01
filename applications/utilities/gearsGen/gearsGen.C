#include<fstream>
#include<cmath>
#include<iostream>

using namespace std;




double xr(double qr, double umin, double lx, double ly, double rt, double rb, double rp, double q0, double u_b, double V, double qrmax, double qamin)
{
    double pi, firQMax, XrVQMax, YrVQMax, RotR, XaQMax, RotA, RotRA, fir, xr, xrv, yrv, YaQMax;

    pi=atan(1)*4;

    firQMax = 1/rp*((-lx+rt*sin(qrmax))/cos(u_b/180*pi)+(V-ly-rt*cos(qrmax))*tan(qrmax)*cos(u_b/180*pi));

    XrVQMax=cos(firQMax)*(-lx+rt*sin(qrmax))/cos(u_b/180*pi)+sin(firQMax)*(-ly-rt*cos(qrmax))+sin(firQMax)*rp+sin(firQMax)*V-rp*firQMax*cos(firQMax);
    
    YrVQMax=-sin(firQMax)*(-lx+rt*sin(qrmax))/cos(u_b/180*pi)+cos(firQMax)*(-ly-rt*cos(qrmax))+cos(firQMax)*rp+cos(firQMax)*V+rp*firQMax *sin(firQMax);

    RotR=(abs(YrVQMax)/(YrVQMax))*atan(abs(YrVQMax/XrVQMax));

    if (XrVQMax<0) RotR=pi-RotR;

    XaQMax=rb*sin(qamin+q0)-rb*qamin*cos(qamin+q0);

    YaQMax=rb*cos(qamin+q0)+rb*qamin*sin(qamin+q0);


    RotA= (abs(YaQMax)/(YaQMax))*atan(abs(YaQMax/XaQMax));

    if (XaQMax<0) RotA=pi-RotA;

    RotRA=RotA-RotR;

    fir=1/rp*((-lx+rt*sin(qr))/cos(u_b/180*pi)+(V-ly-rt*cos(qr))*tan(qr)*cos(u_b/180*pi));

    xrv=cos(fir)*(-lx+rt*sin(qr))/cos(u_b/180*pi)+sin(fir)*(-ly-rt*cos(qr))+sin(fir)*rp+sin(fir)*V-rp*fir*cos(fir);

    yrv=-sin(fir)*(-lx+rt*sin(qr))/cos(u_b/180*pi)+cos(fir)*(-ly-rt*cos(qr))+cos(fir)*rp+cos(fir)*V+rp*fir*sin(fir);

    xr=xrv*cos(RotRA)-yrv*sin(RotRA);

    return xr;
}


double yr
(
    double qr,
    double umin,
    double lx,
    double ly,
    double rt,
    double rb,
    double rp,
    double q0,
    double u_b,
    double V,
    double qrmax,
    double qamin
)
{
    double pi, firQMax, XrVQMax, YrVQMax, RotR, XaQMax, YaQMax, RotA, fir, xrv, yr, yrv, RotRA;
    pi=atan(1)*4;
    firQMax=1/rp*((-lx+rt*sin(qrmax))/cos(u_b/180*pi)+(V-ly-rt*cos(qrmax))*tan(qrmax)*cos(u_b/180*pi));

    XrVQMax=cos(firQMax)*(-lx+rt*sin(qrmax))/cos(u_b/180*pi)+sin(firQMax)*(-ly-rt*cos(qrmax))+sin(firQMax)*rp+sin(firQMax)*V-rp*firQMax* cos(firQMax);
    YrVQMax=-sin(firQMax)*(-lx+rt*sin(qrmax))/cos(u_b/180*pi)+cos(firQMax)*(-ly-rt*cos(qrmax))+cos(firQMax)*rp+cos(firQMax)*V+rp*firQMax*sin(firQMax);

    RotR=(abs(YrVQMax)/(YrVQMax))*atan(abs(YrVQMax/XrVQMax));

    if (XrVQMax<0) RotR=pi-RotR;

    XaQMax=rb*sin(qamin+q0)-rb*qamin*cos(qamin + q0);
    YaQMax=rb*cos(qamin+q0)+rb*qamin*sin(qamin + q0);

    RotA=(abs(YaQMax)/(YaQMax))*atan(abs(YaQMax/XaQMax));
	
    if (XaQMax<0) RotA=pi-RotA;

    RotRA=RotA-RotR;

    fir=1/rp*((-lx+rt*sin(qr))/cos(u_b/180*pi)+(V-ly-rt*cos(qr))*tan(qr)*cos(u_b/180*pi));
    xrv=cos(fir)*(-lx+rt*sin(qr))/cos(u_b/180*pi)+sin(fir)*(-ly-rt*cos(qr))+sin(fir)*rp+sin(fir)*V-rp*fir*cos(fir);
    yrv=-sin(fir)*(-lx+rt*sin(qr))/cos(u_b/180*pi)+cos(fir)*(-ly-rt*cos(qr))+cos(fir)*rp+cos(fir)*V+rp*fir*sin(fir);

    yr=xrv*sin(RotRA)+yrv*cos(RotRA);
    return yr;

}



double xrs
(
    double qr,
    double umin,
    double lx,
    double ly,
    double rt,
    double rb,
    double rp,
    double q0,
    double u_b,
    double V,
    double qrmax,
    double qamin,
    int z
)
{ 
    double pi, xrv, yrv, xrs;

    xrv=xr(qr, umin, lx, ly, rt, rb, rp, q0, u_b, V, qrmax, qamin);

    yrv=yr(qr, umin, lx, ly, rt, rb, rp, q0, u_b, V, qrmax, qamin);

    pi=atan(1)*4;

    xrs=xrv*cos(pi/z)-yrv*sin(pi/z);

    return xrs;
}

double yrs
(
    double qr,
    double umin,
    double lx,
    double ly,
    double rt,
    double rb,
    double rp,
    double q0,
    double u_b,
    double V,
    double qrmax,
    double qamin,
    int z
)
{
    double pi, xrv, yrv, yrs;
    xrv=xr(qr,umin,lx,ly,rt,rb,rp,q0,u_b,V,qrmax,qamin);
    yrv=yr(qr,umin,lx,ly,rt,rb,rp,q0,u_b,V,qrmax, qamin);
    pi = atan(1)*4;
    yrs=xrv*sin(pi/z)+yrv*cos(pi/z);
    return yrs;

}

int main()
{

    int z1_, z2_, i, t, nPoints;
    double pi, mn, an, b, jn1, jn2, ca1_, ca2_, cd1_, cd2_, l1_, l2_, adf, x1_, x2_, s_1, s_2, u, mt, at, dp1_, dp2_, a, bb, xg1, xg2, zn1; double da1_, da2_, rt1, rt2; 
    double zn2, st1, st2, db1_, db2_, atdf, K, D31, E31, df1_, df2_, dtif1, dtif2, dpdf1, dpdf2, bdf, andf, stdf1, stdf2, adf1, adf2, ddf1, ddf2, rb1_, rb2_, q0p, q0r;
    double x1sx, y1sx, x1dx, y1dx;
    double q, qrmax1, u1min, lcx1, lcy1;
    double alpha, gamma, beta, dAngle, x, y, xmod, ymod, angle, r;
    double punti[153][2];

    ofstream outputFile;
    outputFile.open("./punti.dat");
    ofstream outputFile2;
    outputFile2.open("./profilo.dat");
	
    pi=3.141592654;
    nPoints=0;

    // PARAMETRI DA MODIFICARE

    ifstream DatiRuota;
    DatiRuota.open("gearsGenDict");
    cerr <<" ---------------------------------------- "<< endl;
    cerr <<" ----- Gear profile Generation tool ----- "<< endl;
    cerr <<" ---------------------------------------- "<< endl;
    cerr <<" DATI RUOTA "<< endl;

    DatiRuota >> z1_;
    cerr <<" numero di denti [] = "<< z1_ << endl;
    z2_=z1_;
    DatiRuota >> mn;
    cerr <<" modulo normale [mm] = "<< mn << endl;
    DatiRuota >> an;
    cerr <<" angolo di pressione normale [deg] = "<< an << endl;
    DatiRuota >> b;
    cerr <<" angolo d'elica [deg]= "<< b << endl;
    DatiRuota >> jn1;
    cerr <<" contributo gioco normale [mm] = "<< jn1 << endl;
    DatiRuota >> ca1_;
    cerr <<" coefficiente addendum creatore [] = "<< ca1_ << endl;
    DatiRuota >> cd1_;
    cerr <<" coefficiente dedendum creatore [] = "<< cd1_ << endl;
    DatiRuota >> l1_;
    cerr <<" fascia utile [mm] = "<< l1_ << endl;
    DatiRuota >> adf;
    cerr <<" interasse di funzionamento [mm] = "<< adf << endl;
    DatiRuota >> x1_;
    cerr <<" coefficiente spostamento a gioco nullo [] = "<< x1_ << endl;
    DatiRuota >> s_1;
    cerr <<" smusso di testa [] = "<< s_1 << endl;



    /*
    // numero di denti
    z1_=23;
    z2_=23;
    //cerr <<"numero di denti 1 = "<< z1_ << " numero di denti 2 = "<< z2_ << endl;

    // modulo normale
    mn=4;
    //cerr <<"modulo normale = "<< mn << endl;

    // angolo di pressione normale
    an=20;
    //cerr <<"angolo di pressione normale = "<< an << endl;

    // angolo d'elica
    b=0;
    //cerr <<"beta = "<< b << endl;

    // contributo gioco normale
    jn1=0.1;
    //jn2=0.1;
    //cerr <<"contributo gioco normale 1 = "<< jn1 << " contributo gioco normale 2 = "<< jn2 << endl;

    // Coeff add creatore
    ca1_=1.25;
    //ca2_=1.25;
    //cerr <<"Coeff add creatore 1 = "<< ca1_ << " Coeff add creatore 2 = "<< ca2_ << endl;

    // Coefficiente ded creatore
    cd1_=1.51511;
    //cd2_=1.51511;
    //cerr <<"Coefficiente ded creatore 1 = "<< cd1_ << " Coefficiente ded creatore 2 = "<< cd2_ << endl;

    // fascia utile
    l1_=20;
    //l2_=20;
    //cerr <<"fascia utile 1 = "<< l1_ << " fascia utile 2 = "<< l2_ << endl;

    // interasse di funzionamento
    adf=91.5;
    //cerr <<"interasse di funzionamento = "<< adf << endl;

    // coeff so.to a gioco nullo
    x1_=-0.2;
    //x2_=0;
    //cerr <<"coeff so.to a gioco nullo 1 = "<< x1_ << " coeff so.to a gioco nullo 2 = "<< x2_ << endl;

    //smusso di testa
    s_1=0;
    //s_2=0;
    */
    //-----------------------------------------------------------------------------------------------------------------------------------

    // PARAMETRI CALCOLATI

    // Rapporto di trasmissione
    u=z2_/z1_;
    //cerr <<"Rapporto di trasmissione = "<< u << endl;

    // Raggio testa creatore
    rt1=0.3*mn;
    rt2=0.3*mn;
    //cerr <<"Raggio testa creatore 1 = "<< rt1 << " Raggio testa creatore 2 = "<< rt2 << endl;

    // Modulo trasv. teorico
    mt=mn/cos(b/180*pi);
    //cerr <<"modulo trasv. teorico = "<< mt << endl;

    // angolo di pressione trasv. teroico
    at=atan(tan(an/180*pi)/cos(b/180*pi))/pi*180;
    //cerr <<"angolo di pressione trasv. teroico = "<< at << endl;

    // diam. prim. teorico
    dp1_=mt*z1_;
    dp2_=mt*z2_;
    //cerr <<"diam. prim. teorico 1 = "<< dp1_ << " diam. prim. teorico 2 = "<< dp2_ << endl;

    // interasse teorico
    a=(dp1_+dp2_)/2.0;
    //cerr <<"interasse teorico = "<< a << endl;

    // angolo d'elica base
    bb=asin(sin(b/180*pi)*cos(an/180*pi))/pi*180;
    //cerr <<"angolo d'elica base = "<< bb << endl;

    // Coef. sp.to con gioco
    xg1=x1_-jn1/2/mn/sin(an/180*pi);
    xg2=x2_-jn2/2/mn/sin(an/180*pi);

    cout<< "WARNING: uninitialized value of jn2, xg2=0 ALWAYS" << "\n";
    //cerr <<"Coef. sp.to con gioco 1 = "<< xg1 << " Coef. sp.to con gioco 2 = "<< xg2 << endl;

    // Numero denti virtuale
    zn1=z1_/(cos(bb/180*pi)*cos(bb/180*pi))/cos(b/180*pi);
    zn2=z2_/(cos(bb/180*pi)*cos(bb/180*pi))/cos(b/180*pi);
    //cerr <<"Numero denti virtuale 1 = "<< zn1 << " Numero denti virtuale 2 = "<< zn2 << endl;

    // Sp. trasv. su d. p. teo.
    st1=(pi/2+2*xg1*tan(an/180*pi))*mt;
    st2=(pi/2+2*xg2*tan(an/180*pi))*mt;
    //cerr <<"Sp. trasv. su d. p. teo. 1 = "<< st1 << " Sp. trasv. su d. p. teo. 2 = "<< st2 << endl;

    // diametro base
    db1_=dp1_*cos(at/180*pi);
    db2_=dp2_*cos(at/180*pi);
    rb1_=db1_/2;
    rb2_=db2_/2;
    //cerr << "diametro di base 1 = "<< db1_ << "diametro di base 2 = "<< db2_ <<endl;

    // A. di p.  trasv. di f.
    atdf=acos((db1_+db2_)/2/adf)/pi*180;
    //cerr <<"A. di p.  trasv. di f. = "<< atdf << endl;

    // coeff. riavv. interasse
    K=(z1_+z2_)/2*(((tan(atdf/180*pi)-(atdf/180*pi))-(tan(at/180*pi)-(at/180*pi)))/tan(an/180*pi)-(cos(at/180*pi)/cos(atdf/180*pi)-1)/cos(b/180*pi));
    //cerr <<"coeff. riavv. interasse = "<< K << endl;

    //smusso di testa
    D31=s_1; //s per pignone 
    E31=s_2; //s per ruota

    // Diametro fine c.to testa max
    da1_=dp1_+2*mn*(cd1_+x1_)-2*K*mn-2*D31;
    da2_=dp2_+2*mn*(cd2_+x2_)-2*K*mn-2*E31;
    //cerr <<"Diametro fine c.to testa max 1 = "<< da1_ << " Diametro fine c.to testa max 2 = "<< da2_ << endl;
    cerr <<" ------------------------------------- "<< endl;
    cerr <<" diametro fine c.to testa = "<< da1_ << endl;
    // diametro di piede
    df1_=dp1_-2*mn*(ca1_-xg1);
    df2_=dp2_-2*mn*(ca2_-xg2);
    //cerr <<"diametro di piede 1 = "<< df1_ << " diametro di piedex 2 = "<< df2_ << endl;

    // Diam. inizio ev.  teo.
    dtif1=2*sqrt(db1_*db1_/4+(db1_/2*tan(at/180*pi)-(mn*(ca1_-rt1/mn*(1-sin(an/180*pi))-xg1))/sin(at/180*pi))*(db1_/2*tan(at/180*pi)-(mn*(ca1_-rt1/mn*(1-sin(an/180*pi))-xg1))/sin(at/180*pi)));
    dtif2=2*sqrt(db2_*db2_/4+(db2_/2*tan(at/180*pi)-(mn*(ca1_-rt2/mn*(1-sin(an/180*pi))-xg2))/sin(at/180*pi))*(db2_/2*tan(at/180*pi)-(mn*(ca1_-rt2/mn*(1-sin(an/180*pi))-xg2))/sin(at/180*pi)));
    //cerr <<"Diam. inizio ev.  teo. 1 = "<< dtif1 << " Diam. inizio ev.  teo. 2 = "<< dtif2 << endl;

    // diametro primitivo di funzionamento
    dpdf1=2*adf/(u+1);
    dpdf2=2*adf/(u+1);
    //cerr <<"diametro primitivo di funzionamento 1 = "<< dpdf1 << " diametro primitivo di funzionamento 2 = "<< dpdf2 << endl;

    // angolo d'elica di funzionamento
    bdf=atan(tan(b/180*pi)*dpdf1/dp1_)/pi*180;
    //cerr <<"angolo d'elica di funzionamento = "<< bdf << endl;

    // A. di p. norm. di funz.
    andf=atan(tan(atdf/180*pi)*cos(atan(tan(b/180*pi)*dpdf1/dp1_)))/pi*180;
    //cerr <<"A. di p. norm. di funz. = "<< andf << endl;

    // Sp. trasv. al prim. di funz.
    stdf1=dpdf1*(((pi/2+2*xg1*tan(an/180*pi)))/z1_-(tan(atdf/180*pi)-(atdf/180*pi))+(tan(at/180*pi)-(at/180*pi)));
    stdf2=dpdf2*(((pi/2+2*xg2*tan(an/180*pi)))/z2_-(tan(atdf/180*pi)-(atdf/180*pi))+(tan(at/180*pi)-(at/180*pi)));
    //cerr <<"Sp. trasv. al prim. di funz. 1 = "<< stdf1 << " Sp. trasv. al prim. di funz. 2 = "<< stdf2 << endl;

    // Add. di funz.
    adf1=(da1_-dpdf1)/2;
    adf2=(da2_-dpdf2)/2;
    //cerr <<"Add. di funz. 1 = "<< adf1 << " Add. di funz. 2 = "<< adf2 << endl;

    // Ded. di funz.
    ddf1=-(df1_-dpdf1)/2;
    ddf2=-(df2_-dpdf2)/2;
    //cerr <<"Ded. di funz. 1 = "<< ddf1 << " Ded. di funz. 2 = "<< ddf2 << endl;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // PUNTI DEL PROFILO

    // PARTE ATTIVA pignone  ---------------------------------------------------------------------------------------------------
	
    // sfasamento ev. pignone - asse vano
    q0p=(pi*mt/2-2*xg1*mn*tan(at/180*pi))/2/(dp1_/2)-(tan(at/180*pi)-(at/180*pi));

    // sfasamento ev. ruota - asse vano
    q0r=(pi*mt/2-2*xg2*mn*tan(at/180*pi))/2/(dp1_/2)-(tan(at/180*pi)-(at/180*pi));

    // definisco q per il primo punto
    q=sqrt((dtif1/db1_)*(dtif1/db1_)-1);
    //cerr <<"q = "<< q << endl;

    // Ascissa attiva simm dente fianco sx
    x1sx=(-rb1_*sin(q+q0p)+rb1_*q*cos(q+q0p))*cos(-pi/z1_)-(rb1_*cos(q+q0p)+rb1_*q*sin(q+q0p))*sin(-pi/z1_);
    //cerr <<"x1 = "<< x1 << endl;

    // Ascissa attiva simm dente fianco sx
    x1dx=-x1sx;

    // Ordinata attiva simm dente fianco sx
    y1sx=(-rb1_*sin(q+q0p)+rb1_*q*cos(q+q0p))*sin(-pi/z1_)+(rb1_*cos(q+q0p)+rb1_*q*sin(q+q0p))*cos(-pi/z1_);
    //cerr <<"y1 = "<< y1 << endl;

    // Ascissa attiva simm dente fianco sx
    y1dx=y1sx;

    //outputFile << q << " " << x1sx << " " << y1sx <<"\n";
    //outputFile << q << " " << x1dx << " " << y1dx <<"\n";

    outputFile << x1sx << "\t\t" << y1sx <<"\n";
    punti[nPoints][0]=x1sx;
    punti[nPoints][1]=y1sx;

    nPoints=nPoints+1;
    outputFile << x1dx << "\t" << y1dx <<"\n";
    punti[nPoints][0]=x1dx;
    punti[nPoints][1]=y1dx;

    nPoints=nPoints+1;
    for (i=1; i<30;i++)
    {
        // aggiorno q per il prossimo punto
        q=q+(sqrt((da1_/db1_)*(da1_/db1_)-1)-sqrt((dtif1/db1_)*(dtif1/db1_)-1))/29; 
        // Ascissa attiva simm dente fianco sx
        x1sx=(-rb1_*sin(q+q0p)+rb1_*q*cos(q+q0p))*cos(-pi/z1_)-(rb1_*cos(q+q0p)+rb1_*q*sin(q+q0p))*sin(-pi/z1_);
        // Ascissa attiva simm dente fianco sx
        x1dx=-x1sx;
        // Ordinata attiva simm dente fianco sx
        y1sx=(-rb1_*sin(q+q0p)+rb1_*q*cos(q+q0p))*sin(-pi/z1_)+(rb1_*cos(q+q0p)+rb1_*q*sin(q+q0p))*cos(-pi/z1_);
        // Ascissa attiva simm dente fianco sx
        y1dx=y1sx;

        //outputFile << q << " " << x1sx << " " << y1sx <<"\n";
        //outputFile << q << " " << x1dx << " " << y1dx <<"\n";
        outputFile << x1sx << "\t" << y1sx <<"\n";
        punti[nPoints][0]=x1sx;
        punti[nPoints][1]=y1sx;

        nPoints=nPoints+1;
        outputFile << x1dx << "\t" << y1dx <<"\n";
        punti[nPoints][0]=x1dx;
        punti[nPoints][1]=y1dx;
        nPoints=nPoints+1;

    }
    alpha=2*atan(x1sx/y1sx);
    // cerr <<"alpha = "<< alpha << endl;
    // FINE PARTE ATTIVA pignone  -----------------------------------------------------------------------------------------------
    // TRONCATURA ESTERNA pignone  ---------------------------------------------------------------------------------------------------

    q=0;
    x1sx=0;
    y1sx=da1_/2;
    i=0;
    outputFile << x1sx << "\t" << y1sx <<"\n";
    punti[nPoints][0]=x1sx;
    punti[nPoints][1]=y1sx;
    nPoints=nPoints+1;
    // cerr << x1sx << "\t" << y1sx << endl;
    for (i=1; i<10;i++)
    {
	q=q+(alpha/2)/9;
	x1sx=-(da1_/2)*sin(q);
	y1sx=da1_/2-da1_/2*(1-cos(q));
	outputFile << x1sx << "\t" << y1sx <<"\n";
	punti[nPoints][0]=x1sx;
	punti[nPoints][1]=y1sx;
	nPoints=nPoints+1;
	x1dx=+(da1_/2)*sin(q);
	y1dx=da1_/2-da1_/2*(1-cos(q));
	outputFile << x1dx << "\t" << y1dx <<"\n";
	punti[nPoints][0]=x1dx;
	punti[nPoints][1]=y1dx;
	nPoints=nPoints+1;
    }
    // FINE TRONCATURA ESTERNA pignone ------------------------------------------------------------------------------------------------

    // RACCORDO pignone  ---------------------------------------------------------------------------------------------------

    u1min=-(dp1_/2+xg1*mn-rt1-0.5*df1_+1/(1-sin(an/180*pi))*(dp1_-df1_)/2*sin(an/180*pi))*((1-sin(an/180*pi))/cos(an/180*pi));
    //cerr <<"u1min = "<< u1min << endl;
    qrmax1=acos(sin(an/180*pi));
    //cerr <<"qrmax1 = "<< qrmax1 << endl;
    lcx1= (cos(an/180*pi)/(1-sin(an/180*pi))*(dp1_-df1_)/2+u1min)/cos(b/180*pi);
    //cerr <<"lcx1 = "<< lcx1 << endl;
    lcy1=dp1_/2+xg1*mn-rt1-df1_/2;
    //cerr <<"lcy1 = "<< lcy1 << endl;
    i=1;
    q=0;
    x1dx=xrs(q,u1min,lcx1,lcy1,rt1,db1_/2,dp1_/2,q0p,b,mn*xg1,qrmax1,sqrt((dtif1/db1_)*(dtif1/db1_)-1),z1_);
    y1dx=yrs(q,u1min,lcx1,lcy1,rt1,db1_/2,dp1_/2,q0p,b,mn*xg1,qrmax1,sqrt((dtif1/db1_)*(dtif1/db1_)-1),z1_);
    x1sx=-x1dx;
    y1sx=y1dx;
    outputFile << x1sx << "\t" << y1sx <<"\n";
    punti[nPoints][0]=x1sx;
    punti[nPoints][1]=y1sx;
    nPoints=nPoints+1;
    outputFile << x1dx << "\t" << y1dx <<"\n";
    punti[nPoints][0]=x1dx;
    punti[nPoints][1]=y1dx;
    nPoints=nPoints+1;
    beta=atan(-x1dx/y1dx);
    //cerr <<"beta = "<< beta*360/2/pi << endl;
    // cerr <<"x1dx = "<< x1dx << endl;
    // cerr <<"y1dx = "<< y1dx << endl;
    for (i=1; i<30;i++)
    {
	q=q+qrmax1/29;
	x1dx=xrs(q,u1min,lcx1,lcy1,rt1,db1_/2,dp1_/2,q0p,b,mn*xg1,qrmax1,sqrt((dtif1/db1_)*(dtif1/db1_)-1),z1_);
	y1dx=yrs(q,u1min,lcx1,lcy1,rt1,db1_/2,dp1_/2,q0p,b,mn*xg1,qrmax1,sqrt((dtif1/db1_)*(dtif1/db1_)-1),z1_);
	x1sx=-x1dx;
	y1sx=y1dx;
	outputFile << x1sx << "\t" << y1sx <<"\n";
	punti[nPoints][0]=x1sx;
	punti[nPoints][1]=y1sx;
	nPoints=nPoints+1;
	outputFile << x1dx << "\t" << y1dx <<"\n";
	punti[nPoints][0]=x1dx;
	punti[nPoints][1]=y1dx;
	nPoints=nPoints+1;
    }
    // FINE RACCORDO pignone  ---------------------------------------------------------------------------------------------------


    // TRONCATURA INTERNA pignone  ---------------------------------------------------------------------------------------------------

    q=0;
    i=1;
    gamma=pi/z1_;
    //cerr <<"gamma = "<< gamma*360/2/pi << endl;

    for (i=1; i<8;i++)
    {
	q=q+(gamma-beta)/6;
	x1sx=-df1_/2*sin(beta+q);
        // cerr <<"x = "<< x1sx << endl;
	y1sx=df1_/2-df1_/2*(1-cos(beta+q));
        // cerr <<"y = "<< y1sx << endl;
	x1dx=-x1sx;
	y1dx=y1sx;
	outputFile << x1sx << "\t" << y1sx <<"\n";
	punti[nPoints][0]=x1sx;
	punti[nPoints][1]=y1sx;
	nPoints=nPoints+1;
        // cerr <<"q = "<< q*360/2/pi << endl;

	outputFile << x1dx << "\t" << y1dx <<"\n";
	punti[nPoints][0]=x1dx;
	punti[nPoints][1]=y1dx;
	nPoints=nPoints+1;

		
    }

    //cerr <<"nPoints = "<< nPoints << endl;
    // FINE TRONCATURA INTERNA pignone ------------------------------------------------------------------------------------------------

    // SWEEP

    double point[nPoints][2];
    int nPunti;

    nPunti=0;
	
    //ifstream dente;
    //dente.open("puntibis.dat");
    i=0;
    /*	
        while(!dente.eof()) 
        {
		
        dente >> point[i][1] >> point[i][2];
        cerr << point[i][1] << " " << point[i][2] <<"\n";
        i++;
        }
	dente.close();




	for (i=0; i<nPoints; i++)
        {
        dente >> point[i][1] >> point[i][2];
        //cerr << point[i][1] << " " << point[i][2] <<"\n";
        }
    */
    for (t=0; t<z1_; t++) 
    {

        i=0;


        dAngle=2*pi*t/z1_;
        //cerr <<"dAngle = "<< dAngle << endl;

        for (i=0; i<nPoints; i++)
        {
            x = punti[i][0];
            y = punti[i][1];
            //r = sqrt(x*x + y*y);
            //angle = atan2(y,x) + dAngle;
            xmod = x*cos(dAngle) - y*sin(dAngle);
            ymod = x*sin(dAngle) + y*cos(dAngle);

            //cerr <<"angle = "<< angle << endl;
            //xmod = r * cos(angle);
            //ymod = r * sin(angle);
            outputFile2 << xmod << "\t" << ymod <<"\n";
            nPunti++;
        }

    }
    ofstream numPunti;
    numPunti.open("n.dat");
    numPunti << nPunti <<"\n";
    cerr <<" n punti = "<< nPunti << endl;
    cerr <<" ---------------------------------------- "<< endl;



    /////////////////////////////////////////////////////////////////////////////////////

    //Generazione file per pointTOstl

    ofstream solids;
    solids.open("solids");
    solids <<"1 0 "<<l1_<< endl;	
    solids <<"profilo.dat "<< nPunti << " 0 0 0" << endl;

    return 0;
}

