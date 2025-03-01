/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gearProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//const dataType Foam::gearProfile::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::gearProfile::createProfile()
{
    const label constRatio = 0.3;
    
    // transmission ratio
    const scalar u = transmissionRatio_;

    // Raggio testa creatore
    const scalar rt = constRatio * mn_;
    //rt2 = rt

    // Modulo trasv. teorico
    const scalar mt = mn_ / cos(degToRad(b_));

    // angolo di pressione trasv. teorico
    const scalar at = degToRad(tan(an_)/cos(b_));

    // diam. primario teorico
    const scalar dp1_ = mt*z1_;

    // z2 - number of teeth of gear B
    const label z2 = transmissionRatio_*z1_;
    
    const scalar dp2_ = mt*z2;

    // interasse teorico
    const scalar a = (dp1_ + dp2_)/2.0;
    
    // angolo d'elica base [rad]
    const scalar helixAngleBase = degToRad
    (
        asin(sin(b_)*cos(an_))
    );
    const scalar bb = helixAngleBase;

    // Coef. sp.to con gioco
    const scalar xg1 = x1_-0.5*jn1/mn_/sin(an_);
    
    // Numero denti virtuale
    const scalar nVirTeethA = z1_/(sqr(cos(bb))/cos(b_));

    scalar zn1 = nVirTeethA;

    scalar zn2 = z2/(sqr(cos(bb))/cos(b_));
    
    //cerr <<"Numero denti virtuale 1 = "<< zn1 << " Numero denti virtuale 2 = "<< zn2 << endl;
    /*
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

    */
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gearProfile::gearProfile(const dictionary& dict)
:
z1_(dict.lookup<scalar>("numberTeeth")),
transmissionRatio_(dict.lookupOrDefault<scalar>("transmissionRatio", 1)),
mn_(dict.lookup<scalar>("normalModule")),
an_(degToRad(dict.lookup<scalar>("normalPressureAngle"))),
b_(degToRad(dict.lookup<scalar>("helixAngle"))),
jn1(dict.lookup<scalar>("normalBacklashContribution")),
ca1(dict.lookup<scalar>("addendumCoeffCutter")),
cd1(dict.lookup<scalar>("dedendumCoeffCutter")),
l1_ (dict.lookup<scalar>("effectiveFaceWidth")),
adf(dict.lookup<scalar>("operatingCenterDistance")),
x1_(dict.lookup<scalar>("coeffUnderZeroBacklash")),
s_1(dict.lookup<scalar>("tipChamfer"))
{
    createProfile();
}





// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

//Foam::autoPtr<Foam::gearProfile>
//Foam::gearProfile::New()
//{
//    return autoPtr<gearProfile>(new gearProfile);
//}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
