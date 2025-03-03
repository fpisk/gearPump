/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "argList.H"
#include "Time.H"
#include "IOdictionary.H"
#include "OFstream.H"
#include "gearProfile.H"

#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



scalar xr
(
    scalar qr,
    scalar umin,
    scalar lx,
    scalar ly,
    scalar rt,
    scalar rb,
    scalar rp,
    scalar q0,
    scalar u_b,
    scalar V,
    scalar qrmax,
    scalar qamin
)
{


    const scalar pi = constant::mathematical::pi;

    // scalar firQMax, XrVQMax, YrVQMax, RotR, XaQMax, RotA, RotRA, fir, xr, xrv, yrv, YaQMax;

    // const scalar firQMax = 1.0/rp*((-lx+rt*sin(qrmax))/cos(u_b/180*pi)+(V-ly-rt*cos(qrmax))*tan(qrmax)*cos(u_b/180*pi));

    // XrVQMax=cos(firQMax)*(-lx+rt*sin(qrmax))/cos(u_b/180*pi)+sin(firQMax)*(-ly-rt*cos(qrmax))+sin(firQMax)*rp+sin(firQMax)*V-rp*firQMax*cos(firQMax);

    // YrVQMax=-sin(firQMax)*(-lx+rt*sin(qrmax))/cos(u_b/180*pi)+cos(firQMax)*(-ly-rt*cos(qrmax))+cos(firQMax)*rp+cos(firQMax)*V+rp*firQMax *sin(firQMax);

    // RotR=(abs(YrVQMax)/(YrVQMax))*atan(abs(YrVQMax/XrVQMax));

    // if (XrVQMax<0) RotR=pi-RotR;

    // XaQMax=rb*sin(qamin+q0)-rb*qamin*cos(qamin+q0);

    // YaQMax=rb*cos(qamin+q0)+rb*qamin*sin(qamin+q0);


    // RotA= (abs(YaQMax)/(YaQMax))*atan(abs(YaQMax/XaQMax));

    // if (XaQMax<0) RotA=pi-RotA;

    // RotRA=RotA-RotR;

    // fir=1/rp*((-lx+rt*sin(qr))/cos(u_b/180*pi)+(V-ly-rt*cos(qr))*tan(qr)*cos(u_b/180*pi));

    // xrv=cos(fir)*(-lx+rt*sin(qr))/cos(u_b/180*pi)+sin(fir)*(-ly-rt*cos(qr))+sin(fir)*rp+sin(fir)*V-rp*fir*cos(fir);

    // yrv=-sin(fir)*(-lx+rt*sin(qr))/cos(u_b/180*pi)+cos(fir)*(-ly-rt*cos(qr))+cos(fir)*rp+cos(fir)*V+rp*fir*sin(fir);

    // xr=xrv*cos(RotRA)-yrv*sin(RotRA);


    // return xr;
    return 0;
}


// Main program:
int main(int argc, char *argv[])
{
    argList::noParallel();

    argList::addOption
    (
        "dict",
        "file",
        "specify alternative dictionary for the blockMesh description"
    );

    argList args(argc, argv, false, true);

    Time runTime(args.rootPath(), args.caseName());

    word dictName("gearsGenDict");

    // Search for the appropriate blockMesh dictionary....
    fileName dictPath = dictName;

    // Check if the dictionary is specified on the command-line
    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];

        dictPath =
        (
            isDir(dictPath)
          ? dictPath/dictName
          : dictPath
        );
    }


    const gearProfile gear
    (
        IOdictionary
        (
            IOobject
            (
                dictPath,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        )
    );


    // write points
    fileName fName("points");
    OFstream pts(fName);;

    // profile
    fName = "profile";
    OFstream profile(fName);;

    // solids
    fName = "solids";
    OFstream solids(fName);


    Info<< nl
        << "--------------------------" << nl
        << "Addressing" << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
