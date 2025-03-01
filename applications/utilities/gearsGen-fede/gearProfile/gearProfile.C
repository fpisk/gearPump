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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::gearProfile::gearProfile(const dictionary& dict)
:
z1_(dict.lookup<scalar>("numberTeeth")),
z2_(z1_), 
mn_(dict.lookup<scalar>("normalModule")),
an(dict.lookup<scalar>("normalPressureAngle")), 
b(dict.lookup<scalar>("helixAngle")), 
jn1(dict.lookup<scalar>("normalBacklashContribution")),
ca1(dict.lookup<scalar>("addendumCoeffCutter")),  
cd1(dict.lookup<scalar>("dedendumCoeffCutter")), 
l1_ (dict.lookup<scalar>("effectiveFaceWidth")), 
adf(dict.lookup<scalar>("operatingCenterDistance")),
x1_(dict.lookup<scalar>("coeffUnderZeroBacklash")), 
s_1(dict.lookup<scalar>("tipChamfer"))
{


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
