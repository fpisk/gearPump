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

Class
    Foam::gearProfile

Description

SourceFiles
    gearProfile.C

\*---------------------------------------------------------------------------*/

#ifndef gearProfile_H
#define gearProfile_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
// class gearProfile;
// Istream& operator>>(Istream&, gearProfile&);
// Ostream& operator<<(Ostream&, const gearProfile&);


/*---------------------------------------------------------------------------*\
                          Class gearProfile Declaration
\*---------------------------------------------------------------------------*/

class gearProfile
{
    // Private Data

        //- numberTeeth
        const label z1_; 
        const label z2_; 
                
        // normalModule
	const label mn_; 
	// normalPressureAngle [deg]
	const scalar an;

	// helixAngle [deg]
	const scalar b;   

        // normalBacklashContribution [deg] 
        const scalar jn1; 
        
        // addendumCoeffCutter 
        const scalar ca1; 
       
        // dedendumCoeffCutter  
        const scalar cd1; 
        
        // effectiveFaceWidth [mm] 
        const scalar l1_ ; 
       
        // operatingCenterDistance [mm] 
        const scalar adf; 
        
        // coeffUnderZeroBacklash [-] 
        const scalar x1_; 

        // tipChamfer
        const scalar s_1;  
  

    // Private Member Functions


public:

    // Static Data Members

        //- Static data staticData
 


    // Constructors

        //- Construct from components
        gearProfile(const dictionary& dict);

        //- Construct from Istream
        gearProfile(Istream&);

        //- Copy constructor
        // gearProfile(const gearProfile&);

        //- Disallow default bitwise copy construction
        // gearProfile(const gearProfile&) = delete;


    // Selectors

        //- Select null constructed
        //  static autoPtr<gearProfile> New();


    //- Destructor
    ~gearProfile() {};


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

       //  void operator=(const gearProfile&);

        //- Disallow default bitwise assignment
        // void operator=(const gearProfile&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

       // friend Istream& operator>>(Istream&, gearProfile&);
       // friend Ostream& operator<<(Ostream&, const gearProfile&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "gearProfileI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
