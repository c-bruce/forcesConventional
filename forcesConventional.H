/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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
    Foam::functionObjects::FUNCTIONOBJECT

Description
    This function object...

    Example of function object specification:
    \verbatim
    FUNCTIONOBJECT1
    {
        type           FUNCTIONOBJECT;
        libs ("libFUNCTIONOBJECTFunctionObject.so");
        ...
        wordData       someWord;
        scalarData     1.0;
        labelData      1;
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name: FUNCTIONOBJECT | yes      |
        wordData     | some word option...       | no       | defaultWord
        scalarData   | some scalar value...      | yes      |
        labelData    | some label value...       | yes      |
    \endtable

SourceFiles
    FUNCTIONOBJECT.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_forcesConventional_H
#define functionObjects_forcesConventional_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                   Class FUNCTIONOBJECT Declaration
\*---------------------------------------------------------------------------*/

class forcesConventional
:
    public fvMeshFunctionObject,
    public logFiles
{

private:

    // Private Data

        //- Enumeration for ensuring the right file is accessed
        enum class fileID
        {
            mainFile = 0//,
            //binsFile = 1
        };

        // Read from dictionary

            //- Patches to integrate forces over
            labelHashSet patchSet_;

            //- Name of pressure field
            word pName_;

            //- Name of velocity field
            word UName_;

            //- Reference density for incompressible calculations
            scalar rhoRef_;

            //- Reference pressure
            scalar pRef_;
        
        //Forces

            //- Pressure
            vector forceP_;

            //- Viscous
            vector forceV_;

        // Moments

            //- Pressure
            vector momentP_;

            //- Viscous
            vector momentV_;

    // Private Member Functions

        using logFiles::file;

        Ostream& file(const fileID fid)
        {
            return logFiles::file(label(fid));
        }

        //- Create file names
        wordList createFileNames(const dictionary& dict) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Return the effective viscous stress
        tmp<volSymmTensorField> devTau() const;

public:

    //- Runtime type information
    TypeName("forcesConventional");


    // Constructors

        //- Construct from Time and dictionary
        forcesConventional
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        forcesConventional(const forcesConventional&) = delete;


    //- Destructor
    virtual ~forcesConventional();


    // Member Functions

        //- Read the forcesConventional data
        virtual bool read(const dictionary&);

        //- Execute, currently does nothing
        virtual bool execute();

        //- Execute at the final time-loop, currently does nothing
        virtual bool end();

        //- Write the forcesConventional
        virtual bool write();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const forcesConventional&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
