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

\*---------------------------------------------------------------------------*/

#include "forcesConventional.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Time.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(forcesConventional, 0);
    addToRunTimeSelectionTable(functionObject, forcesConventional, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::functionObjects::forcesConventional::createFileNames
(
    const dictionary& dict
) const
{
    Info<< "createFileNames" << nl << endl;
    DynamicList<word> names(1);

    const word forceType(dict.lookup("type"));

    // Name for file(fileID::mainFile=0)
    names.append(forceType);

    return move(names);
}

void Foam::functionObjects::forcesConventional::writeFileHeader(const label i)
{
    Info<< "writeFileHeader" << nl << endl;
    switch (fileID(i))
    {
        case fileID::mainFile:
        {
            // Force data
            writeHeader(file(i), "Forces");

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }

    file(i) << endl;
}

Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forcesConventional::devTau() const
{
    const dictionary& transportProperties = obr_.lookupObject<dictionary>("transportProperties");

    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        transportProperties.lookup("nu")
    );

    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

    return -rhoRef_*nu*dev(twoSymm(fvc::grad(U)));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::forcesConventional::forcesConventional
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    patchSet_(),
    pName_(word::null),
    UName_(word::null),
    rhoRef_(dict.lookup<scalar>("rho")),
    pRef_(0),
    porosity_(false),
    forceP_(),
    forceV_(),
    momentP_(),
    momentV_()
{
    read(dict); // Read dict data
    resetNames(createFileNames(dict)); // Setup files for saving data
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::forcesConventional::~forcesConventional()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::forcesConventional::read(const dictionary& dict)
{
    //----- Run on first timestep -----//

    // Get boundary mesh on patchSet_
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ = pbm.patchSet(wordReList(dict.lookup("patches")));

    // Get p and U field names
    pName_ = dict.lookupOrDefault<word>("p", "p");

    UName_ = dict.lookupOrDefault<word>("U", "U");

    // Get density
    dict.lookup("rho") >> rhoRef_;

    // Reference pressure, 0 by default
    pRef_ = dict.lookupOrDefault<scalar>("pRef", 0.0);

    // Calculate for a porous body?
    dict.readIfPresent("porosity", porosity_);

    return true;
}

bool Foam::functionObjects::forcesConventional::execute()
{
    //----- Run every timestep -----//

    // Set force/moment vectors = Zero
    forceP_ = Zero;
    forceV_ = Zero;

    momentP_ = Zero;
    momentV_ = Zero;

    // Get mesh boundary field
    const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

    // Get pressure field
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    // Get viscous stress tensor field
    tmp<volSymmTensorField> tdevTau = devTau();
    const volSymmTensorField::Boundary& devTaub = tdevTau().boundaryField();

    // Scale pRef by density for incompressible simulations
    scalar pRef = pRef_/rhoRef_;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        label patchi = iter.key();

        // Get boundary field positions
        vectorField Md
        (
            mesh_.C().boundaryField()[patchi]
            //mesh_.C().boundaryField()[patchi] - coordSys_.origin()
        );

        if (!porosity_)
        {
            // Calculate forces/moments for a solid body
            Info << "Calculating forces/moments for a solid body." << endl;

            // Get normal force vector field (pressure)
            vectorField fN
            (
                rhoRef_*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
            );

            // Get tangential force vector field (viscous)
            vectorField fT
            (
                Sfb[patchi] & devTaub[patchi]
            );

            forceP_ += sum(fN);
            forceV_ += sum(fT);

            momentP_ += sum(Md^fN);
            momentV_ += sum(Md^fT);
        }
        else
        {
            // Calculate forces/moments for a porous body
            Info << "Calculating forces/moments for a porous body." << endl;

            // Get normal force vector field (pressure)
            vectorField fN
            (
                rhoRef_*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
            );

            // Get tangential force vector field (viscous)
            vectorField fT
            (
                Sfb[patchi] & devTaub[patchi]
            );

            forceP_ += sum(fN);
            forceV_ += sum(fT);

            momentP_ += sum(Md^fN);
            momentV_ += sum(Md^fT);
        }
    }

    //Info<< "forceP_ = " << forceP_ << nl << forceV_ << endl;

    return true;
}

bool Foam::functionObjects::forcesConventional::end()
{
    //----- Run on final timestep -----//

    return true;
}

bool Foam::functionObjects::forcesConventional::write()
{
    //----- Run every timestep -----//

    // Write out data to file
    if (Pstream::master())
    {
        logFiles::write();

        writeTime(file(fileID::mainFile));

        file(fileID::mainFile) << tab
            << forceP_ << tab
            << forceV_ << tab
            << momentP_ << tab
            << momentV_ << endl;
    }
    
    return true;
}


// ************************************************************************* //
