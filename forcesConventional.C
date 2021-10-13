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

    //return move(names); // OF9
    return names; // OF2106
}

void Foam::functionObjects::forcesConventional::writeFileHeader(const label i)
{
    Info<< "writeFileHeader" << nl << endl;
    switch (fileID(i))
    {
        case fileID::mainFile:
        {
            // Force data
            //writeHeader(file(i), "Forces");
            if (!porosity_)
            {   
                //file(i) // OF9
                files(i) // OF2106
                << "time,"
                << "pressure_fx" << "," << "pressure_fy" << "," << "pressure_fz" << ","
                << "viscous_fx" << "," << "viscous_fy" << "," << "viscous_fz" << ","
                << "pressure_mx" << "," << "pressure_my" << "," << "pressure_mz" << ","
                << "viscous_mx" << "," << "viscous_my" << "," << "viscous_mz";
            }
            else
            {
                //file(i) // OF9
                files(i) // OF2106
                << "time,"
                << "pressure_fx" << "," << "pressure_fy" << "," << "pressure_fz" << ","
                << "viscous_fx" << "," << "viscous_fy" << "," << "viscous_fz" << ","
                << "darcy_fx" << "," << "darcy_fy" << "," << "darcy_fz" << ","
                << "forchheimer_fx" << "," << "forchheimer_fy" << "," << "forchheimer_fz" << ","
                << "pressure_mx" << "," << "pressure_my" << "," << "pressure_mz" << ","
                << "viscous_mx" << "," << "viscous_my" << "," << "viscous_mz" << ","
                << "darcy_mx" << "," << "darcy_my" << "," << "darcy_mz" << ","
                << "forchheimer_mx" << "," << "forchheimer_my" << "," << "forchheimer_mz";
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled file index: " << i
                << abort(FatalError);
        }
    }

    //file(i) << endl; // OF9
    files(i) << endl; // OF2106
}

Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::forcesConventional::devTau() const
{
    const dictionary& transportProperties = obr_.lookupObject<dictionary>("transportProperties");

    dimensionedScalar nu
    (
        "nu",
        dimViscosity,
        transportProperties // OF2106
        //transportProperties.lookup("nu") // OF9
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
    porousZoneSet_(),
    pName_(word::null),
    UName_(word::null),
    K_Name_(word::null),
    //rhoRef_(dict.lookup<scalar>("rho")), // OF9
    rhoRef_(dict.get<scalar>("rho")), // OF2106
    pRef_(0),
    porosity_(false),
    forceP_(),
    forceV_(),
    forceD_(),
    forceF_(),
    momentP_(),
    momentV_(),
    momentD_(),
    momentF_()
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

    // Get cell zone mesh on porousZone
    if (porosity_)
    {
        //const label cellZoneID = mesh_.cellZones().findZoneID(dict.lookup("porousZone")); // OF9
        const label cellZoneID = mesh_.cellZones().findZoneID("porousZone"); // OF2106

        const cellZoneMesh& zoneMesh = mesh_.cellZones()[cellZoneID].zoneMesh();

        porousZoneSet_ = zoneMesh[cellZoneID];
    }

    // Get p and U field names
    pName_ = dict.lookupOrDefault<word>("p", "p");

    UName_ = dict.lookupOrDefault<word>("U", "U");

    K_Name_ = dict.lookupOrDefault<word>("K_", "K_");

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
    forceD_ = Zero;
    forceF_ = Zero;

    momentP_ = Zero;
    momentV_ = Zero;
    momentD_ = Zero;
    momentF_ = Zero;

    // Get mesh boundary field
    const surfaceVectorField::Boundary& Sfb = mesh_.Sf().boundaryField();

    // Get pressure field
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    // Get velocity field
    const volVectorField& U = obr_.lookupObject<volVectorField>(UName_);

    // Get viscous stress tensor field
    tmp<volSymmTensorField> tdevTau = devTau();
    const volSymmTensorField::Boundary& devTaub = tdevTau().boundaryField();

    // Scale pRef by density for incompressible simulations
    scalar pRef = pRef_/rhoRef_;

    forAllConstIter(labelHashSet, patchSet_, iter)
    {
        // Calculate pressure and viscous force/moment contributions
        Info << "Calculating pressure and viscous force/moment contributions." << endl;

        label patchi = iter.key();

        // Get boundary field positions
        vectorField Md
        (
            mesh_.C().boundaryField()[patchi]
            //mesh_.C().boundaryField()[patchi] - coordSys_.origin()
        );

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

    if (porosity_)
    {
        // Calculate porous force/moment contributions
        Info << "Calculating porous force/moment contributions." << endl;

        // Get permeability field
        const volScalarField& K_ = obr_.lookupObject<volScalarField>(K_Name_);
        
        // Get nu and cf
        const dictionary& transportProperties = obr_.lookupObject<dictionary>("transportProperties");
        
        dimensionedScalar nu
        (
            "nu",
            dimViscosity,
            transportProperties // OF2106
            //transportProperties.lookup("nu") // OF9
        );

        dimensionedScalar cf
        (
            "cf",
            dimless,
            transportProperties // OF2106
            //transportProperties.lookup("cf") // OF9
        );

        forAll(porousZoneSet_, i)
        {
            // Get Darcy contribution from ith cell in porousZoneSet_
            vector fD
            (
                rhoRef_
                *mesh_.V()[porousZoneSet_[i]]
                *((nu.value()/K_[porousZoneSet_[i]])*U[porousZoneSet_[i]])
            );

            // Get Forchheimer contribution from ith cell in porousZoneSet_
            vector fF
            (
                rhoRef_
                *mesh_.V()[porousZoneSet_[i]]
                *((cf.value()/sqrt(K_[porousZoneSet_[i]]))*mag(U[porousZoneSet_[i]])*U[porousZoneSet_[i]])
            );

            forceD_ += fD;
            forceF_ += fF;

            momentD_ += mesh_.C()[porousZoneSet_[i]]^fD;
            momentF_ += mesh_.C()[porousZoneSet_[i]]^fF;
        }

        // Parallel
        Pstream::combineGather(forceD_, plusEqOp<vector>());
        Pstream::combineGather(forceF_, plusEqOp<vector>());
        Pstream::combineGather(momentD_, plusEqOp<vector>());
        Pstream::combineGather(momentF_, plusEqOp<vector>());
        Pstream::combineScatter(forceD_);
        Pstream::combineScatter(forceF_);
        Pstream::combineScatter(momentD_);
        Pstream::combineScatter(momentF_);
    }

    // Parallel
    Pstream::combineGather(forceP_, plusEqOp<vector>());
    Pstream::combineGather(forceV_, plusEqOp<vector>());
    Pstream::combineGather(momentP_, plusEqOp<vector>());
    Pstream::combineGather(momentV_, plusEqOp<vector>());
    Pstream::combineScatter(forceP_);
    Pstream::combineScatter(forceV_);
    Pstream::combineScatter(momentP_);
    Pstream::combineScatter(momentV_);

    //Info<< "forceP_ = " << forceP_ << nl << forceP_.x() << endl;
    //Info<< "forceF_ = " << forceF_ << endl;

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

        //writeTime(file(fileID::mainFile)); // OF9
        writeCurrentTime(file(fileID::mainFile)); // OF2106

        if (!porosity_)
        {
            file(fileID::mainFile) << ","
            << forceP_.x() << "," << forceP_.y() << "," << forceP_.z() << ","
            << forceV_.x() << "," << forceV_.y() << "," << forceV_.z() << ","
            << momentP_.x() << "," << momentP_.y() << "," << momentP_.z() << ","
            << momentV_.x() << "," << momentV_.y() << "," << momentV_.z() << endl;
        }
        else
        {
            file(fileID::mainFile) << ","
            << forceP_.x() << "," << forceP_.y() << "," << forceP_.z() << ","
            << forceV_.x() << "," << forceV_.y() << "," << forceV_.z() << ","
            << forceD_.x() << "," << forceD_.y() << "," << forceD_.z() << ","
            << forceF_.x() << "," << forceF_.y() << "," << forceF_.z() << ","
            << momentP_.x() << "," << momentP_.y() << "," << momentP_.z() << ","
            << momentV_.x() << "," << momentV_.y() << "," << momentV_.z() << ","
            << momentD_.x() << "," << momentD_.y() << "," << momentD_.z() << ","
            << momentF_.x() << "," << momentF_.y() << "," << momentF_.z() << endl;
        }
    }
    
    return true;
}


// ************************************************************************* //
