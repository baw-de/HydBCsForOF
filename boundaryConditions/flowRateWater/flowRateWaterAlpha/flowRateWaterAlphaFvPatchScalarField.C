/*---------------------------------------------------------------------------*\
 * 
 * Author
 *   Carsten Thorenz 2024, Federal Waterways Research Institute 
 * 
 * Publication
 *   Thorenz, C. (2024): 'Boundary Conditions for Hydraulic Structures Modelling with OpenFOAM',
 *   10th International Symposium on Hydraulic Structures, Zürich. ISSN 0374-0056 , DOI: 10.3929/ethz-b-000675949
 *
 * License
 *    GPL v3
 * 
\*---------------------------------------------------------------------------*/

#define noDEBUG  // DEBUG or noDEBUG
#define ABS(X)                ( (X) >= 0 ? (X) : -(X) ) 
#define LIMIT(A,B,C)          ( (B) < (A) ? (A) : ( (B) > (C) ? (C) : (B) ) ) 


#include "flowRateWaterAlphaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::flowRateWaterAlphaFvPatchScalarField::flowRateWaterAlphaFvPatchScalarField
(
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF
 )
    :
        mixedFvPatchScalarField(p, iF),
        mode_("sharpen"),
        sharpen_(0.001),
        pressureThreshold_(VGREAT),   
        alphaLowerThreshold_(0.001),   
        alphaUpperThreshold_(0.999)   
{
    refValue() = patchInternalField();
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


    Foam::flowRateWaterAlphaFvPatchScalarField::flowRateWaterAlphaFvPatchScalarField
(
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF,
 const dictionary& dict
 )
    :
        mixedFvPatchScalarField(p, iF),
        mode_(dict.lookupOrDefault<word>("mode", "sharpen")),
        sharpen_(dict.lookupOrDefault<scalar>("sharpen",0.001))    ,
        pressureThreshold_(dict.lookupOrDefault<scalar>("pressureThreshold",VGREAT)),
        alphaLowerThreshold_(dict.lookupOrDefault<scalar>("alphaLowerThreshold",0.001)),    
        alphaUpperThreshold_(dict.lookupOrDefault<scalar>("alphaUpperThreshold",0.999))    
{
    fvPatchField<scalar>::operator=(patchInternalField());
    //    refValue() = *this;

    // Backward compatible to old "refValue"
    sharpen_=dict.lookupOrDefault<scalar>("refValue",sharpen_);
    sharpen_=dict.lookupOrDefault<scalar>("sharpen",sharpen_);

    refValue() = patchInternalField();
    refGrad()=0.0;
    valueFraction()=0.0;
}


    Foam::flowRateWaterAlphaFvPatchScalarField::flowRateWaterAlphaFvPatchScalarField
(
 const flowRateWaterAlphaFvPatchScalarField& ptf,
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF,
 const fvPatchFieldMapper& mapper
 )
    :
        mixedFvPatchScalarField(ptf, p, iF, mapper),
        mode_(ptf.mode_),
        sharpen_(ptf.sharpen_),
        pressureThreshold_(ptf.pressureThreshold_),
        alphaLowerThreshold_(ptf.alphaLowerThreshold_),
        alphaUpperThreshold_(ptf.alphaUpperThreshold_)
{}


    Foam::flowRateWaterAlphaFvPatchScalarField::flowRateWaterAlphaFvPatchScalarField
(
 const flowRateWaterAlphaFvPatchScalarField& tppsf
 )
    :
        mixedFvPatchScalarField(tppsf),
        mode_(tppsf.mode_),
        sharpen_(tppsf.sharpen_),
        pressureThreshold_(tppsf.pressureThreshold_),
        alphaLowerThreshold_(tppsf.alphaLowerThreshold_),
        alphaUpperThreshold_(tppsf.alphaUpperThreshold_)
{}


    Foam::flowRateWaterAlphaFvPatchScalarField::flowRateWaterAlphaFvPatchScalarField
(
 const flowRateWaterAlphaFvPatchScalarField& tppsf,
 const DimensionedField<scalar, volMesh>& iF
 )
    :
        mixedFvPatchScalarField(tppsf, iF),
        mode_(tppsf.mode_),
        sharpen_(tppsf.sharpen_),
        pressureThreshold_(tppsf.pressureThreshold_),
        alphaLowerThreshold_(tppsf.alphaLowerThreshold_),
        alphaUpperThreshold_(tppsf.alphaUpperThreshold_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::flowRateWaterAlphaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const fvPatchField<scalar>&  pp   = patch().lookupPatchField<volScalarField, scalar>("p");
    const fvsPatchField<scalar>& phip = patch().lookupPatchField<surfaceScalarField, scalar>("phi");

    // (1=DirichletBC) or (0=NeumannBC)?
    // Outbound flow: NeumannBC
    this->valueFraction() = neg(phip);

    // Value for inbound flow
    // "sharpen_" controls the sharpening of alpha.water by switching to 0 or 1 instead of Neumann=0
    this->refValue() = 
        sharpen_ * neg(phip)  * pos(this->patchInternalField()-0.5) 
        + (1 - sharpen_ * neg(phip)) *     this->patchInternalField();


    // If the pressure is above the threshold set alpha.water=1
    this->refValue() = pos(pp-pressureThreshold_) + neg(pp-pressureThreshold_) * this->refValue();

    // If alpha < alphaLowerThreshold, set alpha to 0., 
    this->refValue() = pos(this->refValue()-alphaLowerThreshold_) * this->refValue();

    // If alpha > alphaUpperThreshold, set alpha to 1., 
    this->refValue() = pos(this->refValue()-alphaUpperThreshold_) + neg(this->refValue()-alphaUpperThreshold_) * this->refValue();

    mixedFvPatchScalarField::updateCoeffs();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

    void Foam::flowRateWaterAlphaFvPatchScalarField::operator=
(
 const fvPatchField<scalar>& ptf
 )
{
    fvPatchField<scalar>::operator=
        (
         this->refValue() 
        );
}


void Foam::flowRateWaterAlphaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    // os.writeKeyword("mode") << mode_ << token::END_STATEMENT << nl;
    os.writeKeyword("sharpen") << sharpen_ << token::END_STATEMENT << nl;
    os.writeKeyword("pressureThreshold") << pressureThreshold_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaLowerThreshold") << alphaLowerThreshold_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaUpperThreshold") << alphaUpperThreshold_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
        (
         fvPatchScalarField,
         flowRateWaterAlphaFvPatchScalarField
        );
}

// ************************************************************************* //
