/*---------------------------------------------------------------------------*\
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

#include "waterLevel_alpha_prghFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"
#include "fvMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waterLevel_alpha_prghFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waterLevel_alpha_prghFvPatchScalarField::
    waterLevel_alpha_prghFvPatchScalarField
(
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF
 )
    :
        fixedValueFvPatchScalarField(p, iF),
        mode_("dummy"),
        my_data_(),
        flowRateMultiplier_(1.0),
        refValue_(0.0),
        refPoint_(vector::zero),
        relaxationTime_(-100.0),
        dynamicPressureCorrection_(1.0),
        lastTime_(0.0),
        useOldTime_(0.0),
        mean_water_dynamic_pressure_old_(0.0)

{
}

Foam::waterLevel_alpha_prghFvPatchScalarField::
    waterLevel_alpha_prghFvPatchScalarField
(
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF,
 const dictionary& dict
 )
    :
        fixedValueFvPatchScalarField(p, iF),
        mode_(dict.lookupOrDefault<word>("mode", "dummy")),
        my_data_(Function1<scalar>::New("data", dict)),
        flowRateMultiplier_(dict.lookupOrDefault<scalar>("flowRateMultiplier", 1.0)),
        refValue_(dict.lookupOrDefault<scalar>("refValue", 0.0)),
        refPoint_(dict.lookupOrDefault<vector>("refPoint", vector::zero)),
        relaxationTime_(dict.lookupOrDefault<scalar>("dynamicPressureCorrectionRelaxationTime", -100.0)),
        dynamicPressureCorrection_(dict.lookupOrDefault<scalar>("dynamicPressureCorrection", 1.0)),
        lastTime_(0.0),
        useOldTime_(0.0),
        mean_water_dynamic_pressure_old_(dict.lookupOrDefault<scalar>("meanWaterDynamicPressureOld", 0.0))

{
    if (mode_=="ratingCurveFunction")  {
        if(mag(my_data_->value(0.)+my_data_->value(1.)+my_data_->value(2.)+my_data_->value(3.))<SMALL ) {
            FatalErrorIn ("waterLevel_alpha_prgh")
                << " Parameters in table must be set for ratingCurveFunction."  << ".\n"
                << " Error on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }
    } else 
        if (mode_=="ratingCurveTable")  {          
        } else
            if (mode_=="timeSeries")  {          
            } else
                if (mode_=="singleValue")  {          
                } else {
                    FatalErrorIn ("waterLevel_alpha_prgh")
                        << " mode must be singleValue, ratingCurveTable, ratingCurveFunction or timeSeries, not "  << mode_   << ".\n"
                        << " Error on patch " << this->patch().name()
                        << " of field " << this->internalField().name()
                        << " in file " << this->internalField().objectPath()
                        << exit(FatalError);
                }

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
            (
             scalarField("value", dict, p.size())
            );
    }
    else
    {
        FatalErrorIn ("waterLevel_alpha_prgh")
            << " Error on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file: " << this->internalField().objectPath()
            << " value not found! "
            << exit(FatalError);
    }

    // For setFields and for initialization, apply a zero gradient condition:
    fixedValueFvPatchScalarField::operator == ( this->patchInternalField() );

}


Foam::waterLevel_alpha_prghFvPatchScalarField::
    waterLevel_alpha_prghFvPatchScalarField
(
 const waterLevel_alpha_prghFvPatchScalarField& ptf,
 const fvPatch& p,
 const DimensionedField<scalar, volMesh>& iF,
 const fvPatchFieldMapper& mapper
 )
    :
        fixedValueFvPatchScalarField(ptf, p, iF, mapper),
        mode_(ptf.mode_),
        my_data_(ptf.my_data_->clone().ptr()),
        flowRateMultiplier_(ptf.flowRateMultiplier_),
        refValue_(ptf.refValue_),
        refPoint_(ptf.refPoint_),
        relaxationTime_(ptf.relaxationTime_),
        dynamicPressureCorrection_(ptf.dynamicPressureCorrection_),
        lastTime_(ptf.lastTime_),
        useOldTime_(ptf.lastTime_),
        mean_water_dynamic_pressure_old_(ptf.mean_water_dynamic_pressure_old_)


{
}


Foam::waterLevel_alpha_prghFvPatchScalarField::
    waterLevel_alpha_prghFvPatchScalarField
(
 const waterLevel_alpha_prghFvPatchScalarField& ptf
 )
    :
        fixedValueFvPatchScalarField(ptf),
        mode_(ptf.mode_),
        my_data_(ptf.my_data_->clone().ptr()),
        flowRateMultiplier_(ptf.flowRateMultiplier_),
        refValue_(ptf.refValue_),
        refPoint_(ptf.refPoint_),
        relaxationTime_(ptf.relaxationTime_),
        dynamicPressureCorrection_(ptf.dynamicPressureCorrection_),
        lastTime_(ptf.lastTime_),
        useOldTime_(ptf.lastTime_),
        mean_water_dynamic_pressure_old_(ptf.mean_water_dynamic_pressure_old_)

{
}


Foam::waterLevel_alpha_prghFvPatchScalarField::
    waterLevel_alpha_prghFvPatchScalarField
(
 const waterLevel_alpha_prghFvPatchScalarField& ptf,
 const DimensionedField<scalar, volMesh>& iF
 )
    :
        fixedValueFvPatchScalarField(ptf, iF),
        mode_(ptf.mode_),
        my_data_(ptf.my_data_->clone().ptr()),
        flowRateMultiplier_(ptf.flowRateMultiplier_),
        refValue_(ptf.refValue_),
        refPoint_(ptf.refPoint_),
        relaxationTime_(ptf.relaxationTime_),
        dynamicPressureCorrection_(ptf.dynamicPressureCorrection_),
        lastTime_(ptf.lastTime_),
        useOldTime_(ptf.lastTime_),
        mean_water_dynamic_pressure_old_(ptf.mean_water_dynamic_pressure_old_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void Foam::waterLevel_alpha_prghFvPatchScalarField::autoMap
(
 const fvPatchFieldMapper& m
 )
{
}


    void Foam::waterLevel_alpha_prghFvPatchScalarField::rmap
(
 const fvPatchScalarField& ptf,
 const labelList& addr
 )
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::waterLevel_alpha_prghFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const uniformDimensionedVectorField& g = db().time().lookupObject<uniformDimensionedVectorField>("g");
    const fvPatchField<scalar>&  ghp    = this->patch().lookupPatchField<volScalarField, scalar>("gh");
    const fvPatchField<scalar>&  alphap = this->patch().lookupPatchField<volScalarField, scalar>("alpha.water");
    const fvPatchField<scalar>&  rhop   = this->patch().lookupPatchField<volScalarField, scalar>("rho");
    const fvsPatchField<scalar>& phip   = this->patch().lookupPatchField<surfaceScalarField, scalar>("phi");
    const volScalarField& rhoScalarField   = db().lookupObjectRef<volScalarField>("rho");
    const fvMesh & mesh =this->internalField().mesh();

    scalarField physical_pressure_;
    scalarField dynamic_pressure_;
    scalarField alphaTimeWeighed_;

    scalar waterlevel=0.;
    scalar point_offset;
    scalar pressure_offset;
    scalar flowRateAlphaPhi=0.;
    scalar rho1=0.0; // Density of lighter fluid

    scalar mean_water_dynamic_pressure=0.;
    scalar relax=1.;
    vector waterlevelPoint; 

    // How to get the real density of the light fluid? What, if we are not running the solver but e.g. mapFields?
    rho1=gSum(pos0(alphap)*neg(alphap-SMALL) * (rhop ) * patch().magSf()) / ( gSum( pos0(alphap)*neg(alphap-SMALL) *        patch().magSf()) + SMALL);    
    if (rho1 < SMALL) {
        // Info << "Problem in waterLevel_alpha_prgh: No air found on patch " << this->patch().name()
        //     << ". Using global minimum of rho instead." << endl;
        rho1=gMin(rhoScalarField);
    }

    // Compute flowrate of water at the boundary. "Outflow" has a positive sign.
    flowRateAlphaPhi = flowRateMultiplier_ * gSum( (alphap) * phip);

    // Switch based on operation mode
    if (mode_=="singleValue")  {          
        waterlevel = my_data_->value(this->db().time().value()); 
    } else 
        if (mode_=="ratingCurveFunction")  {
            // Compute desired water level from rating curve via formula
            waterlevel = my_data_->value(1.) + sign(flowRateAlphaPhi) * my_data_->value(2.) * pow((my_data_->value(0.) * mag(flowRateAlphaPhi)), my_data_->value(3.));
        } else
            if (mode_=="ratingCurveTable")  {          
                // Compute desired water level from rating curve via table
                waterlevel =  my_data_->value(flowRateAlphaPhi);
            } else
                if (mode_=="timeSeries")  {          
                    // Compute desired water level for current time from curve via table
                    waterlevel =  my_data_->value(this->db().time().value());
                } else {
                    //          Do something silly
                    FatalErrorIn ("waterLevel_alpha_prgh") << " Unexpected mode selection. Something went wrong!" << exit(FatalError);
                }

    // Transfer waterlevel to vector 
    waterlevelPoint = waterlevel * (-g.value()/mag(g.value()));  

    if (this->internalField().name()=="p_rgh") {
        // Geometrical offset in order to match the pressure at the reference point.
        // Necessary in order to avoid pressure differences between boundaries.
        if(mag(refPoint_)<SMALL) {
            // Fallback: Highest point of the mesh is used as reference point  
            refPoint_  =  (min((mesh.C() & g.value()/mag(g.value())))*g.value()/mag(g.value())).value();
        }

        point_offset = ((g.value()&waterlevelPoint)-(g.value()&refPoint_))/mag(g.value());

        // Pressure offset in order to match the pressure at the reference point
        if (point_offset > -1.e99 ) {  // > 0? 
            // pressure_offset = refValue_ + rho1*mag(g.value())*point_offset;  // In air
            pressure_offset = refValue_ + rho1*((g.value()&waterlevelPoint)-(g.value()&refPoint_));  // In air
        } else {
            FatalErrorIn ("waterLevel_alpha_prgh")
                << " refPoint must be above the prescribed water level."  << ".\n"
                << " Error on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalError);
        }

        // Calculate dynamic pressure stabilization for inbound flow.
        dynamic_pressure_  = neg(phip)*0.5*rhop*pow(phip/(patch().magSf()+SMALL),2.0);

        // Do this once per timestep.  
        if (this->db().time().value() - lastTime_ > SMALL) {
            lastTime_ = this->db().time().value();

            Info << "waterLevel BC " <<  this->internalField().name() << " " << this->patch().name() << " Values: Q_out=" << flowRateAlphaPhi << " m³/s; h_set=" << waterlevel  << " m"    << endl;

            // Calculate dynamic pressure stabilization correction for inbound flow.
            mean_water_dynamic_pressure = dynamicPressureCorrection_ * gSum(   (alphap)  * patch().magSf() * dynamic_pressure_) / ( gSum(    (alphap)  * patch().magSf()) + SMALL);


            if (relaxationTime_ > SMALL) {
                // Evaluate relaxation factor 
                relax = min(max(this->db().time().deltaT().value()/relaxationTime_, SMALL), 1.0);
            } 
            if (mag(relaxationTime_) < SMALL) {
                // No relaxation
                relax = 1.;
            } 
            if (relaxationTime_ < -1) {
                // Hard set relaxation. If inflow is there, relaxation MUST be used for stability reasons!
                relax = -1./relaxationTime_;
            } 

            // Initialise with useful value) 
            if (mean_water_dynamic_pressure_old_< SMALL) mean_water_dynamic_pressure_old_ = mean_water_dynamic_pressure;

            // Relax mean_water_dynamic_pressure
            mean_water_dynamic_pressure_old_ = relax * mean_water_dynamic_pressure + (1.0-relax) * mean_water_dynamic_pressure_old_;
            Info << "Dynamic pressure compensation: " << mean_water_dynamic_pressure_old_ << endl;

        }

        // Compute pressure field at the boundary
        //  - Desired hydrostatic profile
        //  - Apply offset necessary to match the reference pressure
        //  - totalPressure stabilization for inbound flow
        //  - correction for totalPressure stabilization for inbound flow
        //    rhop*((g.value()&patch().Cf())-(g.value()&waterlevelPoint))

        physical_pressure_ =
            ( 
             rhop*((g.value()&patch().Cf())-(g.value()&waterlevelPoint))
             + pressure_offset
             - dynamic_pressure_ 
             + alphap * mean_water_dynamic_pressure_old_
            );


        // Compute p_rgh from physical pressure and set on the boundary as a Dirichlet BC
        fixedValueFvPatchScalarField::operator==
            (
             physical_pressure_ 
             - rhop*ghp
            );
    }  
    else 
        if ((this->internalField().name()=="alpha.water")||(this->internalField().name()=="alpha.air")) {
            // Set alpha.water on the boundary to either 0 or 1 as a Dirichlet BC for inbound flow. Use Neumann-Zero for outbound flow.

            // Do this once per timestep.  
            if (this->db().time().value() - lastTime_ > SMALL) {
                lastTime_ = this->db().time().value();
                Info << "waterLevel BC " <<  this->internalField().name() << " " << this->patch().name() << " Values: Q=" << flowRateAlphaPhi << " m³/s; h_set=" << waterlevel  << " m"    << endl;
            }

            if (this->internalField().name()=="alpha.water") 
                fixedValueFvPatchScalarField::operator==
                    (
                     neg0(phip)*(pos((g.value()&patch().Cf())-(g.value()&waterlevelPoint)))
                     + pos(phip)*this->patchInternalField()
                    );

            if (this->internalField().name()=="alpha.air") 
                fixedValueFvPatchScalarField::operator==
                    (
                     neg0(phip)*(neg((g.value()&patch().Cf())-(g.value()&waterlevelPoint)))
                     + pos(phip)*this->patchInternalField()
                    );

        }
        else {
            FatalErrorIn ("waterLevel_alpha_p_rgh")
                << " You should not use this boundary condition with field" << this->internalField().name()
                << " on patch " << this->patch().name()
                << " in file "  << this->internalField().objectPath()
                << exit(FatalError);
        }      


    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::waterLevel_alpha_prghFvPatchScalarField::write
(
 Ostream& os
 ) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("mode") << mode_ << token::END_STATEMENT << nl;
    my_data_->writeData(os);
    os.writeKeyword("flowRateMultiplier") << flowRateMultiplier_ << token::END_STATEMENT << nl;
    os.writeKeyword("useOldTime") << useOldTime_ << token::END_STATEMENT << nl;

    if (this->internalField().name()=="p_rgh") {
        os.writeKeyword("dynamicPressureCorrection") << dynamicPressureCorrection_ << token::END_STATEMENT << nl;
        os.writeKeyword("dynamicPressureCorrectionRelaxationTime") << relaxationTime_  << token::END_STATEMENT << nl;
        os.writeKeyword("refPoint") << refPoint_ << token::END_STATEMENT << nl;
        os.writeKeyword("refValue") << refValue_ << token::END_STATEMENT << nl;
        os.writeKeyword("meanWaterDynamicPressureOld") << mean_water_dynamic_pressure_old_ << token::END_STATEMENT << nl;
    }

    writeEntry("value", os);


}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
        (
         fvPatchScalarField,
         waterLevel_alpha_prghFvPatchScalarField
        );
}

// ************************************************************************* //

