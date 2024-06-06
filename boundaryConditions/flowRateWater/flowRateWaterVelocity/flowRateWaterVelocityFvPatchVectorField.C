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
#define FLOWRATE_CHECKER
#define noFLOWRATE_CORRECTOR


#include "flowRateWaterVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::flowRateWaterVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flowRateWaterVelocityFvPatchVectorField::
    flowRateWaterVelocityFvPatchVectorField
(
 const fvPatch& p,
 const DimensionedField<vector, volMesh>& iF
 )
    :
        fixedValueFvPatchVectorField(p, iF),
        flowRateWater_(),
        inletDir_(0, 0, 0),    
        flowRateCorrector_(0.),
        relaxationTime_(0.0),
        avgU_old_(0.),
        lastTime_(0.),
        alphaLowerThreshold_(-GREAT),   
        alphaUpperThreshold_(0.999)   


{
}


Foam::flowRateWaterVelocityFvPatchVectorField::
    flowRateWaterVelocityFvPatchVectorField
(
 const fvPatch& p,
 const DimensionedField<vector, volMesh>& iF,
 const dictionary& dict
 )
    :
        fixedValueFvPatchVectorField(p, iF),
        flowRateWater_(Function1<scalar>::New("flowRateWater", dict)),
        inletDir_(dict.lookupOrDefault("inletDir", vector::zero)),
        flowRateCorrector_(dict.lookupOrDefault<scalar>("flowRateCorrector",1.0)),
        relaxationTime_(dict.lookupOrDefault<scalar>("relaxationTime",0.0)),
        avgU_old_(dict.lookupOrDefault<scalar>("avgU",0.0)),
        lastTime_(dict.lookupOrDefault<scalar>("lastTime",0.0)),
        alphaLowerThreshold_(dict.lookupOrDefault<scalar>("alphaLowerThreshold",-GREAT)),    
        alphaUpperThreshold_(dict.lookupOrDefault<scalar>("alphaUpperThreshold",0.999))    


{
    // fixedValueFvPatchVectorField::evaluate();
    //Initialise with the value entry if evaluation is not possible
    fvPatchVectorField::operator=
        (
         vectorField("value", dict, p.size())
        );

}


Foam::flowRateWaterVelocityFvPatchVectorField::
    flowRateWaterVelocityFvPatchVectorField
(
 const flowRateWaterVelocityFvPatchVectorField& ptf,
 const fvPatch& p,
 const DimensionedField<vector, volMesh>& iF,
 const fvPatchFieldMapper& mapper
 )
    :
        fixedValueFvPatchVectorField(ptf, p, iF, mapper),
        flowRateWater_(ptf.flowRateWater_->clone().ptr()),
        inletDir_(ptf.inletDir_),
        flowRateCorrector_(ptf.flowRateCorrector_),
        relaxationTime_(ptf.relaxationTime_),
        avgU_old_(ptf.avgU_old_),
        lastTime_(ptf.lastTime_),
        alphaLowerThreshold_(ptf.alphaLowerThreshold_),
        alphaUpperThreshold_(ptf.alphaUpperThreshold_)

{
    //     timeVsData_(ptf.timeVsData_.clone()),
}


Foam::flowRateWaterVelocityFvPatchVectorField::
    flowRateWaterVelocityFvPatchVectorField
(
 const flowRateWaterVelocityFvPatchVectorField& ptf
 )
    :
        fixedValueFvPatchVectorField(ptf),
        flowRateWater_(ptf.flowRateWater_->clone().ptr()),
        inletDir_(ptf.inletDir_),
        flowRateCorrector_(ptf.flowRateCorrector_),
        relaxationTime_(ptf.relaxationTime_),
        avgU_old_(ptf.avgU_old_),
        lastTime_(ptf.lastTime_),
        alphaLowerThreshold_(ptf.alphaLowerThreshold_),
        alphaUpperThreshold_(ptf.alphaUpperThreshold_)

{
}


Foam::flowRateWaterVelocityFvPatchVectorField::
    flowRateWaterVelocityFvPatchVectorField
(
 const flowRateWaterVelocityFvPatchVectorField& ptf,
 const DimensionedField<vector, volMesh>& iF
 )
    :
        fixedValueFvPatchVectorField(ptf, iF),
        flowRateWater_(ptf.flowRateWater_->clone().ptr()),
        inletDir_(ptf.inletDir_),
        flowRateCorrector_(ptf.flowRateCorrector_),
        relaxationTime_(ptf.relaxationTime_),
        avgU_old_(ptf.avgU_old_),
        lastTime_(ptf.lastTime_),
        alphaLowerThreshold_(ptf.alphaLowerThreshold_),
        alphaUpperThreshold_(ptf.alphaUpperThreshold_)

{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    void Foam::flowRateWaterVelocityFvPatchVectorField::autoMap
(
 const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
}


    void Foam::flowRateWaterVelocityFvPatchVectorField::rmap
(
 const fvPatchVectorField& ptf,
 const labelList& addr
 )
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);
}


void Foam::flowRateWaterVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;

    const fvPatchField<scalar>& alphap = patch().lookupPatchField<volScalarField, scalar>("alpha.water");
    scalar flowrate = flowRateWater_->value(this->db().time().value());

    scalar A1;

    // Inlet direction is normal vector of each patch face, normalized, inbound positive
    vectorField n{-patch().nf()};

    // If inlet direction is prescribed use it.
    if (mag(inletDir_)>SMALL) { n = inletDir_ / mag(inletDir_); }

    // If the inlet direcion is not pointing in the right direction, set it to zero.
    n = neg(n & patch().nf()) * n;


    // Sum up the wetted boundary area, projected on inlet direction
    // Weighted with alpha filling and thresholds
    A1 = gSum(mag(patch().Sf() & n) * alphap * pos(neg(alphap-alphaLowerThreshold_)+pos(alphap-alphaUpperThreshold_)-SMALL) );

    if(A1<SMALL) {
        FatalErrorIn
            (
             "flowRateWaterVelocityFvPatchVectorField::updateCoeffs()"
            )   << " Wetted area on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()  << nl 
            << " is zero or inletDir is wrong!" << nl 
            << exit(FatalError);
    }


    // Evaluate mean velocity
    scalar avgU1 = flowrate/max(A1,SMALL);

    // relaxationTime_ = 10.; // Estimate 1/6*L/sqrt(H)

    // Evaluate relaxation factor 
    scalar relax = 1.0;
    if (relaxationTime_>SMALL) {
        relax = min(max(this->db().time().deltaT().value()/relaxationTime_, SMALL), 1.0);

        if (avgU_old_ < SMALL) avgU_old_  = avgU1;

        // Relax U
        avgU1  = relax * avgU1 + (1.0-relax) * avgU_old_;

    }

    // Do this once per timestep.  
    if (this->db().time().value() - lastTime_ > SMALL) {
        lastTime_ = this->db().time().value();
        avgU_old_ = avgU1;
    }


    // Check the flow rate
    fvPatchField<vector> U=this->patch().lookupPatchField<volVectorField, vector>("U"); 
    scalar realFlowRate = -gSum( U & patch().Sf() * alphap) + SMALL;

#ifdef FLOWRATE_CHECKER
    if(mag(realFlowRate-flowrate)/mag(realFlowRate+flowrate+SMALL)>0.001) {
        Info << "Warning: Flowrate on patch " << this->patch().name() << " is wrong by " << mag(200.*(realFlowRate-flowrate)/mag(realFlowRate+flowrate+SMALL)) << "%." << nl;
        if(relaxationTime_>SMALL) {
            Info << "Hint: flowRateRelaxationTime is set to " << relaxationTime_ << " s" << nl;
        }
    } 
#endif

#ifdef FLOWRATE_CORRECTOR
    if((mag(flowRateCorrector)<0.5)||(mag(flowRateCorrector)>2.)) {
        Info << "########################################################################################" <<  nl;
        Info << "Problem at the bounadry with flowRateWaterVelocity. Should vanish after a few timesteps." <<  nl;
        Info << "flowRateCorrector = "<<flowRateCorrector <<  nl;
        Info << "########################################################################################" <<  nl;
        flowRateCorrector=1.;
    } else {
        if (flowrate/realFlowRate > 1.1)  
            flowRateCorrector=flowRateCorrector*1.1;
        else if (flowrate/realFlowRate < 0.9)  
            flowRateCorrector=flowRateCorrector*0.9;
        else
            flowRateCorrector=flowRateCorrector*flowrate/realFlowRate;
    }
    // Correct U iteratively
    avgU1 = avgU1*flowRateCorrector;
#endif

    vectorField U_Boundary {n * avgU1};

    if (flowrate < 0.0) {
        // Is flow outbound? If yes, get rid of transversal momentum!
        // Flow direction in the domain
        vectorField U_Innen {this->patchInternalField()};
        // For outbound flow, add transversal velocity components to the BC
        U_Boundary  =  U_Boundary + mag(n) * ( U_Innen - (U_Innen & patch().nf()) * patch().nf());        
    } 

    // Take into account thresholds
    U_Boundary = U_Boundary * (pos(neg(alphap-alphaLowerThreshold_)+pos(alphap-alphaUpperThreshold_)-SMALL)+SMALL);

    fixedValueFvPatchVectorField::operator==
        (
         U_Boundary
        );

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::flowRateWaterVelocityFvPatchVectorField::write
(
 Ostream& os
 ) const
{
    fvPatchVectorField::write(os);
    flowRateWater_->writeData(os);
    os.writeKeyword("inletDir") << inletDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationTime") << relaxationTime_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaLowerThreshold") << alphaLowerThreshold_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphaUpperThreshold") << alphaUpperThreshold_ << token::END_STATEMENT << nl;
    os.writeKeyword("avgU") << avgU_old_ << token::END_STATEMENT << nl;

    writeEntry("value", os);

}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
        (
         fvPatchVectorField,
         flowRateWaterVelocityFvPatchVectorField
        );
}

// ************************************************************************* //
