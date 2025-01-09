/*---------------------------------------------------------------------------*\
  Author
    Carsten Thorenz 2024, Federal Waterways Research Institute 
  
  Publication
    Thorenz, C. (2024): 'Boundary Conditions for Hydraulic Structures Modelling with OpenFOAM',
    10th International Symposium on Hydraulic Structures, Zürich. ISSN 0374-0056 , DOI: 10.3929/ethz-b-000675949
 
  License
     GPL v3
 
  Description
    This boundary condition provides a water level condition, which prescribes
    values both for alpha.water and p_rgh in order to maintain a certain water
    level. For U a boundary condition with "Neumann zero" character should be 
    used. "pressureInletOutletVelocity" is recommended for U.
   
    The prescribed water level can be:
    
    - a fixed single value
    - a table with a time series of values
    - a series of values which depend on the flow rate of water ("rating curve"), defined by a table 
    - a series of values which depend on the water flown out ("storage table"), defined by a table 
        
    The boundary condition is primarily intended for outbound flow. For inbound flow,
    special treatment is necessary to stabilize the solution. This introduces errors, 
    which can mostly be compensated with the "dynamicPressureCorrection" parameter. If 
    dynamicPressureCorrection is set to zero, the boundary condition will prescribe the 
    total pressure (energy height) for inbound flow and no longer the water level.
    
    For inbound flow, the parameter "dynamicPressureCorrectionRelaxationTime" must be set to 
    a useful value in order to avoid self-exciting water level oscillations. Estimate: 
    
       dynamicPressureCorrectionRelaxationTime=2*L/sqrt(g*h) 
       
    with L as length from the BC patch to the next major wave reflecting boundary and h as 
    the water depth.        
     
   
Usage
    \table
        Property      | Description                   | Req'd | Default
        mode          | Mode of operation             | yes | wordDefault
        data          | Value or table for water level| yes |
        flowRateMultiplier| Multiply Q before h evaluation| no  |
        refPoint      | Reference point for pressure  | no  |
        refValue      | Reference value for pressure  | no  |
        dynamicPressureCorrection| Correction for inbound flow stabilization| no  | 1.0
        dynamicPressureCorrectionRelaxationTime| Relaxation for inbound flow   | no  |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            waterLevel_alpha_prgh;
        mode            timeSeries;         //  singleValue, ratingCurveTable, timeSeries or waterStorageTable 
        // data            constant 26.25;  // water level for single value
         
         
        data table (     
                  (0.    123.)   
                  (100.  123.4)   
                  (2000. 123.)   
             );                             // Table with data for time series ([s] [m]), rating curve ([m^3/s] [m]), or storage curve  ([m^3] [m]),   
        flowRateMultiplier  1.;             // Optional: Multiplier for the computed flowrate before evaluating flow rate dependent level
        refPoint        (0 0 0);            // Optional: Point in the air, where pressure has to have a certain value
        refValue        0;                  // Optional: Pressure in the air at refPoint
        dynamicPressureCorrection 1;        // Optional: Amount of correction of inbound flow pressure reduction
        dynamicPressureCorrectionRelaxationTime  100; // Optional: Relax correction of inbound flow pressure reduction. Set to period of long waves. 
        value           nonuniform List<scalar>
    }
    \endverbatim

SourceFiles
    waterLevel_alpha_prghFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef waterLevel_alpha_prghFvPatchScalarField_H
#define waterLevel_alpha_prghFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class waterLevel_alpha_prghFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class waterLevel_alpha_prghFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data
        // - Mode of operation
        word mode_;

        // qh relation table or h(t)-table or single value or ...
        autoPtr<Function1<scalar> > my_data_;

        //- Multiplier for measured flowrate before using qh-relation
        scalar flowRateMultiplier_;

        //- Reference pressure
        scalar refValue_;

        //- Reference pressure location
        vector refPoint_;
        
        //- Time span for smearing dynamicPressureCorrection
        scalar dynamicPressureCorrectionRelaxationTime_;

        //- Factor of total pressure stabilization correction
        scalar dynamicPressureCorrection_;
        
        //- Last execution time
        scalar lastTime_;
        
        //- Sum of water that left over the boundary
        scalar waterStored_;
        
        // Stored dynamic pressure
        scalar mean_water_dynamic_pressure_old_;

    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("waterLevel_alpha_prgh");


    // Constructors

        //- Construct from patch and internal field
        waterLevel_alpha_prghFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        waterLevel_alpha_prghFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        waterLevel_alpha_prghFvPatchScalarField
        (
            const waterLevel_alpha_prghFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        waterLevel_alpha_prghFvPatchScalarField
        (
            const waterLevel_alpha_prghFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new waterLevel_alpha_prghFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        waterLevel_alpha_prghFvPatchScalarField
        (
            const waterLevel_alpha_prghFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new waterLevel_alpha_prghFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
