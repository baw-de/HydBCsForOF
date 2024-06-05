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
 * Description
 *    Describes a volumetric flowrate of water in an interFoam simulation.
 * 
 *    Example of the boundary condition specification:
 *    @verbatim
 *    inlet
 *    {
 *        type            flowRateWaterVelocity;
 *        flowRate        0.2;        // Volumetric flow rate [m3/s]. Can also be a table ( (t1  q1)  (t2  q2)  (t3 q3) )
 * 
 *        // OPTIONAL ENTRIES:
 *        // inletDir     (1. 2. 3.);    // Vector for prescribed direction ([m] [m] [m])
 *        // relaxationTime    0.0;      // Relax adapation of the prescribed velocity [s]. Helpful against upstream waves. 
 *                                       // Can be estimated from: upstreamLength/sqrt(g * waterDepth)        
 *        // alphaLowerThreshold -9999.; // For all cells with alphaLowerThreshold < alpha < alphaUpperThreshold, the 
 *        // alphaUpperThreshold  0.999; //   velocity will be set to zero. Default: Stop motion for alpha < 0.999
 *        
 *        value           uniform (0 0 0); // placeholder
 *    }
 *    @endverbatim
 * 
 * Note
 *    - The value is positive inwards
 * 
 * \*---------------------------------------------------------------------------*/

#ifndef flowRateWaterVelocityFvPatchVectorField_H
#define flowRateWaterVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
     *     Class flowRateWaterVelocityFvPatchVectorField Declaration
     * \*---------------------------------------------------------------------------*/

    class flowRateWaterVelocityFvPatchVectorField
        :
            public fixedValueFvPatchVectorField
    {
        // Private Data


        //  table or single value for water flowrate
        autoPtr<Function1<scalar> > flowRateWater_;

        //- Inlet direction for water
        vector inletDir_;

        //- Inlet flow rate correction value
        scalar flowRateCorrector_;

        //- Time span for smearing the computed velocity
        scalar relaxationTime_;

        //- Old velocity
        scalar avgU_old_;

        //- Last execution time
        scalar lastTime_;

        // Alpha value above which U = 0
        scalar alphaLowerThreshold_;        

        // Alpha value below which U = 0
        scalar alphaUpperThreshold_;        



        // Private Member Functions

        //- Return current time
        scalar t() const;


        public:

        //- Runtime type information
        TypeName("flowRateWaterVelocity");


        // Constructors

        //- Construct from patch and internal field
        flowRateWaterVelocityFvPatchVectorField
            (
             const fvPatch&,
             const DimensionedField<vector, volMesh>&
            );

        //- Construct from patch, internal field and dictionary
        flowRateWaterVelocityFvPatchVectorField
            (
             const fvPatch&,
             const DimensionedField<vector, volMesh>&,
             const dictionary&
            );

        //- Construct by mapping onto a new patch
        flowRateWaterVelocityFvPatchVectorField
            (
             const flowRateWaterVelocityFvPatchVectorField&,
             const fvPatch&,
             const DimensionedField<vector, volMesh>&,
             const fvPatchFieldMapper&
            );

        //- Copy construct
        flowRateWaterVelocityFvPatchVectorField
            (
             const flowRateWaterVelocityFvPatchVectorField&
            );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
                (
                 new flowRateWaterVelocityFvPatchVectorField(*this)
                );
        }

        //- Construct as copy setting internal field reference
        flowRateWaterVelocityFvPatchVectorField
            (
             const flowRateWaterVelocityFvPatchVectorField&,
             const DimensionedField<vector, volMesh>&
            );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
            (
             const DimensionedField<vector, volMesh>& iF
            ) const
            {
                return tmp<fvPatchVectorField>
                    (
                     new flowRateWaterVelocityFvPatchVectorField
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
             const fvPatchVectorField&,
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


