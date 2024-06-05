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
 *    Describes a boundary conditionfor alpha.water in an interFoam simulation. This should be used for flow situations where the
 *    the flowrate or velocity is prescribed and the waterlevel is free.
 *    
 *    The basic behaviour is like zeroGradient. But due to numerical issues, zeroGradient produces sometimes instabilities
 *    for inbound flows. Thus, here addons are provided to avoid this.
 * 
 *    Example of the boundary condition specification:
 *    @verbatim
 *    inlet
 *    {
 *      type            flowRateWaterAlpha;    
 * 
 *      // OPTIONAL ENTRIES:
 *      // pressureThreshold    1.e99;          // For all boundary cells with pressure p above pressureThreshold, alpha.water=1
 *      //                                      // is prescribed. This avoids bubbles for cells with p > pressureThreshold.
 *      // alphaLowerThreshold  0.001;          // For all boundary cells with alpha below alphaUpperThreshold, alpha.water=0
 *      //                                      // is prescribed. This is a numerical trick to avoid rain drops. Low values increase
 *      //                                      // bubble risk, high values increase risk for jumps in waterlevel. Default 0.999
 *      // alphaUpperThreshold  0.999;          // For all boundary cells with alpha above alphaUpperThreshold, alpha.water=1
 *      //                                      // is prescribed. This is a numerical trick to avoid bubbles or water "creeping"
                                                // up a slope. High values increase bubble risk, low values increase risk for 
                                                // jumps in waterlevel. Default 0.999
 *      // sharpen              0.001;          // Weight between a threshold for alpha and a zeroGradient. This is a numerical
 *                                              // trick to avoid bubbles. Low values increase bubble risk, high values increase
 *                                              // risk for jumps in waterlevel. Default 0.001
 * 
 *      value           uniform 0.;
 *    }
 *    @endverbatim
 * 
 * 
 * \*---------------------------------------------------------------------------*/

#ifndef flowRateWaterAlphaFvPatchScalarField_H
#define flowRateWaterAlphaFvPatchScalarField_H

#include "mixedFvPatchFields.H"


namespace Foam
{

    /*---------------------------------------------------------------------------*\
     *                           Class flowRateWaterAlphaFvPatch Declaration
     * \*---------------------------------------------------------------------------*/

    class flowRateWaterAlphaFvPatchScalarField
        :
            public mixedFvPatchScalarField
    {
        // - Mode of operation (UNUSED!)
        word mode_;

        // Private data

        // Gradually replace zeroGradient values with 0 or 1
        scalar sharpen_;        

        // Pressure value to determine alpha.water=1
        scalar pressureThreshold_;        

        // Alpha value to determine alpha.water=0
        scalar alphaLowerThreshold_;        

        // Alpha value to determine alpha.water=1
        scalar alphaUpperThreshold_;        


        public:

        //- Runtime type information
        TypeName("flowRateWaterAlpha");


        // Constructors

        //- Construct from patch and internal field
        flowRateWaterAlphaFvPatchScalarField
            (
             const fvPatch&,
             const DimensionedField<scalar, volMesh>&
            );

        //- Construct from patch, internal field and dictionary
        flowRateWaterAlphaFvPatchScalarField
            (
             const fvPatch&,
             const DimensionedField<scalar, volMesh>&,
             const dictionary&
            );

        //- Construct by mapping given flowRateWaterAlphaFvPatchScalarField
        //  onto a new patch
        flowRateWaterAlphaFvPatchScalarField
            (
             const flowRateWaterAlphaFvPatchScalarField&,
             const fvPatch&,
             const DimensionedField<scalar, volMesh>&,
             const fvPatchFieldMapper&
            );

        //- Construct as copy
        flowRateWaterAlphaFvPatchScalarField
            (
             const flowRateWaterAlphaFvPatchScalarField&
            );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
                (
                 new flowRateWaterAlphaFvPatchScalarField(*this)
                );
        }

        //- Construct as copy setting internal field reference
        flowRateWaterAlphaFvPatchScalarField
            (
             const flowRateWaterAlphaFvPatchScalarField&,
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
                     new flowRateWaterAlphaFvPatchScalarField(*this, iF)
                    );
            }


        // Member functions

        // Evaluation functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;

        // Member operators

        virtual void operator=(const fvPatchField<scalar>& pvf);
    };


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //


