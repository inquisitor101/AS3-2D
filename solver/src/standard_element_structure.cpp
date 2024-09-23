#include "standard_element_structure.hpp"



//-----------------------------------------------------------------------------------
// CStandardElement member functions.
//-----------------------------------------------------------------------------------


CStandardElement::CStandardElement
(
 ETypeDOF       typeDOFs,
 unsigned short nPoly,
 unsigned short nInt 
)
	:
		mTypeDOFsSol(typeDOFs), mNPolySol(nPoly), mNInt1D(nInt)
 /*
	* Constructor for the standard element, which initiates an element in parametric
	* space for an interface boundary.
	*/
{
	// Deduce the number os solution and integration points.
	mNSol1D = mNPolySol+1;
	mNSol2D = mNSol1D*mNSol1D;
	mNInt2D = mNInt1D*mNInt1D;

	// Compute the solution DOFs in 1D and in parametric space.
	ComputeLocationDOFs1D(mTypeDOFsSol, mNSol1D, mRSol1D);

	// Compute the integration data.
	InitializeQuadrature(mNInt1D, mRInt1D, mWInt1D, mWInt2D);

	// Initialize the Vandermonde matrices in 1D. Note, the bases are the orthonormal Legendre polynomials.
	NPolynomialUtility::OrthonormalLegendreBasis1D( mNSol1D, mRSol1D, mVandermondeSol1D ); // At solution    points.
	NPolynomialUtility::OrthonormalLegendreBasis1D( mNSol1D, mRInt1D, mVandermondeInt1D ); // At integration points.

	// Initialize the derivative of the Vandermonde matrices in 1D.
	NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D( mNSol1D, mRSol1D, mDerVandermondeSol1D ); // Solution    points.
	NPolynomialUtility::DerivativeOrthonormalLegendreBasis1D( mNSol1D, mRInt1D, mDerVandermondeInt1D ); // Integration points.

	// Initialize the Lagrange polynomials in 1D.
	NPolynomialUtility::LagrangeBasis1D( mVandermondeInt1D, mVandermondeSol1D, mLagrangeInt1D );

	// Initialize the derivative of the Lagrange polynomials in 1D.
	NPolynomialUtility::DerivativeLagrangeBasis1D( mDerVandermondeSol1D, mVandermondeSol1D, mDerLagrangeSol1D );
	NPolynomialUtility::DerivativeLagrangeBasis1D( mDerVandermondeInt1D, mVandermondeSol1D, mDerLagrangeInt1D );

	// Initialize the derivative of the Lagrange polynomials at the solution min and max faces in 1D.  
	NPolynomialUtility::DerivativeLagrangeBasisFace1D( mDerLagrangeSol1D, mDerLagrangeMinFace1D, mDerLagrangeMaxFace1D );

	// Compute the transpose of the Lagrange polynomial at integration points in 1D.
	NLinearAlgebra::TransposeAS3Matrix(mLagrangeInt1D, mLagrangeInt1DTrans);
	// Compute the transpose of the derivative Lagrange polynomial at integration points in 1D.
	NLinearAlgebra::TransposeAS3Matrix(mDerLagrangeInt1D, mDerLagrangeInt1DTrans);
	// Compute the transpose of the derivative Lagrange polynomial at solution points in 1D.
	NLinearAlgebra::TransposeAS3Matrix(mDerLagrangeSol1D, mDerLagrangeSol1DTrans);
}

//-----------------------------------------------------------------------------------

CStandardElement::~CStandardElement
(
 void
)
 /*
	* Destructor, which cleans up after the standard element class.
	*/
{

}

//-----------------------------------------------------------------------------------

void CStandardElement::InitializeQuadrature
(
 unsigned int           nInt1D,
 CMatrixAS3<as3double> &rInt1D,
 CMatrixAS3<as3double> &wInt1D,
 CMatrixAS3<as3double> &wInt2D
)
 /*
	* Function that computes the integration properties on the reference element.
	*/
{
	// Allocate memory for the integration points and weights in 1D.
	rInt1D.resize(nInt1D);
	wInt1D.resize(nInt1D);

	// Allocate memory for temporary vectors, since they are double-precision only in 
	// the quadrature structure. NOTE, CGaussJacobiQuadrature must be refactored.
	as3vector1d<double> rtmp(nInt1D);
	as3vector1d<double> wtmp(nInt1D);

	// Compute the quadrature points and weights in 1D.
	CGaussJacobiQuadrature::GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rtmp, wtmp);

	// Copy the temporary data to their correct vector, in the proper precision.
	for(size_t i=0; i<nInt1D; i++)
	{
		rInt1D[i] = static_cast<as3double>( rtmp[i] );
		wInt1D[i] = static_cast<as3double>( wtmp[i] );
	}

	// Allocate the memory for the integration weights in 2D.
	wInt2D.resize(mNInt2D);

	// Form the integration weights in 2D, which are tensor-products of 1D vectors.
	size_t s = 0;
	for(size_t j=0; j<nInt1D; j++)
		for(size_t i=0; i<nInt1D; i++)
			wInt2D[s++] = wInt1D[i]*wInt1D[j];
}

//-----------------------------------------------------------------------------------

void CStandardElement::ComputeLocationDOFs1D
(
 ETypeDOF               typeDOFs,
 unsigned int           nDOFs1D,
 CMatrixAS3<as3double> &rDOFs1D
)
 /*
	* Function that computes the locations of the DOFs in 1D, over the interval [-1, +1].
	*/
{
	// Allocate memory for the nodes.
	rDOFs1D.resize(nDOFs1D);

	// Relative tolerance.
	const as3double TOL = 1.0e-10;

	// Determine what type of distribution the DOFs obey.
	switch( typeDOFs )
	{
		// Equidistant distribution.
		case(ETypeDOF::EQD):
		{
			// Step size.
			const as3double dh = 2.0/( (as3double) nDOFs1D - 1.0 );

			// Starting position.
			const as3double r0 = -1.0;

			for(unsigned int i=0; i<nDOFs1D; i++)
				rDOFs1D[i] = (as3double) i*dh + r0;

			// Hard-enforce last node.
			rDOFs1D[nDOFs1D-1] = 1.0;

			break;
		}

		// Legendre-Gauss-Lobatto distribution.
		case(ETypeDOF::LGL):
		{
			switch( nDOFs1D ){
				case 1:
				{
					rDOFs1D[0] =  0.0;
					break;
				}

				case 2:
				{
					rDOFs1D[0] = -1.0; 
					rDOFs1D[1] =  1.0;
					break;
				}

				case 3:
				{
					rDOFs1D[0] = -1.0; 
					rDOFs1D[1] =  0.0; 
					rDOFs1D[2] =  1.0;
					break;
				}

				case 4:
				{
					const as3double t0 = std::sqrt(1.0/5.0);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t0; 
					rDOFs1D[2] =  t0;
					rDOFs1D[3] =  1.0;
					break;
				}

				case 5:
				{
					const as3double t0 = std::sqrt(3.0/7.0);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t0; 
					rDOFs1D[2] =  0.0; 
					rDOFs1D[3] =  t0;
					rDOFs1D[4] =  1.0;
					break;
				}

				case 6:
				{
					const as3double t1 = 2.0*std::sqrt(7.0)/21.0;
					const as3double t2 = std::sqrt(1.0/3.0 + t1);
					const as3double t3 = std::sqrt(1.0/3.0 - t1);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t2; 
					rDOFs1D[2] = -t3;
					rDOFs1D[3] =  t3; 
					rDOFs1D[4] =  t2;
					rDOFs1D[5] =  1.0;
					break;
				}

				case 7:
				{
					const as3double t1 = 2.0*std::sqrt(5.0/3.0)/11.0;
					const as3double t2 = std::sqrt(5.0/11.0 + t1);
					const as3double t3 = std::sqrt(5.0/11.0 - t1);
					rDOFs1D[0] = -1.0;
					rDOFs1D[1] = -t2; 
					rDOFs1D[2] = -t3;
					rDOFs1D[3] =  0.0;
					rDOFs1D[4] =  t3; 
					rDOFs1D[5] =  t2;
					rDOFs1D[6] =  1.0;

					break;
				}

				case 8:
				{
					rDOFs1D[0] = -1.0;
          rDOFs1D[1] = -0.8717401485096066153375;
          rDOFs1D[2] = -0.5917001814331423021445;
          rDOFs1D[3] = -0.2092992179024788687687;
          rDOFs1D[4] =  0.2092992179024788687687;
          rDOFs1D[5] =  0.5917001814331423021445;
          rDOFs1D[6] =  0.8717401485096066153375;
          rDOFs1D[7] =  1.0;

					break;
				}

				case 9:
				{
          rDOFs1D[0] = -1.0;
          rDOFs1D[1] = -0.8997579954114601573124;
          rDOFs1D[2] = -0.6771862795107377534459;
          rDOFs1D[3] = -0.3631174638261781587108;
          rDOFs1D[4] =  0.0;
          rDOFs1D[5] =  0.3631174638261781587108;
          rDOFs1D[6] =  0.6771862795107377534459;
          rDOFs1D[7] =  0.8997579954114601573124;
          rDOFs1D[8] =  1.0;

					break;
				}

				case 10:
				{
          rDOFs1D[0] = -1.0;
          rDOFs1D[1] = -0.9195339081664588138289;
          rDOFs1D[2] = -0.7387738651055050750031;
          rDOFs1D[3] = -0.4779249498104444956612;
          rDOFs1D[4] = -0.1652789576663870246262;
          rDOFs1D[5] =  0.1652789576663870246262;
          rDOFs1D[6] =  0.4779249498104444956612;
          rDOFs1D[7] =  0.7387738651055050750031;
          rDOFs1D[8] =  0.9195339081664588138289;
          rDOFs1D[9] =  1.0;

					break;
				}

				case 11:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.93400143040805913433227413609938;
					rDOFs1D[2]  = -0.78448347366314441862241781610846;
					rDOFs1D[3]  = -0.56523532699620500647096396947775;
 					rDOFs1D[4]  = -0.29575813558693939143191151555906;
          rDOFs1D[5]  =  0.0;
					rDOFs1D[6]  =  0.29575813558693939143191151555906;
					rDOFs1D[7]  =  0.56523532699620500647096396947775;
					rDOFs1D[8]  =  0.78448347366314441862241781610846;
					rDOFs1D[9]  =  0.93400143040805913433227413609938;
					rDOFs1D[10] =  1.0;

					break;
				}

				case 12:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.94489927222288222340758013830322;
					rDOFs1D[2]  = -0.81927932164400667834864158171690;
					rDOFs1D[3]  = -0.63287615303186067766240485444366;
					rDOFs1D[4]  = -0.39953094096534893226434979156697;
					rDOFs1D[5]  = -0.13655293285492755486406185573969;
					rDOFs1D[6]  =  0.13655293285492755486406185573969;
					rDOFs1D[7]  =  0.39953094096534893226434979156697;
					rDOFs1D[8]  =  0.63287615303186067766240485444366;
					rDOFs1D[9]  =  0.81927932164400667834864158171690;
					rDOFs1D[10] =  0.94489927222288222340758013830322;
					rDOFs1D[11] =  1.0;

					break;
				}

				case 13:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.95330984664216391189690546475545;
					rDOFs1D[2]  = -0.84634756465187231686592560709875;
					rDOFs1D[3]  = -0.68618846908175742607275903956636;
					rDOFs1D[4]  = -0.48290982109133620174693723363693;
					rDOFs1D[5]  = -0.24928693010623999256867370037423;
					rDOFs1D[6]  =  0.0;
					rDOFs1D[7]  =  0.24928693010623999256867370037423;
					rDOFs1D[8]  =  0.48290982109133620174693723363693;
					rDOFs1D[9]  =  0.68618846908175742607275903956636;
					rDOFs1D[10] =  0.84634756465187231686592560709875;
					rDOFs1D[11] =  0.95330984664216391189690546475545;
					rDOFs1D[12] =  1.0;

					break;
				}

				case 14:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.95993504526726090135510016201542;
					rDOFs1D[2]  = -0.86780105383034725100022020290826;
					rDOFs1D[3]  = -0.72886859909132614058467240052088;
					rDOFs1D[4]  = -0.55063940292864705531662270585908;
					rDOFs1D[5]  = -0.34272401334271284504390340364167;
					rDOFs1D[6]  = -0.11633186888370386765877670973616;
					rDOFs1D[7]  =  0.11633186888370386765877670973616;
					rDOFs1D[8]  =  0.34272401334271284504390340364167;
					rDOFs1D[9]  =  0.55063940292864705531662270585908;
					rDOFs1D[10] =  0.72886859909132614058467240052088;
					rDOFs1D[11] =  0.86780105383034725100022020290826;
					rDOFs1D[12] =  0.95993504526726090135510016201542;
					rDOFs1D[13] =  1.0;

					break;
				}

				case 15:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.96524592650383857279585139206960;
					rDOFs1D[2]  = -0.88508204422297629882540163148223;
					rDOFs1D[3]  = -0.76351968995181520070411847597629;
					rDOFs1D[4]  = -0.60625320546984571112352993863673;
					rDOFs1D[5]  = -0.42063805471367248092189693873858;
					rDOFs1D[6]  = -0.21535395536379423822567944627292;
					rDOFs1D[7]  =  0.0;
					rDOFs1D[8]  =  0.21535395536379423822567944627292;
					rDOFs1D[9]  =  0.42063805471367248092189693873858;
					rDOFs1D[10] =  0.60625320546984571112352993863673;
					rDOFs1D[11] =  0.76351968995181520070411847597629;
					rDOFs1D[12] =  0.88508204422297629882540163148223;
					rDOFs1D[13] =  0.96524592650383857279585139206960;
					rDOFs1D[14] =  1.0;

					break;
				}

				case 16:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.96956804627021793295224273836746;
					rDOFs1D[2]  = -0.89920053309347209299462826151985;
					rDOFs1D[3]  = -0.79200829186181506393108827096315;
					rDOFs1D[4]  = -0.65238870288249308946788321964058;
					rDOFs1D[5]  = -0.48605942188713761178189078584687;
					rDOFs1D[6]  = -0.29983046890076320809835345472230;
					rDOFs1D[7]  = -0.10132627352194944784303300504592;
					rDOFs1D[8]  =  0.10132627352194944784303300504592;
					rDOFs1D[9]  =  0.29983046890076320809835345472230;
					rDOFs1D[10] =  0.48605942188713761178189078584687;
					rDOFs1D[11] =  0.65238870288249308946788321964058;
					rDOFs1D[12] =  0.79200829186181506393108827096315;
					rDOFs1D[13] =  0.89920053309347209299462826151985;
					rDOFs1D[14] =  0.96956804627021793295224273836746;
					rDOFs1D[15] =  1.0;

					break;
				}

				case 17:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.97313217663141831415697950187372;
					rDOFs1D[2]  = -0.91087999591557359562380250639773;
					rDOFs1D[3]  = -0.81569625122177030710675055323753;
					rDOFs1D[4]  = -0.69102898062768470539491935737245;
					rDOFs1D[5]  = -0.54138539933010153912373340750406;
					rDOFs1D[6]  = -0.37217443356547704190723468073526;
					rDOFs1D[7]  = -0.18951197351831738830426301475311;
					rDOFs1D[8]  =  0.0;
					rDOFs1D[9]  =  0.18951197351831738830426301475311;
					rDOFs1D[10] =  0.37217443356547704190723468073526;
					rDOFs1D[11] =  0.54138539933010153912373340750406;
					rDOFs1D[12] =  0.69102898062768470539491935737245;
					rDOFs1D[13] =  0.81569625122177030710675055323753;
					rDOFs1D[14] =  0.91087999591557359562380250639773;
					rDOFs1D[15] =  0.97313217663141831415697950187372;
					rDOFs1D[16] =  1.0;

					break;
				}

				case 18:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.976105557412198542864518924341700;
					rDOFs1D[2]  = -0.920649185347533873837854625431280;
					rDOFs1D[3]  = -0.835593535218090213713646362327940;
					rDOFs1D[4]  = -0.723679329283242681306210365302070;
					rDOFs1D[5]  = -0.588504834318661761173535893193560;
					rDOFs1D[6]  = -0.434415036912123975342287136740670;
					rDOFs1D[7]  = -0.266362652878280984167665332025600;
					rDOFs1D[8]  = -0.089749093484652111022645010088562;
					rDOFs1D[9]  =  0.089749093484652111022645010088562;
					rDOFs1D[10] =  0.266362652878280984167665332025600;
					rDOFs1D[11] =  0.434415036912123975342287136740670;
					rDOFs1D[12] =  0.588504834318661761173535893193560;
					rDOFs1D[13] =  0.723679329283242681306210365302070;
					rDOFs1D[14] =  0.835593535218090213713646362327940;
					rDOFs1D[15] =  0.920649185347533873837854625431280;
					rDOFs1D[16] =  0.976105557412198542864518924341700;
					rDOFs1D[17] =  1.0;

					break;
				}

				case 19:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.97861176622208009515263406311022;
					rDOFs1D[2]  = -0.92890152815258624371794025879655;
					rDOFs1D[3]  = -0.85246057779664609308595597004106;
					rDOFs1D[4]  = -0.75149420255261301416363748963394;
					rDOFs1D[5]  = -0.62890813726522049776683230622873;
					rDOFs1D[6]  = -0.48822928568071350277790963762492;
					rDOFs1D[7]  = -0.33350484782449861029850010384493;
					rDOFs1D[8]  = -0.16918602340928157137515415344488;
					rDOFs1D[9]  =  0.0;
					rDOFs1D[10] =  0.16918602340928157137515415344488;
					rDOFs1D[11] =  0.33350484782449861029850010384493;
					rDOFs1D[12] =  0.48822928568071350277790963762492;
					rDOFs1D[13] =  0.62890813726522049776683230622873;
					rDOFs1D[14] =  0.75149420255261301416363748963394;
					rDOFs1D[15] =  0.85246057779664609308595597004106;
					rDOFs1D[16] =  0.92890152815258624371794025879655;
					rDOFs1D[17] =  0.97861176622208009515263406311022;
					rDOFs1D[18] =  1.0;

					break;
				}

				case 20:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.980743704893914171925446438584230;
					rDOFs1D[2]  = -0.935934498812665435716181584930630;
					rDOFs1D[3]  = -0.866877978089950141309847214616290;
					rDOFs1D[4]  = -0.775368260952055870414317527594690;
					rDOFs1D[5]  = -0.663776402290311289846403322971160;
					rDOFs1D[6]  = -0.534992864031886261648135961828980;
					rDOFs1D[7]  = -0.392353183713909299386474703815820;
					rDOFs1D[8]  = -0.239551705922986495182401356927090;
					rDOFs1D[9]  = -0.080545937238821837975944518159554;
					rDOFs1D[10] =  0.080545937238821837975944518159554;
					rDOFs1D[11] =  0.239551705922986495182401356927090;
					rDOFs1D[12] =  0.392353183713909299386474703815820;
					rDOFs1D[13] =  0.534992864031886261648135961828980;
					rDOFs1D[14] =  0.663776402290311289846403322971160;
					rDOFs1D[15] =  0.775368260952055870414317527594690;
					rDOFs1D[16] =  0.866877978089950141309847214616290;
					rDOFs1D[17] =  0.935934498812665435716181584930630;
					rDOFs1D[18] =  0.980743704893914171925446438584230;
					rDOFs1D[19] =  1.0;

					break;
				}

				case 21:
				{
					rDOFs1D[0]  = -1.0;
					rDOFs1D[1]  = -0.98257229660454802823448127655541;
					rDOFs1D[2]  = -0.94197629695974553429610265066144;
					rDOFs1D[3]  = -0.87929475532359046445115359630494;
					rDOFs1D[4]  = -0.79600192607771240474431258966036;
					rDOFs1D[5]  = -0.69405102606222323262731639319467;
					rDOFs1D[6]  = -0.57583196026183068692702187033809;
					rDOFs1D[7]  = -0.44411578327900210119451634960735;
					rDOFs1D[8]  = -0.30198985650876488727535186785875;
					rDOFs1D[9]  = -0.15278551580218546600635832848567;
					rDOFs1D[10] =  0.0;
					rDOFs1D[11] =  0.15278551580218546600635832848567;
					rDOFs1D[12] =  0.30198985650876488727535186785875;
					rDOFs1D[13] =  0.44411578327900210119451634960735;
					rDOFs1D[14] =  0.57583196026183068692702187033809;
					rDOFs1D[15] =  0.69405102606222323262731639319467;
					rDOFs1D[16] =  0.79600192607771240474431258966036;
					rDOFs1D[17] =  0.87929475532359046445115359630494;
					rDOFs1D[18] =  0.94197629695974553429610265066144;
					rDOFs1D[19] =  0.98257229660454802823448127655541;
					rDOFs1D[20] =  1.0;

					break;
				}

				default: ERROR("LGL distribution for given polynomial order not (yet) implemented!");
			}

			break;
		}

		default: ERROR("Unknown type of DOF distribution!");
	}

	// Make sure DOFs are correct.
	as3double tmp = 0.0;
	for(unsigned int i=0; i<nDOFs1D; i++) tmp += rDOFs1D[i];
	if( fabs(tmp) > TOL ) ERROR("Precision errors in nodal point computation.");
}


