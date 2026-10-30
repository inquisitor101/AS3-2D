#include "test_standard_element.hpp"


//-----------------------------------------------------------------------------------
// CTest_SE member functions.
//-----------------------------------------------------------------------------------


void CTest_SE::SetUp
(
 void
)
 /*
	* Function that serves as a constructor for the object.
	*/
{
	// Suppress cout
	auto buffer = std::cout.rdbuf(nullptr);

	mConfig  = std::make_unique<CConfig>("param_t1.cfg");
	ASSERT_GT(mConfig->GetnZone(), 0);

	mElement.resize(mConfig->GetnZone());
	for(unsigned short iZone=0; iZone<mConfig->GetnZone(); iZone++) 
	{
		mElement[iZone] = CGenericFactory::CreateStandardElement(mConfig.get(), iZone);
	}

	// Restore cout.
	std::cout.rdbuf(buffer);
}

//-----------------------------------------------------------------------------------

as3double CTest_SE::ComputeLagrangePolynomial
(
 as3double               x,
 std::vector<as3double> &r,
 size_t                  i
)
 /*
	* Function that evaluates a single entry (x) on a
	* given basis (r) using the Lagrange polynomial at
	* degree (i).
	*/
{
  // Number of basis points.
  size_t n = r.size();
	// Temporary value to store results in.
	as3double ell = 1.0;
	for(size_t j=0; j<n; j++)
		if( j != i ) ell *= (x - r[j])/(r[i]-r[j]);

	return ell;
}

//-----------------------------------------------------------------------------------

void CTest_SE::CheckBasisLocation1D
(
 CStandardElement *element
)
	/*
	 * Function that ensures the DOF location is correct.
	 */
{
	CMatrixAS3<as3double> A = element->GetrSol1D();

	CMatrixAS3<as3double> A0; 
	ComputeLocationDOFs1D(element->GetnSol1D(),
			                  element->GetTypeDOFsSol(),
												A0);

	NTest_LA::MatrixErrorNormLinf(A, A0, 1.0e-12);
}

//-----------------------------------------------------------------------------------

void CTest_SE::CheckQuadrature
(
 CStandardElement *element
)
	/*
	 * Function that ensures the quadrature nodes and weights are correct.
	 */
{
	CMatrixAS3<as3double> A = element->GetrInt1D();
	CMatrixAS3<as3double> B = element->GetwInt1D();

	CMatrixAS3<as3double> A0;
	CMatrixAS3<as3double> B0;
	ComputeQuadrature1D(element->GetnInt1D(), A0, B0);

	NTest_LA::MatrixErrorNormLinf(A, A0, 1.0e-12);
	NTest_LA::MatrixErrorNormLinf(B, B0, 1.0e-12);

	CMatrixAS3<as3double> C = element->GetwInt2D();

	CMatrixAS3<as3double> C0;

	CTest_TP::KroneckerProduct(B0, B0, C0);

	NTest_LA::MatrixErrorNormLinf(C, C0, 1.0e-12);
}

//-----------------------------------------------------------------------------------

void CTest_SE::ComputeLocationDOFs1D
(
 size_t                 nDOFs1D, 
 ETypeDOF               typeDOFs, 
 CMatrixAS3<as3double> &rDOFs1D
)
	/*
	 * Function that computes the solution DOFs in 1D on a standard element in [-1,+1].
	 * Note, these are based on the values used by Edwin van der Weide in VCP3D.
	 */
{
  // Allocate the memory for rDOFs1D.
  rDOFs1D.resize(nDOFs1D);

  // Make a distinction between the possibilities for the DOF locations.
  switch( typeDOFs )
  {
		case ETypeDOF::EQD:
    {
      // Equidistant spacing is used.
      const as3double dh = 2.0/(nDOFs1D-1);

      for(int i=0; i<nDOFs1D; ++i)
        rDOFs1D[i] = -1.0 + i*dh;

      break;
    }

    //--------------------------------------------------------------------------

		case ETypeDOF::LGL:
    {
      switch( nDOFs1D )
      {
        case 1:
        {
          rDOFs1D[0] = 0.0;
          break;
        }

        case 2:
        {
          rDOFs1D[0] = -1.0; rDOFs1D[1] = 1.0;
          break;
        }

        case 3:
        {
          rDOFs1D[0] = -1.0; rDOFs1D[1] = 0.0; rDOFs1D[2] = 1.0;
          break;
        }

        case 4:
        {
          rDOFs1D[0] = -1.0;           rDOFs1D[1] = -std::sqrt(1.0/5.0);
          rDOFs1D[2] = std::sqrt(1.0/5.0); rDOFs1D[3] = 1.0;
          break;
        }

        case 5:
        {
          rDOFs1D[0] = -1.0;              rDOFs1D[1] = -std::sqrt(3.0/7.0); rDOFs1D[2] = 0.0;
          rDOFs1D[3] = std::sqrt(3.0/7.0); rDOFs1D[4] = 1.0;
          break;
        }

        case 6:
        {
          const as3double t1 = 2.0*std::sqrt(7.0)/21.0;
          const as3double t2 = std::sqrt(1.0/3.0 + t1);
          const as3double t3 = std::sqrt(1.0/3.0 - t1);
          rDOFs1D[0] = -1.0; rDOFs1D[1] = -t2; rDOFs1D[2] = -t3;
          rDOFs1D[3] =  t3;  rDOFs1D[4] =  t2; rDOFs1D[5] =  1.0;
          break;
        }

        case 7:
        {
          const as3double t1 = 2.0*std::sqrt(5.0/3.0)/11.0;
          const as3double t2 = std::sqrt(5.0/11.0 + t1);
          const as3double t3 = std::sqrt(5.0/11.0 - t1);
          rDOFs1D[0] = -1.0; rDOFs1D[1] = -t2; rDOFs1D[2] = -t3;  rDOFs1D[3] = 0.0;
          rDOFs1D[4] =  t3;  rDOFs1D[5] =  t2; rDOFs1D[6] =  1.0;
          break;
        }

        case 8:
        {
          rDOFs1D[0] = -1.0;
          rDOFs1D[1] = (as3double) -0.8717401485096066153375;
          rDOFs1D[2] = (as3double) -0.5917001814331423021445;
          rDOFs1D[3] = (as3double) -0.2092992179024788687687;
          rDOFs1D[4] = (as3double)  0.2092992179024788687687;
          rDOFs1D[5] = (as3double)  0.5917001814331423021445;
          rDOFs1D[6] = (as3double)  0.8717401485096066153375;
          rDOFs1D[7] =  1.0;
          break;
        }

        case 9:
        {
          rDOFs1D[0] = -1.0;
          rDOFs1D[1] = (as3double) -0.8997579954114601573124;
          rDOFs1D[2] = (as3double) -0.6771862795107377534459;
          rDOFs1D[3] = (as3double) -0.3631174638261781587108;
          rDOFs1D[4] =  0.0;
          rDOFs1D[5] = (as3double)  0.3631174638261781587108;
          rDOFs1D[6] = (as3double)  0.6771862795107377534459;
          rDOFs1D[7] = (as3double)  0.8997579954114601573124;
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

        default: FAIL();
      }

      break;
    }

    //--------------------------------------------------------------------------

    default: FAIL();
  }
}

//-----------------------------------------------------------------------------------

void CTest_SE::ComputeQuadrature1D
(
 size_t                   nInt1D,
 CMatrixAS3<as3double>   &rInt1D,
 CMatrixAS3<as3double>   &wInt1D
)
	/*
	 * Function that computes the quadrature points and weights in 1D on a standard element in [-1,+1].
	 * Note, these have been computed via: 
	 * https://nl.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
	 */
{
  // Allocate the memory.
  rInt1D.resize(nInt1D);
	wInt1D.resize(nInt1D);

	// Check which set of coefficients for GL quadrature to use.
	switch(nInt1D)
	{
		case(2):
		{
			rInt1D[0] = -0.57735026918962576450914878050196;
    	rInt1D[1] =  0.57735026918962576450914878050196;
		
			wInt1D[0] =  1.0;
      wInt1D[1] =  1.0;

			break;
		}

		case(4):
		{
			rInt1D[0] = -0.86113631159405246151550272770692;
		  rInt1D[1] = -0.33998104358485625731134405214107;
		  rInt1D[2] =  0.33998104358485625731134405214107;
		  rInt1D[3] =  0.86113631159405246151550272770692;
		
			wInt1D[0] =  0.34785484513745440482423987305083;
			wInt1D[1] =  0.65214515486254598375381874575396;
			wInt1D[2] =  0.65214515486254598375381874575396;
			wInt1D[3] =  0.34785484513745440482423987305083;

			break;
		}

		case(5):
		{
			rInt1D[0] = -0.90617984593866385267801888403483;
			rInt1D[1] = -0.53846931010568299669216685288120;
			rInt1D[2] =  0.0;
			rInt1D[3] =  0.53846931010568299669216685288120;
			rInt1D[4] =  0.90617984593866385267801888403483;
		
			wInt1D[0] =  0.23692688505618911265493409246119;
			wInt1D[1] =  0.47862867049936624885830838138645;
			wInt1D[2] =  0.56888888888888888888888888888889;
			wInt1D[3] =  0.47862867049936624885830838138645;
			wInt1D[4] =  0.23692688505618911265493409246119;

			break;
		}

		case(7):
		{
			rInt1D[0] =  -0.9491079123427583752459213428665;
			rInt1D[1] =  -0.7415311855993944600839995473506;
			rInt1D[2] =	 -0.4058451513773971286447306283662;
			rInt1D[3] =	  0.0;
			rInt1D[4] =		0.4058451513773971286447306283662;
			rInt1D[5] =		0.7415311855993944600839995473506;
			rInt1D[6] =		0.9491079123427583752459213428665;
     
		
			wInt1D[0] =  0.12948496616886964738490917170566;
			wInt1D[1] =  0.27970539148927675565659001222230;
			wInt1D[2] =	 0.38183005050511903410992431417981;
			wInt1D[3] =	 0.41795918367346938775510204081633;
			wInt1D[4] =	 0.38183005050511903410992431417981;
			wInt1D[5] =	 0.27970539148927675565659001222230;
			wInt1D[6] =	 0.12948496616886964738490917170566;

			break;
		}

		case(8):
		{
			rInt1D[0] = -0.96028985649753639819437012192793; 
			rInt1D[1] = -0.79666647741362672796583410672611;
			rInt1D[2] = -0.52553240991632899081764662696514;
			rInt1D[3] = -0.18343464249564983559181996497500;
			rInt1D[4] =  0.18343464249564983559181996497500;
			rInt1D[5] =  0.52553240991632899081764662696514;
			rInt1D[6] =  0.79666647741362672796583410672611;
			rInt1D[7] =  0.96028985649753639819437012192793; 
		
			wInt1D[0] =  0.10122853629037681377766944024188;
			wInt1D[1] =  0.22238103445337442654050619239570;
			wInt1D[2] =  0.31370664587788732458051299545332;
			wInt1D[3] =  0.36268378337836210123512614700303;
			wInt1D[4] =  0.36268378337836210123512614700303;
			wInt1D[5] =  0.31370664587788732458051299545332;
			wInt1D[6] =  0.22238103445337442654050619239570;
			wInt1D[7] =  0.10122853629037681377766944024188;

			break;
		}

		case(10):
		{
			rInt1D[0] = -0.97390652851717174343093574861996; 
			rInt1D[1] = -0.86506336668898464736798814556096;
			rInt1D[2] = -0.67940956829902443558921731892042;
			rInt1D[3] = -0.43339539412924721339948064269265;
			rInt1D[4] = -0.14887433898163116019475182838505;
			rInt1D[5] =  0.14887433898163116019475182838505;
			rInt1D[6] =  0.43339539412924721339948064269265;
			rInt1D[7] =  0.67940956829902443558921731892042;
			rInt1D[8] =  0.86506336668898464736798814556096;
			rInt1D[9] =  0.97390652851717174343093574861996; 
		
			wInt1D[0] =  0.06667134430868802696945607522138;
			wInt1D[1] =  0.14945134915058055913306134243612;
			wInt1D[2] =  0.21908636251598206934332324635761;
			wInt1D[3] =  0.26926671930999618309598986343190;
			wInt1D[4] =  0.29552422471475292553577673970722;
			wInt1D[5] =  0.29552422471475292553577673970722;
			wInt1D[6] =  0.26926671930999618309598986343190;
			wInt1D[7] =  0.21908636251598206934332324635761;
			wInt1D[8] =  0.14945134915058055913306134243612;
			wInt1D[9] =  0.06667134430868802696945607522138; 
	
			break;
		}

		case(11):
		{
			rInt1D[ 0] = -0.97822865814605686196614442451391;  
			rInt1D[ 1] = -0.88706259976809542777687056513969;
			rInt1D[ 2] = -0.73015200557404935644001398031833;
			rInt1D[ 3] = -0.51909612920681169612180383410305;
			rInt1D[ 4] = -0.26954315595234490388065751176327;
			rInt1D[ 5] =  0.0;
			rInt1D[ 6] =  0.26954315595234490388065751176327;
			rInt1D[ 7] =  0.51909612920681169612180383410305;
			rInt1D[ 8] =  0.73015200557404935644001398031833;
			rInt1D[ 9] =  0.88706259976809542777687056513969;
			rInt1D[10] =  0.97822865814605686196614442451391; 

			wInt1D[ 0] =  0.05566856711617353820065190461718;
			wInt1D[ 1] =  0.12558036946490408469756516751659;
			wInt1D[ 2] =	0.18629021092773401235831443045754;
			wInt1D[ 3] =	0.23319376459199023243762383117428;
			wInt1D[ 4] =	0.26280454451024670703418451012112;
			wInt1D[ 5] =	0.27292508677790061621948325409903;
			wInt1D[ 6] =	0.26280454451024670703418451012112;
			wInt1D[ 7] =	0.23319376459199023243762383117428;
			wInt1D[ 8] =	0.18629021092773401235831443045754;
			wInt1D[ 9] =	0.12558036946490408469756516751659;
			wInt1D[10] =	0.05566856711617353820065190461718; 

			break;
		}

		case(13):
		{
			rInt1D[ 0] = -0.98418305471858813504582030873280; 
			rInt1D[ 1] = -0.91759839922297792291772111639148;
			rInt1D[ 2] = -0.80157809073330987814642867306247;
			rInt1D[ 3] = -0.64234933944034011688017926644534;
			rInt1D[ 4] = -0.44849275103644692386239967163419;
			rInt1D[ 5] = -0.23045831595513477374481681181351;
			rInt1D[ 6] =  0.0;
			rInt1D[ 7] =  0.23045831595513477374481681181351;
			rInt1D[ 8] =  0.44849275103644692386239967163419;
			rInt1D[ 9] =  0.64234933944034011688017926644534;
			rInt1D[10] =  0.80157809073330987814642867306247;
			rInt1D[11] =  0.91759839922297792291772111639148;
			rInt1D[12] =  0.98418305471858813504582030873280;  

			wInt1D[ 0] =  0.04048400476531581471117959836192; 
			wInt1D[ 1] =  0.09212149983772843775398087018402;
			wInt1D[ 2] =  0.13887351021978713849769349053531;
			wInt1D[ 3] =  0.17814598076194568254670969054132;
			wInt1D[ 4] =  0.20781604753688834308356092606118;
			wInt1D[ 5] =  0.22628318026289709341547506937786;
			wInt1D[ 6] =  0.23255155323087389751535170034913;
			wInt1D[ 7] =  0.22628318026289709341547506937786;
			wInt1D[ 8] =  0.20781604753688834308356092606118;
			wInt1D[ 9] =  0.17814598076194568254670969054132;
			wInt1D[10] =  0.13887351021978713849769349053531;
			wInt1D[11] =  0.09212149983772843775398087018402;
			wInt1D[12] =  0.04048400476531581471117959836192; 

			break;
		}

		case(14):
		{
			rInt1D[ 0] = -0.98628380869681242515412122884300;
			rInt1D[ 1] = -0.92843488366357362906455819029361;
			rInt1D[ 2] = -0.82720131506976501967187687114347;
			rInt1D[ 3] = -0.68729290481168536786071854294278;
			rInt1D[ 4] = -0.51524863635815409956819621584145;
			rInt1D[ 5] = -0.31911236892788968910750213581196;
			rInt1D[ 6] = -0.10805494870734366763542766420869;
			rInt1D[ 7] =  0.10805494870734366763542766420869;
			rInt1D[ 8] =  0.31911236892788968910750213581196;
			rInt1D[ 9] =  0.51524863635815409956819621584145;
			rInt1D[10] =  0.68729290481168536786071854294278;
			rInt1D[11] =  0.82720131506976501967187687114347;
			rInt1D[12] =  0.92843488366357362906455819029361;
			rInt1D[13] =  0.98628380869681242515412122884300; 

			wInt1D[ 0] =  0.03511946033175199904929897343208;
			wInt1D[ 1] =  0.08015808715976004139580624041627;
			wInt1D[ 2] =  0.12151857068790319904572072573501;
			wInt1D[ 3] =  0.15720316715819357411554335612891;
			wInt1D[ 4] =  0.18553839747793784975549158389185;
			wInt1D[ 5] =  0.20519846372129568745634742299444;
			wInt1D[ 6] =  0.21526385346315768387626121693756;
			wInt1D[ 7] =  0.21526385346315768387626121693756;
			wInt1D[ 8] =  0.20519846372129568745634742299444;
			wInt1D[ 9] =  0.18553839747793784975549158389185;
			wInt1D[10] =  0.15720316715819357411554335612891;
			wInt1D[11] =  0.12151857068790319904572072573501;
			wInt1D[12] =  0.08015808715976004139580624041627;
			wInt1D[13] =  0.03511946033175199904929897343208; 

			break;
		}

		//--------------------------------------------------------------------------
		
		default: FAIL();
	}
}







