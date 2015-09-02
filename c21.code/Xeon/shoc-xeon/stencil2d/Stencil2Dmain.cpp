#include <iostream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include "omp.h"

#include "OptionParser.h"
#include "ResultDatabase.h"
#include "Timer.h"
#include "BadCommandLine.h"
#include "InvalidArgValue.h"
#include "Matrix2D.h"
#include "HostStencilFactory.h"
#include "HostStencil.h"
#include "MICStencilFactory.h"
#include "MICStencil.h"
#include "InitializeMatrix2D.h"
#include "ValidateMatrix2D.h"
#include "StencilUtil.h"

#include "InitializeMatrix2D.cpp"
#include "ValidateMatrix2D.cpp"
#include "StencilUtil.cpp"
#include "StencilFactory.cpp"
#include "CommonMICStencilFactory.cpp"
#include "HostStencil.cpp"
#include "MICStencil.cpp"
#include "HostStencilFactory.cpp"
#include "MICStencilFactory.cpp"

#include <sys/mman.h>

// prototypes of auxiliary functions defined in this file or elsewhere
void CheckOptions( const OptionParser& opts );
void EnsureStencilInstantiation( void );

#define LINESIZE 64

template<class T> void
MICValidate(const Matrix2D<T>& s, const Matrix2D<T>& t, double valErrThreshold, unsigned int nValErrsToPrint)
{
    assert( (s.GetNumRows() == t.GetNumRows()) && (s.GetNumColumns() == t.GetNumColumns()) );
    unsigned int uHaloWidth = LINESIZE / sizeof(T);

	#if 1

		for( unsigned int i = uHaloWidth; i < s.GetNumRows() - uHaloWidth; i++ )
		{
			for( unsigned int j = uHaloWidth; j < s.GetNumColumns() - uHaloWidth; j++ )
			{
				T expVal 	= s.GetConstData()[i][j];
				T actualVal = t.GetConstData()[i][j];
				T delta 	= fabsf( actualVal - expVal );
				T relError 	= (expVal != 0.0f) ? delta / expVal : 0.0f;

				if( relError > valErrThreshold )
				{
					// printf("i = %d, j = %d\n", i, j);
					std::cout<<"Failed\n";
					return;
				}
			}
		}

		std::cout<<"Passed\n";
	#endif // if 1

	#if 0
		std::cout<<"Expected Value \n";
		for( unsigned int i = 0; i < s.GetNumRows(); i++ )
		{
			for( unsigned int j = 0; j < s.GetNumColumns(); j++ )
			{
				T expVal = s.GetConstData()[i][j];
				std::cout<<expVal<<" ";
			}
			std::cout<<endl;
		}

		std::cout<<"Calculated vaue \n";
		for( unsigned int i = 0; i < s.GetNumRows(); i++ )
		{
			for( unsigned int j = 0; j < s.GetNumColumns(); j++ )
			{
				T expVal = t.GetConstData()[i][j];
				std::cout<<expVal<<" ";
			}
			std::cout<<endl;
		}
	#endif // if 0
}

template<class T> void
DoTest( const char* timerDesc, ResultDatabase& resultDB, OptionParser& opts )
{
    StencilFactory<T>* 	stdStencilFactory 	= NULL;
    Stencil<T>* 		stdStencil 			= NULL;
    StencilFactory<T>* 	testStencilFactory 	= NULL;
    Stencil<T>* 		testStencil 		= NULL;

	stdStencilFactory 	= new HostStencilFactory<T>;
    testStencilFactory 	= new MICStencilFactory<T>;
    assert( (stdStencilFactory != NULL) && (testStencilFactory != NULL) );

    // Do a sanity check on option values
    CheckOptions( opts );
    stdStencilFactory->CheckOptions( opts );
    testStencilFactory->CheckOptions( opts );

    // Extract and validate options
    std::vector<long long> arrayDims = opts.getOptionVecInt( "customSize" );
    if( arrayDims.size() != 2 )
    	cerr << "Dim size: " << arrayDims.size() << "\n";

	if (arrayDims[0] == 0) // User has not specified a custom size
	{
		const int probSizes[4] = { 768, 1408, 2048, 4096 };
		int sizeClass = opts.getOptionInt("size");

		if (!(sizeClass >= 0 && sizeClass < 5))
        {
			//throw InvalidArgValue( "Size class must be between 1-4" );
        }

		arrayDims[0] = arrayDims[1] = probSizes[sizeClass - 1];
	}

	long int seed 					= (long)opts.getOptionInt( "seed" );
	bool beVerbose 					= opts.getOptionBool( "verbose" );
	unsigned int nIters 			= (unsigned int)opts.getOptionInt( "num-iters" );
	double valErrThreshold 			= (double)opts.getOptionFloat( "val-threshold" );
	unsigned int nValErrsToPrint 	= (unsigned int)opts.getOptionInt( "val-print-limit" );

	// Define Halo
	unsigned int haloWidth 	= LINESIZE / sizeof(T);
	float haloVal 			= (float)opts.getOptionFloat( "haloVal" );

	// Build a description of this experiment
	std::ostringstream experimentDescriptionStr;
	experimentDescriptionStr<< nIters << ':'<< arrayDims[0] << 'x' << arrayDims[1];

	unsigned int 	nPasses = (unsigned int)opts.getOptionInt( "passes" );
	unsigned long 	npts 	= (arrayDims[0] + 2 * haloWidth - 2) * (arrayDims[1] + 2*haloWidth - 2);
	unsigned long 	nflops 	= npts * 11 * nIters;
	cout<<"FLOP are = "<< nflops <<endl;

	Matrix2D<T> 	exp(arrayDims[0] + 2 * haloWidth, arrayDims[1] + 2 * haloWidth);
	Initialize<T> 	init(seed, haloWidth, haloVal);

	init(exp);
	if(beVerbose)
		std::cout << "initial state:\n" << exp << std::endl;

	stdStencil = stdStencilFactory->BuildStencil(opts);

	(*stdStencil)(exp, nIters);

	if( beVerbose )
		std::cout << "expected result:\n" << exp << std::endl;

	// Compute the result on the Xeon Phi device
	Matrix2D<T> data(arrayDims[0] + 2 * haloWidth, arrayDims[1] + 2 * haloWidth);
	testStencil = testStencilFactory->BuildStencil( opts );

	//printf("1.pIn[1300] = %f \n", data.GetFlatData()[1300]);
	std::cout<<"Passes:"<<nPasses<<endl;
	for( unsigned int pass = 0; pass < nPasses; pass++ )
	{
		init(data);

		double start 		= curr_second();
		(*testStencil)(data, nIters);
		double elapsedTime 	= curr_second() - start;

		double gflopsPCIe 	= (nflops / elapsedTime) / 1e9;

		resultDB.AddResult(timerDesc, experimentDescriptionStr.str(), "GFLOPS_PCIe", gflopsPCIe);

		if( beVerbose )
			std::cout << "observed result, pass " << pass << ":\n"<< data<< std::endl;

		MICValidate(exp, data, valErrThreshold, nValErrsToPrint);
	}

    // clean up - normal termination
    delete stdStencil;
    delete stdStencilFactory;
    delete testStencil;
    delete testStencilFactory;
}

void RunBenchmark(OptionParser& opts, ResultDatabase& resultDB )
{
//    std::cout << mlockall(MCL_FUTURE) << std::endl;
       std::cout << "Running Single Precision test :" << std::endl;
       DoTest<float>( "SP_Sten2D", resultDB, opts);

       std::cout << "Running Double Precision test :" << std::endl;
       DoTest<double>( "DP_Sten2D", resultDB, opts );
//   std::cout << munlockall() << std::endl;
}

// Adds command line options to given OptionParser
void addBenchmarkSpecOptions( OptionParser& opts )
{
    opts.addOption( "customSize", 		OPT_VECINT, "0,0", 		"specify custom problem size");
    opts.addOption( "num-iters", 		OPT_INT, 	"1000", 	"number of stencil iterations" );
    opts.addOption( "weight-center", 	OPT_FLOAT, 	"0.25", 	"center value weight" );
    opts.addOption( "weight-cardinal", 	OPT_FLOAT, 	"0.15", 	"cardinal values weight" );
    opts.addOption( "weight-diagonal", 	OPT_FLOAT, 	"0.05", 	"diagonal values weight" );
    opts.addOption( "seed", 			OPT_INT, 	"71594", 	"random number generator seed" );
    opts.addOption( "val-threshold", 	OPT_FLOAT, 	"0.1", 	"validation error threshold" );
    opts.addOption( "val-print-limit", 	OPT_INT, 	"15", 		"number of validation errors to print" );
    opts.addOption( "haloVal", 			OPT_FLOAT, 	"0.0", 		"value to use for halo data" );
}

// validate stencil-independent values
void CheckOptions( const OptionParser& opts )
{
    // check matrix dimensions - must be 2d, must be positive
    std::vector<long long> arrayDims = opts.getOptionVecInt( "customSize" );
    if( arrayDims.size() != 2 )
    {
        throw InvalidArgValue( "overall size must have two dimensions" );
    }

    if( (arrayDims[0] < 0) || (arrayDims[1] < 0) )
    {
        throw InvalidArgValue( "each size dimension must be positive" );
    }

    // validation error threshold must be positive
    float valThreshold = opts.getOptionFloat( "val-threshold" );
    if( valThreshold <= 0.0f )
    {
        throw InvalidArgValue( "validation threshold must be positive" );
    }

    // number of validation errors to print must be non-negative
    int nErrsToPrint = opts.getOptionInt( "val-print-limit" );
    if( nErrsToPrint < 0 )
    {
        throw InvalidArgValue( "number of validation errors to print must be non-negative" );
    }
}
