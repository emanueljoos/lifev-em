/* -*- Mode : c++; c-tab-always-indent: t; indent-tabs-mode: nil; -*-

  <short description here>

  Gilles Fourestey gilles.fourestey@epfl.ch

*/
/** \file DataFSI.cpp

*/


#include <life/lifesolver/DataFSI.hpp>


namespace LifeV
{


// ===================================================
// Constructors
// ===================================================

DataFSI::DataFSI( ) :
    M_dataFluid                     ( new dataFluid_Type() ),
    M_dataSolid                     ( new dataSolid_Type() ),
    M_maxSubIterationNumber         (),
    M_absoluteTolerance             (),
    M_relativeTolerance             (),
    M_errorTolerance                (),
    M_linesearch                    (),
    M_preconditioner                (),
    M_DDNpreconditioner             (),
    M_DDBlockPreconditioner         (),
    M_method                        (),
    M_algorithm                     (),
    M_defaultOmega                  (),
    M_updateEvery                   (),
    M_fluidInterfaceFlag            (),
    M_solidInterfaceFlag            (),
    M_structureInterfaceFlag        (),
    M_harmonicInterfaceFlag         (),
    M_interfaceTolerance            (),
    M_RobinNeumannCoupling          (),
    M_RobinNeumannFluidCoefficient  (),
    M_RobinNeumannSolidCoefficient  ()
{
}


DataFSI::DataFSI( const DataFSI& DataFSI ) :
    M_dataFluid                     ( DataFSI.M_dataFluid ),
    M_dataSolid                     ( DataFSI.M_dataSolid ),
    M_maxSubIterationNumber         ( DataFSI.M_maxSubIterationNumber ),
    M_absoluteTolerance             ( DataFSI.M_absoluteTolerance ),
    M_relativeTolerance             ( DataFSI.M_relativeTolerance ),
    M_errorTolerance                ( DataFSI.M_errorTolerance ),
    M_linesearch                    ( DataFSI.M_linesearch ),
    M_preconditioner                ( DataFSI.M_preconditioner ),
    M_DDNpreconditioner             ( DataFSI.M_DDNpreconditioner ),
    M_DDBlockPreconditioner         ( DataFSI.M_DDBlockPreconditioner ),
    M_method                        ( DataFSI.M_method ),
    M_algorithm                     ( DataFSI.M_algorithm ),
    M_defaultOmega                  ( DataFSI.M_defaultOmega ),
    M_updateEvery                   ( DataFSI.M_updateEvery ),
    M_fluidInterfaceFlag            ( DataFSI.M_fluidInterfaceFlag ),
    M_solidInterfaceFlag            ( DataFSI.M_solidInterfaceFlag ),
    M_structureInterfaceFlag        ( DataFSI.M_structureInterfaceFlag ),
    M_harmonicInterfaceFlag         ( DataFSI.M_harmonicInterfaceFlag ),
    M_interfaceTolerance            ( DataFSI.M_interfaceTolerance ),
    M_RobinNeumannCoupling          ( DataFSI.M_RobinNeumannCoupling ),
    M_RobinNeumannFluidCoefficient  ( DataFSI.M_RobinNeumannFluidCoefficient ),
    M_RobinNeumannSolidCoefficient  ( DataFSI.M_RobinNeumannSolidCoefficient )
{
}


// ===================================================
// Methods
// ===================================================


DataFSI&
DataFSI::operator=( const DataFSI& DataFSI )
{
    if ( this != &DataFSI )
    {
        M_dataFluid                     = DataFSI.M_dataFluid;
        M_dataSolid                     = DataFSI.M_dataSolid;
        M_maxSubIterationNumber         = DataFSI.M_maxSubIterationNumber;
        M_absoluteTolerance             = DataFSI.M_absoluteTolerance;
        M_relativeTolerance             = DataFSI.M_relativeTolerance;
        M_errorTolerance                = DataFSI.M_errorTolerance;
        M_linesearch                    = DataFSI.M_linesearch;
        M_preconditioner                = DataFSI.M_preconditioner;
        M_DDNpreconditioner             = DataFSI.M_DDNpreconditioner;
        M_DDBlockPreconditioner         = DataFSI.M_DDBlockPreconditioner;
        M_method                        = DataFSI.M_method;
        M_algorithm                     = DataFSI.M_algorithm;
        M_defaultOmega                  = DataFSI.M_defaultOmega;
        M_updateEvery                   = DataFSI.M_updateEvery;
        M_fluidInterfaceFlag            = DataFSI.M_fluidInterfaceFlag;
        M_solidInterfaceFlag            = DataFSI.M_solidInterfaceFlag;
        M_structureInterfaceFlag        = DataFSI.M_structureInterfaceFlag;
        M_harmonicInterfaceFlag         = DataFSI.M_harmonicInterfaceFlag;
        M_interfaceTolerance            = DataFSI.M_interfaceTolerance;
        M_RobinNeumannCoupling          = DataFSI.M_RobinNeumannCoupling;
        M_RobinNeumannFluidCoefficient  = DataFSI.M_RobinNeumannFluidCoefficient;
        M_RobinNeumannSolidCoefficient  = DataFSI.M_RobinNeumannSolidCoefficient;
    }

	return *this;
}


void
DataFSI::setup( const GetPot& dataFile, const std::string& section )
{
    M_dataFluid->setup( dataFile );
    M_dataSolid->setup( dataFile );

    // Problem - Non Linear Richardson Parameters
    M_maxSubIterationNumber = dataFile( ( section + "/maxSubIter" ).data(), 300 );
    M_absoluteTolerance = dataFile( ( section + "/abstol" ).data(), 1.e-07 );
    M_relativeTolerance = dataFile( ( section + "/reltol" ).data(), 1.e-04 );
    M_errorTolerance = dataFile( ( section + "/etamax" ).data(), 1.e-03 );
    M_linesearch = static_cast<Int> ( dataFile( ( section + "/linesearch" ).data(), 0 ) );

    // Problem - Preconditioner
    M_preconditioner = static_cast<Preconditioner> ( dataFile( ( section + "/precond" ).data(), DIRICHLET_NEUMANN ) );
    M_DDNpreconditioner = static_cast<DDNPreconditioner> ( dataFile( ( section + "/DDNprecond" ).data(), DDN_DIRICHLET_NEUMANN ) );

    // Problem - Methods
    M_method = dataFile( ( section + "/method" ).data(), "steklovPoincare" );
    M_algorithm = dataFile( ( section + "/algorithm" ).data(), "DirichletNeumann" );

    // Problem - FixPoint / EJ
    M_defaultOmega = dataFile( ( section + "/defOmega" ).data(), 0.001);
    M_updateEvery = dataFile( ( section + "/updateEvery" ).data(), 1);

    // Interface
    M_fluidInterfaceFlag = dataFile( "interface/fluid_flag",     1 );
    M_solidInterfaceFlag = dataFile( "interface/solid_flag",     M_fluidInterfaceFlag );
    M_structureInterfaceFlag = dataFile( "interface/structure_flag", M_fluidInterfaceFlag );
    M_harmonicInterfaceFlag = dataFile( "interface/harmonic_flag",  M_fluidInterfaceFlag );
    M_interfaceTolerance = dataFile( "interface/tolerance",      0. );

    // Interface - Monolithic
    M_DDBlockPreconditioner = dataFile( "interface/DDBlockPrec", 0 );
    M_RobinNeumannCoupling  = dataFile( "interface/robinNeumannCoupling", false );
    M_RobinNeumannFluidCoefficient = dataFile( "interface/alphaf", 0.5 );
    M_RobinNeumannSolidCoefficient = dataFile( "interface/alphas", 0.5 );
}


bool
DataFSI::isMonolithic()
{
    return !( M_method.compare( "monolithic" ) && M_method.compare( "fullMonolithic" ) );
}


void
DataFSI::showMe( std::ostream& output )
{
    output << "\n*** Values for data fluid\n\n";
    M_dataFluid->showMe();

    output << "\n*** Values for data solid\n\n";
    M_dataSolid->showMe();

    output << "\n*** Values for problem\n\n";
    output << "Max subiteration number          = " << M_maxSubIterationNumber << std::endl;
    output << "Absolute tolerance               = " << M_absoluteTolerance << std::endl;
    output << "Relative tolerance               = " << M_relativeTolerance << std::endl;
    output << "Max error tolerance              = " << M_errorTolerance << std::endl;
    output << "Linesearch                       = " << M_linesearch << std::endl;

    output << "Preconditioner                   = " << M_preconditioner << std::endl;
    output << "DDNPreconditioner                = " << M_DDNpreconditioner << std::endl;
    output << "DDBlockPreconditioner            = " << M_DDBlockPreconditioner << std::endl;

    output << "Method                           = " << M_method << std::endl;
    output << "Algorithm                        = " << M_algorithm << std::endl;

    output << "Default Omega                    = " << M_defaultOmega << std::endl;
    output << "Update every                     = " << M_updateEvery << std::endl;

    output << "\n*** Values for interface\n\n";
    output << "Interface fluid                  = " << M_fluidInterfaceFlag << std::endl;
    output << "Interface solid                  = " << M_solidInterfaceFlag << std::endl;
    output << "Interface structure              = " << M_structureInterfaceFlag << std::endl;
    output << "Interface harmonic               = " << M_harmonicInterfaceFlag << std::endl;
    output << "Interface tolerance              = " << M_interfaceTolerance << std::endl;

    output << "Robin-Neumann coupling           = " << M_RobinNeumannCoupling << std::endl;
    output << "Robin-Neumann fluid coefficient  = " << M_RobinNeumannFluidCoefficient << std::endl;
    output << "Robin-Neumann solid coefficient  = " << M_RobinNeumannSolidCoefficient << std::endl;
}

}