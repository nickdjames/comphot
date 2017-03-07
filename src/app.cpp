#include <string>
#include <iostream>

#include <boost/program_options.hpp>

extern "C"
{
#include "comphot.h"
}

namespace bpo = boost::program_options;
namespace
{
// Program option data parsed from command line
struct Options
{
    std::string offsetimage;
    std::string fixedimage;
    std::string flatimage;
    int cenx = 0;
    int ceny = 0;
    int border = 5; // ignore pixels this close to edge for estimates
    float apradius = -1; // Default is to use auto
};

// Valid options must contain an offset and fixed image
bool isValid( const Options& options )
{
    return !options.offsetimage.empty() &&
        !options.fixedimage.empty();
}

std::ostream& printUsage(
    std::ostream& out,
    const char* name,
    const bpo::options_description& options )
{
    out << "Usage: " << name << " [options] offsetimage fixedimage cenx ceny:\n"
        << "\n  offsetimage/fixedimage Offset and fixed images (FIT)"
        << "\n  cenx/centy             Photocentre (pixels)"
        << "\n\nOptions:\n"
        << options << '\n';

    return out;
}

// Parse options from command line with the help of boost
// to handle Windows/Linux and other platforms
Options parseOptions( int argc, char** argv )
{
    Options options;

    // Required options that are expected to be given on
    // the command line without any option name
    bpo::options_description required;
    required.add_options()
    (
        "offsetimage",
        bpo::value< std::string >( &options.offsetimage )->required()
    )
    (
        "fixedimage",
        bpo::value< std::string >( &options.fixedimage )->required()
    )
    (
        "cenx",
        bpo::value< int >( &options.cenx )->required()
    )
    (
        "ceny",
        bpo::value< int >( &options.ceny )->required()
    );

    // Optional options that can be specified on the command line
    // using their given option names
    bpo::options_description optional;
    optional.add_options()
    (
        "apmax,r",
        bpo::value< float >( &options.apradius ),
        "Optional photometric aperture radius (arcsec)"
    )
    (
        "border,b",
        bpo::value< int >( &options.border ),
        "Ignore border pixels (default 5)"
    )
    (
        "flatimage,f",
        bpo::value< std::string >( &options.flatimage ),
        "Optional flat normalization image (FIT)"
    )
    (
        "help,h",
        "Usage info"
    );

    bpo::options_description all;
    all.add( required ).add( optional );

    // Set up the options that we expect without any option names
    bpo::positional_options_description pos;
    pos.add( "offsetimage", 1 );
    pos.add( "fixedimage", 1 );
    pos.add( "cenx", 1 );
    pos.add( "ceny", 1 );

    bpo::variables_map vm;

    try
    {
        bpo::store(
            bpo::command_line_parser( argc, argv )
                .options( all )
                .positional( pos )
                .run(), vm );

        if ( vm.count( "help" ) )
        {
            printUsage( std::cout, argv[ 0 ], optional );
        }
        else
        {
            // Check validity of all given options
            bpo::notify( vm );
        }
    }
    catch ( const bpo::error& e )
    {
        std::cout << "Error: " << e.what() << '\n';
        printUsage( std::cout, argv[ 0 ], optional );
    }

    return options;
}

// Convert option data to comphot internal format
ComphotConfig createComphotConfig( const Options& options )
{
    ComphotConfig config;
    config.offsetimage = options.offsetimage.c_str();
    config.fixedimage = options.fixedimage.c_str();
    config.cenx = options.cenx;
    config.ceny = options.ceny;
    config.border = options.border;
    config.apradius = options.apradius;

    if ( !options.flatimage.empty() )
    {
        config.flatimage = options.flatimage.c_str();
    }
    else
    {
        config.flatimage = nullptr;
    }

    return config;
}
} // namespace

int main( int argc, char* argv[] )
{
    Options options = parseOptions( argc, argv );

    if ( isValid( options ) )
    {
        ComphotConfig config = createComphotConfig( options );
        // call comphot to process
        process( &config );

        return 0;
    }

    return 1;
}

