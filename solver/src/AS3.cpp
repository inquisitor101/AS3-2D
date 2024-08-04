#include "AS3.hpp"



int main(int argc, char **argv)
{

  // Check if correct usage is executed.
  std::ostringstream message;
  if (argc != 2)
	{
    // Create instructions on how to use the program.
		message << "Usage: \n" << argv[0]
            << "\n <input config file>";
    ERROR(message.str());
  }

	// Instantiate a driver class and provide it with the configuration file.
	std::unique_ptr<CDriver> driver_container = std::make_unique<CDriver>( argv[1] );

	// Initialize the entire data for the simulation.
	driver_container->InitializeData();

	// Exit happily.
	return 0;
}

