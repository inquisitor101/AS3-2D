


template<class InputPrecision>
void NImportFile::NAS3BinaryFile::ReadElementCoordinates
(
 const NImportFile::NAS3BinaryFile::CBinaryFile &file,
 CSinglezoneGeometry                            &grid
)
 /*
  * Read directly into element storage when precisions match.
  * Otherwise convert through one reusable element-sized buffer.
  */
{
  const std::size_t nDim    = 2;
  const std::size_t nValues = nDim * std::size_t(grid.GetnNodeGrid2D());
  
  as3vector1d<InputPrecision> buffer;
  if constexpr( !std::is_same_v<InputPrecision, as3double> ) buffer.resize(nValues);
  
  // Linear order matches the writer's j-outer, i-inner element loops.
  for( std::size_t iElem = 0; iElem < grid.GetnElem(); iElem++ )
  {
    auto& coordinates = grid.GetElementGeometry(iElem)->GetCoordSolDOFs();
  
    if constexpr( std::is_same_v<InputPrecision, as3double> )
    {
      if( std::fread(coordinates.data(), sizeof(InputPrecision),
                     nValues, file.mHandle.get()) != nValues )
      {
        ERROR("Could not read element coordinates: " + file.mFilename);
      }
  
      if( file.mSwap )
      {
        NInputUtility::SwapBytes( coordinates.data(), sizeof(InputPrecision), nValues );
      }
    }
    else
    {
      if( std::fread(buffer.data(), sizeof(InputPrecision),
                     nValues, file.mHandle.get()) != nValues )
      {
        ERROR("Could not read element coordinates: " + file.mFilename);
      }
  
      if( file.mSwap )
      {
        NInputUtility::SwapBytes( buffer.data(), sizeof(InputPrecision), nValues );
      }
  
      for( std::size_t i = 0; i < nValues; i++ )
      {
        coordinates[i] = static_cast<as3double>(buffer[i]);
      }
    }
  }
}
